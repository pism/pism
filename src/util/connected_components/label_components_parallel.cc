/* Copyright (C) 2023 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "label_components.hh"
#include "label_components_impl.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/Scalar.hh"

#include <cmath>  // pow, ceil, log10
#include <map>
#include <mpi.h>
#include <queue>

namespace pism {

namespace details {

//! West boundary of a sub-domain.
static bool west_boundary(const Grid &grid, int i) {
  int i_first = grid.xs();
  return (i == i_first) and (i != 0);
};

//! East boundary of a sub-domain.
static bool east_boundary(const Grid &grid, int i) {
  int i_last = grid.xs() + grid.xm() - 1;
  return (i == i_last) and (i != (int)grid.Mx() - 1);
};

//! North boundary of a sub-domain.
static bool north_boundary(const Grid &grid, int j) {
  int j_last = grid.ys() + grid.ym() - 1;
  return (j == j_last) and (j != (int)grid.My() - 1);
};

//! South boundary of a sub-domain.
static bool south_boundary(const Grid &grid, int j) {
  int j_first = grid.ys();
  return (j == j_first) and (j != 0);
};

/*!
 * Inspect sub-domain edges to detect connections between patches owned by individual
 * sub-domains.
 */
static std::map<int, std::set<int> > detect_connections(array::Scalar1 &mask) {
  auto grid = mask.grid();
  std::map<int, std::set<int> > graph;

  mask.update_ghosts();

  array::AccessScope list{ &mask };

  auto maybe_add_edge = [&](int a, int b) {
    if (a > 0 and b > 0) {
      graph[std::min(a, b)].insert(std::max(a, b));
    }
  };

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int M = mask.as_int(i, j);

    // add an edge (M, M) to ensure that the graph contains information about patches that
    // don't touch sub-domain edges. This is needed to be able to re-label them properly
    // (i.e. to make labels consecutive).
    //
    // Note that this is incompatible with the optimization mentioned in
    // label_components_parallel()
    maybe_add_edge(M, M);

    if (north_boundary(*grid, j)) {
      maybe_add_edge(M, mask.as_int(i, j + 1));
    }

    if (east_boundary(*grid, i)) {
      maybe_add_edge(M, mask.as_int(i + 1, j));
    }

    if (south_boundary(*grid, j)) {
      maybe_add_edge(M, mask.as_int(i, j - 1));
    }

    if (west_boundary(*grid, i)) {
      maybe_add_edge(M, mask.as_int(i - 1, j));
    }
  } // end of the loop over grid points

  return graph;
}

/*!
 * Pack information about connections (graph edges) into the form needed for communication.
 */
static std::vector<int> pack_data(const std::map<int, std::set<int> > &graph) {
  // FIXME (maybe): we could reduce the amount of data we have to send by being more
  // clever about packing data.
  std::vector<int> result{};
  for (const auto &e : graph) {
    int a = e.first;
    for (auto b : e.second) {
      result.push_back(a);
      result.push_back(b);
    }
  }
  return result;
}

/*!
 * Unpack data gathered from all the ranks that contain connected components.
 */
static std::map<int, std::set<int> > unpack_data(const std::vector<int> &buffer) {
  std::map<int, std::set<int> > result;
  for (size_t n = 0; n < buffer.size() / 2; ++n) {
    int a = buffer[2 * n + 0];
    int b = buffer[2 * n + 1];

    // insert *both* (a,b) and (b,a) because we want to be able to easily get a list
    // of nodes adjacent to a given node
    result[a].insert(b);
    result[b].insert(a);
  }

  return result;
}

/*!
 * Assemble complete graph info by gathering edges from all MPI ranks that have at least
 * one patch.
 */
static std::map<int, std::set<int> > assemble_graph(MPI_Comm comm,
                                                    const std::map<int, std::set<int> > &subgraph) {
  int size = 0;
  MPI_Comm_size(comm, &size);
  int rank = 0;
  MPI_Comm_rank(comm, &rank);

  auto message = pack_data(subgraph);

  // gather message lengths to prepare the buffer:
  std::vector<int> message_length(size);
  message_length[rank] = (int)message.size();
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, message_length.data(), 1, MPI_INT, comm);

  // compute offsets
  std::vector<int> offsets(size);
  offsets[0] = 0;
  for (int k = 1; k < size; ++k) {
    offsets[k] = offsets[k - 1] + message_length[k - 1];
  }

  // allocate the buffer
  int buffer_size = offsets.back() + message_length.back();
  std::vector<int> buffer(buffer_size);

  // Gather adjacency info:
  MPI_Allgatherv(message.data(), (int)message.size(), MPI_INT, buffer.data(), message_length.data(),
                 offsets.data(), MPI_INT, comm);

  // unpack messages:
  return unpack_data(buffer);
}

/*!
 * Traverse the `graph` using a version of the "breadth first search" algorithm and find
 * final labels for all of its nodes.
 *
 * If `mark_isolated_patches` is true, isolated patches get labeled with `1`, the rest
 * with `0`.
 *
 * If `mark_isolated_patches` is false, each connected component in the graph gets a
 * unique label starting from `1` and using consecutive integers.
 */
static std::map<int, int> resolve_labels(const std::map<int, std::set<int> > &graph,
                                         bool mark_isolated_patches) {

  // visited[X] == true if we already visited the node with the ID equal to X
  std::map<int, bool> visited;
  // map from node IDs to labels
  std::map<int, int> result;

  // label used for "attached" cells
  const int attached = 1;

  // Loop over all edges in the graph and start the breadth-first search starting from the
  // first node of each edge. This is sufficient since if (a,b) is in the graph then (b,a)
  // is also in the graph (by construction).
  for (const auto &edge : graph) {
    int u = edge.first;

    if (visited[u]) {
      continue;
    }

    // select a provisional label:
    int label = u;

    std::set<int> equivalent{};

    std::queue<int> queue;
    queue.push(u);
    while (not queue.empty()) {
      int v = queue.front();
      queue.pop();

      if (visited[v]) {
        continue;
      }

      // Odd labels are used to mark "attached" cells; replace this label with "1". Note
      // that this is the smallest label we use, so `min(label, v)` will preserve it.
      label = (v % 2 == 1) ? attached : std::min(label, v);

      equivalent.insert(v);
      visited[v] = true;

      auto it = graph.find(v);
      if (it != graph.end()) {
        for (int w : it->second) {
          queue.push(w);
        }
      }
    }

    // record the label for all the "equivalent" nodes (i.e. all the nodes in a particular
    // connected component)
    for (int v : equivalent) {
      result[v] = label;
    }
  }

  if (mark_isolated_patches) {
    // mark "attached" cells with `0` and the rest with `1`:
    for (auto &p : result) {
      p.second = (p.second == attached) ? 0 : 1;
    }
  } else {
    // make labels consecutive:
    int current_label = 1;
    std::map<int, int> new_labels;
    for (auto &p : result) {
      if (new_labels.find(p.second) == new_labels.end()) {
        new_labels[p.second] = current_label++;
      }
      p.second = new_labels[p.second];
    }
  }

  return result;
}

//! Assign new labels to elements of `mask`. Does not touch background grid cells.
/*!
 * This function assumes that `mask` was labeled using
 * `connected_components::details::label()`
 */
void relabel(array::Scalar &mask, const std::map<int, int> &labels) {

  if (labels.empty()) {
    return;
  }

  array::AccessScope list{ &mask };

  for (Points p(*mask.grid()); p; p.next()) {
    const int i = p.i(), j = p.j();

    int old_label = mask.as_int(i, j);

    auto it = labels.find(old_label);
    if (it == labels.end()) {
      continue;
    }

    mask(i, j) = it->second;
  }
}

std::map<int, int> final_labels(array::Scalar1 &input, bool subdomain_is_not_empty,
                                bool mark_isolated_patches) {
  // Iterate over the grid to find connections between patches owned by individual
  // sub-domains. Updates the ghosts of `mask`.
  auto edges = details::detect_connections(input);

  // We create a sub-communicator containing *only* the ranks with at least one patch to
  // reduce the amount of data sent around by the rest of this algorithm. We may be able
  // to reduce this even further by excluding sub-domains with patches that are *entirely*
  // contained within a sub-domain. On the other hand, this would make it impossible to
  // ensure that patches are labeled with consecutive numbers starting from 1 (this
  // simplifies the code elsewhere).
  //
  // This optimization *is* an option if `mark_isolated_patches` is set since in that case
  // the only values used are `1` for isolated patches and `0` for the rest of the domain.
  std::map<int, int> labels;
  {
    MPI_Comm comm = MPI_COMM_NULL;
    MPI_Comm_split(input.grid()->com, subdomain_is_not_empty ? 1 : 0, 0, &comm);

    // the following block is executed by ranks that have at least one patch
    if (subdomain_is_not_empty) {
      auto graph = details::assemble_graph(comm, edges);

      labels = details::resolve_labels(graph, mark_isolated_patches);
    }

    MPI_Comm_free(&comm);
  }
  return labels;
}

//! Compute the first label a particular rank can use. This is supposed to ensure that all
//! labels are unique (i.e. different ranks do not use the same label).
int first_label(const Grid &grid) {

  // Here we set `first_label` to ensure that all ranks use *different* labels.
  // find the smallest N such that grid.xm()*grid.ym() < 10^N
  int exponent        = static_cast<int>(std::ceil(std::log10(grid.max_patch_size())));
  int labels_per_rank = static_cast<int>(std::pow(10, exponent));

  // FIXME: this is very unlikely, but I should still check for integer overflow.

  return grid.rank() * labels_per_rank + 1;
}

} // end of namespace details

namespace connected_components {

using Array = connected_components::details::PISMArray;
using Mask = connected_components::details::Mask<Array>;

void label(array::Scalar1 &mask) {

  Array input(mask);

  bool mark_isolated_patches = false;
  label_components_impl(Mask(input, -1), mark_isolated_patches, mask);
}

void label_isolated(array::Scalar1 &mask, int reachable) {

  Array input(mask);

  bool mark_isolated_patches = true;
  label_components_impl(Mask(input, reachable), mark_isolated_patches, mask);
}

} // end of namespace connected_components

} // end of namespace pism
