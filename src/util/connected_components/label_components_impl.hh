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

#include "pism/util/array/Scalar.hh"
#include "connected_components_impl.hh"
#include <map>

namespace pism {

namespace details {

int first_label(const Grid &grid);

//! Assign new labels to elements of `mask`. Does not touch background grid cells.
void relabel(array::Scalar &mask, const std::map<int, int> &labels);

std::map<int, int> final_labels(array::Scalar1 &input, bool subdomain_is_not_empty,
                                bool mark_isolated_patches);

} // end of namespace details

/*!
 * Labels connected components in parallel, using a custom definition of "foreground" and
 * "attached" grid cells.
 *
 * Note that `connected_components::details::label()` below is a template function: this
 * is why we have to include the implementation here.
 *
 * The type `T` has to implement `is_foreground(row, col)` and `is_attached(row, col)`.
 *
 * Here's the idea:
 *
 * 1. Identify connected components in each sub-domain.
 *
 * 2. "Update ghosts" of `output`, then iterate over sub-domain edges to identify
 *    connections between patches in sub-domains that make up connected components
 *    spanning multiple sub-domains. This defines a graph: each patch on a sub-domain is a
 *    node, two nodes are connected by an edge if and only if they "touch".
 *
 * 3. Gather graph description on all sub-domains that have at least one patch.
 *
 * 4. Use breadth-first search to traverse this graph and compute final labels.
 *
 * 5. Apply final labels.
 *
 * This method communicates ghosts (once), number of graph edges per sub-domain (once) and
 * then all graph edges (once, only to sub-domains that have at least one patch).
 *
 * Graph traversal is done "redundantly", i.e. each participating sub-domain traverses the
 * whole graph even if it contains one isolated node. This is needed to ensure that
 * resulting labels use consecutive numbers. (Consecutive labels are useful for indexing
 * elsewhere in the code).
 *
 * We /could/ gather graph information on one MPI processor, traverse the graph to compute
 * final labels, then scatter final labels. It is not clear if this would be better, though.
 *
 * The current implementation uses
 *
 * - MPI_Allgather() to gather number of graph edges per subdomain (send 4 bytes, receive
     4 bytes per subdomain).
 *
 * - MPI_Allgatherv() to gather edges to all participating subdomains (send 8 bytes per
     local edge, receive 8-16 bytes per edge).
 *
 * The alternative implementation would use
 *
 * - MPI_Gather() to gather number of graph edges per sub-domain to *one* of sub-domains.
 *   (each sub-domain sends 4 bytes, one sub-domain receives 4 bytes per sub-domain).
 *
 * - MPI_Gatherv() to gather edges from all participating sub-domains. (All sub-domains send
 *   8 bytes per local edge, one sub-domain receives 8 bytes per edge in the whole graph.)
 *
 * - MPI_Bcast() to scatter the mapping from old labels to new labels (8 bytes per local
 *   sub-domain).
 *
 * TO DO: run benchmarks!
 */
template <typename T>
void label_components_impl(const T &input, bool mark_isolated_patches, array::Scalar1 &output) {

  // 1. Label patches owned by individual sub-domains (independently and in parallel):
  bool non_empty_subdomain = false;
  {
    connected_components::details::PISMArray out(output);

    using namespace connected_components::details;
    bool assign_final_labels = false;
    non_empty_subdomain = label(input, mark_isolated_patches, details::first_label(*output.grid()),
                                assign_final_labels, out);
  }

  // 2. Resolve labels:
  auto labels = details::final_labels(output, non_empty_subdomain, mark_isolated_patches);

  // 3. Apply final labels:
  details::relabel(output, labels);
}

} // end of namespace pism
