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

#include <array>
#include <vector>

/*!
 * This file contains the implementation of the *serial* connected component labeling
 * algorithm using run-length encoding and the disjoint-set data structure
 * (https://en.wikipedia.org/wiki/Disjoint-set_data_structure).
 *
 */

namespace pism {
namespace connected_components {
namespace details {

inline int resolve_label(const std::vector<int> &labels, int provisional_label) {
  int final_label = provisional_label;

  while (final_label != labels[final_label]) {
    final_label = labels[final_label];
  }

  return final_label;
}

struct Run {
  int row, col, length, label;
};

/*!
 * Helper class wrapping a C-style 2D array using the row-major storage order. It is used
 * to implement the interface required by connected_components::details::label().
 */
template<typename T>
class CArray {
public:
  CArray(T *array, int nrows, int ncols) : m_array(array), m_nrows(nrows), m_ncols(ncols) {
    // empty
  }

  T &operator()(int r, int c) {
    return m_array[c + m_ncols * r];
  }

  T &operator()(int r, int c) const {
    return m_array[c + m_ncols * r];
  }

  std::array<int, 2> shape() const {
    return { m_nrows, m_ncols };
  }

private:
  T *m_array;
  int m_nrows;
  int m_ncols;
};

/*!
 * Helper class wrapping `pism::array::Scalar` to implement the interface required to use
 * connected_components::details::label().
 */
class PISMArray {
public:
  PISMArray(pism::array::Scalar &array)
    : m_array(array) {

    auto grid = array.grid();

    m_nrows = grid->ym();
    m_ncols = grid->xm();

    m_c0 = grid->xs();
    m_r0 = grid->ys();

    m_array.begin_access();
  }

  ~PISMArray() {
    m_array.end_access();
  }

  inline double &operator()(int r, int c) {
    return m_array(m_c0 + c, m_r0 + r);
  }

  inline double &operator()(int r, int c) const {
    return m_array(m_c0 + c, m_r0 + r);
  }

  std::array<int,2> shape() const {
    return {m_nrows, m_ncols};
  }

private:
  array::Scalar &m_array;
  int m_nrows;
  int m_ncols;
  int m_r0;
  int m_c0;
};

//! Adds "foregrdound" and "attached" concepts to an "array".
template <typename ARRAY>
class Mask {
public:
  Mask(ARRAY &data, int reachable_value = -1) : m_data(data), m_reachable(reachable_value) {
    // empty
  }

  std::array<int,2> shape() const {
    return m_data.shape();
  }

  inline bool is_foreground(int i, int j) const {
    return m_data(i, j) > 0;
  }

  inline bool is_attached(int i, int j) const {
    return static_cast<int>(m_data(i, j)) == m_reachable;
  }

private:
  ARRAY &m_data;
  int m_reachable;
};

/*!
 * Label connected components using a *serial* algorithm. Designed to be used on its own
 * or as a part of the "parallel" implementation.
 *
 * Uses `input` to generate the mask of patches "on the fly", saves results to `output`.
 *
 * This is implemented as a template function to allow using different "2D array" types.
 *
 * The "mask" type has to implement methods
 *
 * - `bool is_foreground(row, col)` indicating if a grid cell is a "foreground" cell (part
 *   of a "patch" or a "component") or "background".
 *
 * - `bool is_attached(row, col)` indicating if a grid cell is "attached". Patches that do
 *   not contain any "attached" cells are considered "isolated". This method is called
 *   only if `identify_isolated_patches` is `true`.
 *
 * The "array" type has to implement `number_t& operator(row, col)` for some numeric type
 * `number_t` (int, double, etc). The argument `output` is "write only".
 *
 * If `assign_final_labels` is true, then set labels to "final" values. This is
 * appropriate when this function is used on its own (serial algorithm).
 *
 * Specifically:
 *
 * - If `assign_final_labels` is true and `identify_isolated_patches` is true, then
 *   isolated patches are marked with `1` and the rest are set to `0`.
 *
 * - If `assign_final_labels` is true and `identify_isolated_patches` is false, then
     patches are labeled with consecutive numbers starting from `min_label`.
 *
 * - If `assign_final_labels` is false then patches "attached" to cells where
 *   `is_attached(i, j)` is true are marked with the smallest odd number that is greater
 *   than or equal to `min_label`. All other patches get consecutive *even* labels. This
 *   labeling scheme is used by the parallel version of this code.
 *
 * Uses labels starting from `min_label`. In the parallel implementation each MPI rank
 * (each sub-domain) has to use *unique* labels; `min_labels` is used to guarantee this.
 *
 * Note: only foreground cells of `output` are updated.
 */
template <typename mask, typename array>
bool label(const mask &input, bool identify_isolated_patched, int min_label,
           bool assign_final_labels, array &output) {

  // iteration limits:
  const int col_min = 0;
  const int row_min = 0;
  auto shape = input.shape();
  const int nrows = shape[0];
  const int ncols = shape[1];
  // Labels: the "attached" label should be the smallest one.
  const int background            = 0;
  const int attached_label        = 1;
  const int min_provisional_label = 2;

  // This array encodes equivalences between labels: label k is
  // equivalent to labels[k], which is equivalent to
  // labels[labels[k]], etc. The smallest label from a set of
  // equivalent ones is used as a "representative" label for the set.
  // By design labels[k] <= k.
  std::vector<int> labels = {background, attached_label};
  int provisional_label = min_provisional_label;

  std::vector<Run> runs;

  // storage for labels assigned to pixels in the row above the current one
  std::vector<int> row_above(ncols, background);

  for (int r = row_min; r < nrows; ++r) {

    int c = col_min;
    while (c < ncols) {
      if (input.is_foreground(r, c)) {
        // we're looking at a foreground pixel that must be the
        // beginning of a run
        int L = provisional_label;
        int c_start = c;

        runs.push_back({ r /* row */, c_start /* column */, 0 /* length */, L /* label */});
        Run &current_run = runs.back();

        // Iterate over all the pixels in this run
        while (c < ncols and input.is_foreground(r, c)) {

          if (identify_isolated_patched and input.is_attached(r, c)) {
            // looking at a pixel attached to a "grounded" area
            if (L != provisional_label) {
              labels[L] = attached_label;
            }
            L = attached_label;
          }

          int T = row_above[c];
          c += 1;

          if (T != background) {
            // foreground pixel in the row above
            if (T < L) {
              if (L != provisional_label) {
                labels[L] = T;
              }
              L = T;
            } else if (T > L) {
              labels[T] = L;
            }
          }
        } // end of the loop over pixels in a run

        if (L == provisional_label) {
          // Failed to assign a label by looking at the row above: we
          // need to add this label to the array encoding
          // equivalences.
          labels.push_back(L);
          provisional_label += 1;
        }

        // Done with a run: record the length and the label.
        current_run.length = c - c_start;
        current_run.label = L;

        // Record pixel labels in the row above.
        for (int n = 0; n < current_run.length; ++n) {
          row_above[c_start + n] = L;
        }
      } else {
        // background pixel
        row_above[c] = background;
        c += 1;
      }
    } // end of the loop over columns
  } // end of the loop over rows

  auto N_labels = static_cast<int>(labels.size());

  // Flatten the table of equivalences
  {
    for (int k = 0; k < N_labels; ++k) {
      labels[k] = resolve_label(labels, k);
    }
  }

  // Rename labels.
  //
  if (assign_final_labels) {
    // this block is executed by the serial version
    if (identify_isolated_patched) {
      // isolated patches are labeled with "1", the rest with "0":
      for (int k = 0; k < N_labels; ++k) {
        labels[k] = (labels[k] == attached_label) ? 0 : 1;
      }
    } else {
      // Use consecutive numbers starting from `min_label`.
      int L = min_label;
      for (int k = min_provisional_label; k < N_labels; ++k) {
        if (labels[k] == k) {
          labels[k] = L;
          L += 1;
        } else {
          labels[k] = labels[labels[k]];
        }
      }
    }
  } else {
    // This block is executed by the parallel version (we need to use unique labels to
    // build the graph describing connections between patches on different sub-domains).
    //
    //
    // Use the smallest allowed odd number to indicate "attached" cells:
    const int attached = min_label + (1 - (min_label % 2));

    // Use consecutive even numbers for the rest:
    int L = attached + 1;
    for (int k = min_provisional_label; k < N_labels; ++k) {
      if (labels[k] == attached_label) {
        labels[k] = attached;
      } else if (labels[k] == k) {
        labels[k] = L;
        L += 2;
      } else {
        labels[k] = labels[labels[k]];
      }
    }
  }

  // Second scan: assign labels. Note that this does not touch background pixels (we
  // iterate over runs).
  for (int k = 0; k < (int)runs.size(); ++k) {
    auto &r = runs[k];
    int L = labels[r.label];
    for (int n = 0; n < r.length; ++n) {
      output(r.row, r.col + n) = L;
    }
  }

  return (not runs.empty());
}

} // end of namespace details
} // end of namespace connected_components
} // end of namespace pism
