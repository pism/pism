/* Copyright (C) 2013, 2014, 2016, 2021 PISM Authors
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

#include <vector>
#include <algorithm>            // min, max, abs

int resolve_label(const std::vector<int> &labels, int provisional_label) {
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
 * Label connected components. Modifies `image` in place.
 */
void label_connected_components(double *image, int nrows, int ncols, bool identify_icebergs, int mask_grounded, int first_label) {

  // storage for labels assigned to pixels in the row above the current one
  std::vector<int> row_above(ncols, 0);

  // Labels: the "grounded" label should be the smallest one.
  int grounded_label = 1;
  int provisional_label = 2;
  // This array encodes equivalences between labels: label k is
  // equivalent to labels[k], which is equivalent to
  // labels[labels[k]], etc. The smallest label from a set of
  // equivalent ones is used as a "representative" label for the set.
  // By design labels[k] <= k.
  std::vector<int> labels = {0, grounded_label};

  std::vector<Run> runs;

  for (int r = 0; r < nrows; ++r) {
    const double *row = &image[r * ncols];

    int c = 0;
    while (c < ncols) {
      if (row[c] > 0) {
        // we're looking at a foreground pixel that must be the
        // beginning of a run
        int L = provisional_label;
        int c_start = c;

        runs.push_back({r, c_start, 0, L});
        Run &current_run = runs.back();

        // Iterate over all the pixels in this run
        while (c < ncols and row[c] > 0) {

          if (identify_icebergs and (int)row[c] == mask_grounded) {
            // looking at a "grounded" pixel
            if (L != provisional_label) {
              labels[L] = grounded_label;
            }
            L = grounded_label;
          }

          int T = row_above[c];
          c += 1;

          if (T > 0) {
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
        row_above[c] = 0;
        c += 1;
      }
    } // end of the loop over columns
  } // end of the loop over rows

  int N_labels = labels.size();

  // Flatten the table of equivalences
  {
    for (int k = 0; k < N_labels; ++k) {
      labels[k] = resolve_label(labels, k);
    }
  }

  // Rename labels to make them consecutive
  {
    int L = first_label;
    for (int k = 1; k < N_labels; ++k) {
      if (labels[k] == k) {
        labels[k] = L;
        L += 1;
      } else {
        labels[k] = labels[labels[k]];
      }
    }
  }

  if (identify_icebergs) {
    // Blobs connected to grounded areas have the label "1", icebergs
    // have labels 2 and greater.
    for (int k = 1; k < N_labels; ++k) {
      labels[k] = labels[k] > 1;
    }
  } else {
    // Here we subtract 1 because provisional labels start at 2 (1 is
    // a special label used to identify blobs connected to "grounded"
    // areas).
    for (int k = 1; k < N_labels; ++k) {
      labels[k] -= 1;
    }
  }

  // Second scan: assign labels
  for (int k = 0; k < (int)runs.size(); ++k) {
    auto &r = runs[k];
    int run_start = r.row * ncols + r.col;
    for (int n = 0; n < r.length; ++n) {
      image[run_start + n] = labels[r.label];
    }
  }
}
