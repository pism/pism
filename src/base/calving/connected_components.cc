/* Copyright (C) 2013, 2014 PISM Authors
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
#include <cmath>

void run_union(std::vector<unsigned int> &parent, unsigned int run1, unsigned int run2) {
  if (parent[run1] == run2 || parent[run2] == run1)
    return;

  while (parent[run1] != 0)
    run1 = parent[run1];

  while (parent[run2] != 0)
    run2 = parent[run2];

  if (run1 > run2)
    parent[run1] = run2;
  else if (run1 < run2)
    parent[run2] = run1;
  else
    return;

}

void cc(double *image, unsigned int n_rows, unsigned int n_cols, bool identify_icebergs, double mask_grounded) {
  unsigned int max_runs = 2*n_rows;
  const double eps = 1e-6;

  std::vector<unsigned int> parents(max_runs), lengths(max_runs),
    rows(max_runs), columns(max_runs), mask(max_runs);

  unsigned int run_number = 0, r, c, parent;

  // First scan
  for (r = 0; r < n_rows; ++r) {
    for (c = 0; c < n_cols; ++c) {

      if (image[r*n_cols + c] > 0.0) {
        // looking at a foreground pixel
        if (c > 0 && image[r*n_cols + (c-1)] > 0.0)
          // one to the left is also a foreground pixel; continue the run
          lengths[run_number] += 1;
        else {
          // one to the left is a background pixel (or this is column 0); start a new run

          // set the run just above as a parent (if present)
          if (r > 0 && image[(r-1)*n_cols + c] > 0.0)
            parent = (unsigned int)image[(r-1)*n_cols + c];
          else
            parent = 0;

          run_number += 1;

          // allocate more storage (if needed)
          if (run_number == max_runs) {
            max_runs += n_rows;
            parents.resize(max_runs);
            lengths.resize(max_runs);
            rows.resize(max_runs);
            columns.resize(max_runs);
            mask.resize(max_runs);
          }

          // Record this run
          rows[run_number]    = r;
          columns[run_number] = c;
          lengths[run_number] = 1;
          parents[run_number] = parent;
          mask[run_number]    = 0;
        }

        if (r > 0 && image[(r-1)*n_cols + c] > 0.0)
          run_union(parents, (unsigned int)image[(r-1)*n_cols + c], run_number);

        if (mask[run_number] == 0 && fabs(image[r*n_cols + c] - mask_grounded) < eps)
          mask[run_number] = 1;

        image[r*n_cols + c] = run_number;
      }
    }
  }

  // Assign labels to runs.
  // This uses the fact that children always follow parents,
  // so we can do just one sweep: by the time we get to a node (run),
  // its parent already has a final label.
  //
  // We use "parents" to store labels here, because once a run's label is computed
  // we don't need to know its parent run any more.

  unsigned int label = 0;
  std::vector<unsigned int> grounded(run_number + 1);
  for (r = 0; r <= run_number; ++r) {
    if (parents[r] == 0) {
      parents[r] = label;
      label += 1;
    } else {
      parents[r] = parents[parents[r]];
    }

    // remember current blob (parents[r]) as "grounded" if the current run is
    // "grounded"
    if (mask[r] == 1)
      grounded[parents[r]] = 1;
  }

  // Second scan (re-label)
  if (identify_icebergs) {
    for (r = 0; r <= run_number; ++r) {
      for (c = 0; c < lengths[r]; ++c)
        image[rows[r]*n_cols + columns[r] + c] = 1 - grounded[parents[r]];
    }
  } else {
    for (r = 0; r <= run_number; ++r) {
      for (c = 0; c < lengths[r]; ++c)
        image[rows[r]*n_cols + columns[r] + c] = parents[r];
    }
  }

  // Done!
}
