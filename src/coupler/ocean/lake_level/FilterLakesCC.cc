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

#include "FilterLakesCC.hh"

namespace pism {

FilterLakesCC::FilterLakesCC(IceGrid::ConstPtr g, double fill_value)
  : ValidCC<ConnectedComponents>(g),
    m_fill_value(fill_value) {
  //empty
}

void FilterLakesCC::filter_map(int n_filter, IceModelVec2S &lake_level) {
  // m_mask_run will be used in compute_runs() via ForegroundCond()
  prepare_mask(lake_level, m_mask_run);

  set_mask_validity(n_filter, lake_level, m_mask_validity);

  VecList lists;
  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(lists, max_items);

  int run_number = 1;

  compute_runs(run_number, lists, max_items);

  labelMap(run_number, lists, lake_level);
}

bool FilterLakesCC::ForegroundCond(int i, int j) const {
  int mask = m_mask_run.as_int(i, j);

  return (mask > 1);
}

void FilterLakesCC::labelMap(int run_number, const VecList &lists, IceModelVec2S &result) const {
  IceModelVec::AccessList list{&result};

  const auto
    &i_vec      = lists.find("i")->second,
    &j_vec      = lists.find("j")->second,
    &len_vec    = lists.find("lengths")->second,
    &parents    = lists.find("parents")->second,
    &valid_list = lists.find("valid")->second;

  for(int k = 0; k <= run_number; ++k) {
    const int label = trackParentRun(k, parents);
    const bool valid = (valid_list[label] > 0);
    if (not valid) {
      const int j = j_vec[k];
      for(int n = 0; n < len_vec[k]; ++n) {
        const int i = i_vec[k] + n;
        result(i, j) = m_fill_value;
      }
    }
  }
}

void FilterLakesCC::prepare_mask(const IceModelVec2S &lake_level, IceModelVec2Int &result) {
  IceModelVec::AccessList list{ &result, &lake_level};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    //Set sink, where pism_mask is ocean or at margin of computational domain
    if (isLake(lake_level(i, j))) {
      result(i, j) = 2;
    } else {
      result(i, j) = 0;
    }
  }
  result.update_ghosts();
}

void FilterLakesCC::set_mask_validity(int threshold, const IceModelVec2S &lake_level,
                                      IceModelVec2Int &result) {

  IceModelVec2S ll_tmp(m_grid, "temp_lake_level", WITH_GHOSTS, 1);
  ll_tmp.copy_from(lake_level);

  IceModelVec::AccessList list{ &result, &ll_tmp };
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    int n_neighbors = 0;
    if (ll_tmp(i, j) != m_fill_value) {
      auto level = ll_tmp.star(i, j);
      for (auto direction : { North, East, South, West }) {
        if (level[direction] != m_fill_value) {
          ++n_neighbors;
        }
      }
    }

    // Set cell valid if the number of neighbors exceeds threshold
    if (n_neighbors >= threshold) {
      result(i, j) = 1;
    } else {
      result(i, j) = 0;
    }
  } // end of the loop over grid points

  result.update_ghosts();
}

} // end of namespace pism