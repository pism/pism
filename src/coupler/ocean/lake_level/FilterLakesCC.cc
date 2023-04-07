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
  {
    IceModelVec::AccessList list{ &m_mask_run, &lake_level };
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      //Set sink, where pism_mask is ocean or at margin of computational domain
      if (lake_level(i, j) != m_fill_value) {
        m_mask_run(i, j) = 2;
      } else {
        m_mask_run(i, j) = 0;
      }
    }
    m_mask_run.update_ghosts();
  }

  connected_components::set_validity_mask(n_filter,
                                          lake_level,
                                          [this](double value) { return value != m_fill_value; },
                                          m_mask_validity);

  VecList lists;
  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(lists, max_items);

  int run_number = 1;

  compute_runs(run_number, lists, max_items);

  {
    connected_components::set_labels(run_number, lists, lake_level);

    const auto &valid = lists.find("valid")->second;
    auto condition    = [valid](double label) { return not (valid[(int)label] > 0); };
    connected_components::replace_values(lake_level, condition, m_fill_value);
  }
}

bool FilterLakesCC::ForegroundCond(int i, int j) const {
  return (m_mask_run(i, j) > 1);
}

} // end of namespace pism
