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

#include "IsolationCC.hh"

namespace pism {

IsolationCC::IsolationCC(IceGrid::ConstPtr g, const IceModelVec2S &thk,
                         const double thk_theshold)
  : SinkCC(g),
    m_thk_threshold(thk_theshold),
    m_thk(&thk) {
  {
    IceModelVec::AccessList list{ &m_mask_run };
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      //Set not isolated at margin
      m_mask_run(i, j) = grid_edge(*m_grid, i, j) ? 1 : 0;
    }
    m_mask_run.update_ghosts();
  }

  m_fields.push_back(m_thk);
}

void IsolationCC::find_isolated_spots(IceModelVec2Int &result) {
  VecList lists;
  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(lists, max_items);

  int run_number = 1;

  compute_runs(run_number, lists, max_items);

  connected_components::set_labels(run_number, lists, result);
  connected_components::replace_values(result,
                                       [](double label) { return label == 1; },
                                       1);
}

bool IsolationCC::ForegroundCond(int i, int j) const {
  double thk = (*m_thk)(i, j);
  int mask = m_mask_run.as_int(i, j);
  return ((thk < m_thk_threshold) or (mask > 0));
}

} // end of namespace pism
