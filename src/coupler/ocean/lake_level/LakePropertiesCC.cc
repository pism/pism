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

#include "LakePropertiesCC.hh"

namespace pism {

using connected_components::is_valid;

LakePropertiesCC::LakePropertiesCC(IceGrid::ConstPtr g, const double fill_value, const IceModelVec2S &target_level,
                                   const IceModelVec2S &lake_level)
  : ConnectedComponents(g), m_fill_value(fill_value), m_target_level(&target_level),
    m_current_level(&lake_level) {

  m_min_lakelevel.create(m_grid, "min_ll_mask", WITH_GHOSTS, 1);
  m_min_lakelevel.set(m_fill_value);

  m_max_lakelevel.create(m_grid, "max_ll_mask", WITH_GHOSTS, 1);
  m_max_lakelevel.set(m_fill_value);

  m_masks.push_back(&m_min_lakelevel);
  m_masks.push_back(&m_max_lakelevel);

  m_fields.push_back(m_target_level);
  m_fields.push_back(m_current_level);
}

void LakePropertiesCC::getLakeProperties(IceModelVec2S &min_level, IceModelVec2S &max_level) {
  VecList lists;
  unsigned int max_items = 2 * m_grid->ym();
  init_VecList(lists, max_items);

  int run_number = 1;

  compute_runs(run_number, lists, max_items);

  min_level.copy_from(m_min_lakelevel);
  max_level.copy_from(m_max_lakelevel);
}

void LakePropertiesCC::init_VecList(VecList &lists, const unsigned int length) {
  ConnectedComponents::init_VecList(lists, length);

  std::vector<double> min_ll_list(length), max_ll_list(length);
  lists["min_ll"]  = min_ll_list;
  lists["max_ll"]  = max_ll_list;

  for (unsigned int k = 0; k < 2; ++k) {
    lists["min_ll"][k]  = m_fill_value;
    lists["max_ll"][k]  = m_fill_value;
  }
}

void LakePropertiesCC::setRunMinLevel(double level, int run, VecList &lists) {
  if (run == 0) {
    return;
  }

  run = connected_components::trackParentRun(run, lists["parents"]);
  if (is_valid(level)) {
    if (is_valid(lists["min_ll"][run])) {
      level = std::min(level, lists["min_ll"][run]);
    }
    lists["min_ll"][run] = level;
  }
}

void LakePropertiesCC::setRunMaxLevel(double level, int run, VecList &lists) {
  if (run == 0) {
    return;
  }

  run = connected_components::trackParentRun(run, lists["parents"]);
  if (is_valid(level)) {
    if (is_valid(lists["max_ll"][run])) {
      level = std::max(level, lists["max_ll"][run]);
    }
    lists["max_ll"][run] = level;
  }
}


bool LakePropertiesCC::ForegroundCond(int i, int j) const {
  return is_valid((*m_target_level)(i, j));
}

void LakePropertiesCC::labelMask(int run_number, const VecList &lists) {

  connected_components::set_labels(run_number, lists, m_mask_run);

  IceModelVec::AccessList list{&m_mask_run, &m_min_lakelevel, &m_max_lakelevel};

  const auto
    &min_ll  = lists.find("min_ll")->second,
    &max_ll  = lists.find("max_ll")->second;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto label = static_cast<int>(m_mask_run(i, j));
    m_min_lakelevel(i, j) = min_ll[label];
    m_max_lakelevel(i, j) = max_ll[label];
  }
}

static double min(double a, double b) {
  return (is_valid(a) and ((a < b) or not is_valid(b))) ? a : b;
}

static double max(double a, double b) {
  return (is_valid(a) and ((a > b) or not is_valid(b))) ? a : b;
}

void LakePropertiesCC::treatInnerMargin(int i, int j,
                                        bool isNorth, bool isEast, bool isSouth, bool isWest,
                                        VecList &lists, bool &changed) {
  ConnectedComponents::treatInnerMargin(i, j, isNorth, isEast, isSouth, isWest, lists, changed);

  int run = m_mask_run.as_int(i, j);
  if (run <= 0) {
    return;
  }

  StarStencil<bool> flags;
  flags.e = isEast;
  flags.w = isWest;
  flags.n = isNorth;
  flags.s = isSouth;

  auto min_ll = m_min_lakelevel.star(i, j);
  auto max_ll = m_max_lakelevel.star(i, j);

  double min_level = min_ll.ij, max_level = max_ll.ij;

  for (auto d : { North, East, South, West }) {
    if (flags[d]) {
      min_level = min(min_ll[d], min_level);
      max_level = max(min_ll[d], min_level);
    }
  }

  if (min_level != min_ll.ij) {
    setRunMinLevel(min_level, run, lists);
    changed = true;
  }

  if (max_level != max_ll.ij) {
    setRunMaxLevel(max_level, run, lists);
    changed = true;
  }
}

void LakePropertiesCC::startNewRun(int i, int j, int &run_number, int &parent, VecList &lists) {
  ConnectedComponents::startNewRun(i, j, run_number, parent, lists);

  lists["min_ll"][run_number]  = (*m_current_level)(i, j);
  lists["max_ll"][run_number]  = (*m_current_level)(i, j);
}

void LakePropertiesCC::continueRun(int i, int j, int &run_number, VecList &lists) {
  ConnectedComponents::continueRun(i, j, run_number, lists);

  setRunMinLevel((*m_current_level)(i, j), run_number, lists);
  setRunMaxLevel((*m_current_level)(i, j), run_number, lists);
}

} // end of namespace pism
