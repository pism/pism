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
  if (isLake(level)) {
    if (isLake(lists["min_ll"][run])) {
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
  if (isLake(level)) {
    if (isLake(lists["max_ll"][run])) {
      level = std::max(level, lists["max_ll"][run]);
    }
    lists["max_ll"][run] = level;
  }
}


bool LakePropertiesCC::ForegroundCond(int i, int j) const {
  const double target  = (*m_target_level)(i, j);

  return isLake(target);
}

void LakePropertiesCC::labelMask(int run_number, const VecList &lists) {
  IceModelVec::AccessList list;
  list.add(m_masks.begin(), m_masks.end());

  const auto
    &i_vec   = lists.find("i")->second,
    &j_vec   = lists.find("j")->second,
    &len_vec = lists.find("lengths")->second,
    &parents = lists.find("parents")->second,
    &min_ll  = lists.find("min_ll")->second,
    &max_ll  = lists.find("max_ll")->second;

  for (int k = 0; k <= run_number; ++k) {
    const int label = connected_components::trackParentRun(k, parents);
    const double min_ll_label = min_ll[label],
                 max_ll_label = max_ll[label];
    int j = j_vec[k];
    for (unsigned int n = 0; n < len_vec[k]; ++n) {
      const int i = i_vec[k] + n;

      m_mask_run(i, j) = label;
      m_min_lakelevel(i, j) = min_ll_label;
      m_max_lakelevel(i, j) = max_ll_label;
    }
  }
}

void LakePropertiesCC::treatInnerMargin(int i, int j,
                                        const bool isNorth, const bool isEast, const bool isSouth, const bool isWest,
                                        VecList &lists, bool &changed) {
  ConnectedComponents::treatInnerMargin(i, j, isNorth, isEast, isSouth, isWest, lists, changed);

  int run = m_mask_run.as_int(i, j);
  if (run > 0) {
    StarStencil<double> min_ll_star  = m_min_lakelevel.star(i, j),
                        max_ll_star  = m_max_lakelevel.star(i, j);

    double min_level = min_ll_star.ij,
           max_level = max_ll_star.ij;

    if (isWest) {
      if (isLake(min_ll_star.w) and ((min_ll_star.w < min_level) or not isLake(min_level))) {
        min_level = min_ll_star.w;
      }
      if (isLake(max_ll_star.w) and ((max_ll_star.w > max_level) or not isLake(max_level))) {
        max_level = max_ll_star.w;
      }
    }
    if (isNorth) {
      if (isLake(min_ll_star.n) and ((min_ll_star.n < min_level) or not isLake(min_level))) {
        min_level = min_ll_star.n;
      }
      if (isLake(max_ll_star.n) and ((max_ll_star.n > max_level) or not isLake(max_level))) {
        max_level = max_ll_star.n;
      }
    }
    if (isEast) {
      if (isLake(min_ll_star.e) and ((min_ll_star.e < min_level) or not isLake(min_level))) {
        min_level = min_ll_star.e;
      }
      if (isLake(max_ll_star.e) and ((max_ll_star.e > max_level) or not isLake(max_level))) {
        max_level = max_ll_star.e;
      }
    }
    if (isSouth) {
      if (isLake(min_ll_star.s) and ((min_ll_star.s < min_level) or not isLake(min_level))) {
        min_level = min_ll_star.s;
      }
      if (isLake(max_ll_star.s) and ((max_ll_star.s > max_level) or not isLake(max_level))) {
        max_level = max_ll_star.s;
      }
    }
    if (min_level != min_ll_star.ij) {
      setRunMinLevel(min_level, run, lists);
      changed = true;
    }
    if (max_level != max_ll_star.ij) {
      setRunMaxLevel(max_level, run, lists);
      changed = true;
    }
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
