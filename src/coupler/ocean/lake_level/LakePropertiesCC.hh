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

#ifndef PISM_LAKEPROPERTIESCC_H
#define PISM_LAKEPROPERTIESCC_H

#include "pism/util/connected_components_lakecc.hh"

namespace pism {

/*!
 * LakePropertiesCC collects the minimum and maximum current water level of each lake
 * basin.
 */
class LakePropertiesCC : public ConnectedComponents {
public:
  LakePropertiesCC(IceGrid::ConstPtr g, const IceModelVec2S &target_level,
                   const IceModelVec2S &lake_level);
  virtual ~LakePropertiesCC() = default;

  void getLakeProperties(IceModelVec2S &min_level, IceModelVec2S &max_level);

private:
  const IceModelVec2S *m_target_level, *m_current_level;
  IceModelVec2S m_min_lakelevel, m_max_lakelevel;

  void setRunMinLevel(double level, int run, VecList &lists);
  void setRunMaxLevel(double level, int run, VecList &lists);

  void init_VecList(VecList &lists, unsigned int length);
  bool ForegroundCond(int i, int j) const;
  void labelMask(int run_number, const VecList &lists);
  void treatInnerMargin(int i, int j,
                                bool isNorth, bool isEast, bool isSouth, bool isWest,
                                VecList &lists, bool &changed);
  void startNewRun(int i, int j, int &run_number, int &parent, VecList &lists);
  void continueRun(int i, int j, int &run_number, VecList &lists);
};

} // end of namespace pism

#endif /* PISM_LAKEPROPERTIESCC_H */
