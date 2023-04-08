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

#ifndef PISM_LAKEACCUMULATOR_H
#define PISM_LAKEACCUMULATOR_H

#include "pism/util/connected_components_lakecc.hh"

namespace pism {

class LakeAccumulatorCCSerial : public ConnectedComponentsSerial {
public:
  LakeAccumulatorCCSerial(IceGrid::ConstPtr g);
  virtual ~LakeAccumulatorCCSerial() = default;
  void init(const IceModelVec2S &lake_level);
  void accumulate(const IceModelVec2S &in, IceModelVec2S &result);

protected:
  virtual bool ForegroundCond(int i, int j) const;

private:
  bool m_initialized;
  VecList m_lists;
  int m_run_number;

  void prepare_mask(const IceModelVec2S &lake_level);
};

} // end of namespace pism

#endif /* PISM_LAKEACCUMULATOR_H */
