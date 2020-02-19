/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
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

#ifndef _POLAKECC_H_
#define _POLAKECC_H_

#include "pism/util/Mask.hh"
#include "pism/coupler/LakeLevel.hh"

namespace pism {
namespace ocean {
namespace lake_level {
class LakeCC : public LakeLevel {
public:
  LakeCC(IceGrid::ConstPtr g);
  ~LakeCC();

protected:
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void update_impl(const Geometry &geometry, double t, double dt);
  virtual void init_impl(const Geometry &geometry);
  virtual bool expandMargins_impl() const;

  IceModelVec2S m_target_level;

private:
  GeometryCalculator m_gc;
  std::string m_option;

  double m_next_update_time;
  int m_max_update_interval_years,
      m_patch_iter;

  bool checkOceanBasinsVanished(const IceModelVec2S &bed,
                                const IceModelVec2S &old_sl,
                                const IceModelVec2S &new_sl);
  bool iterativelyPatchTargetLevel(const IceModelVec2S &bed,
                                   const IceModelVec2S &thk,
                                   const IceModelVec2S &sl);
  unsigned int patch_lake_levels(const IceModelVec2S &bed,
                                 const IceModelVec2S &thk,
                                 const IceModelVec2S &sea_level,
                                 IceModelVec2S &lake_level);
};

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism
#endif /* _POLAKECC_H_ */