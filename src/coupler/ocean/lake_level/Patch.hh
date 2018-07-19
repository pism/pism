/* Copyright (C) 2018 PISM Authors
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

#ifndef _LAKE_LEVEL_PATCH_H_
#define _LAKE_LEVEL_PATCH_H_

#include "pism/coupler/LakeLevel.hh"

namespace pism {
namespace ocean {
namespace lake_level {

class Patch : public LakeLevel {
public:
  Patch(IceGrid::ConstPtr g, std::shared_ptr<LakeLevel> in);
  ~Patch();

protected:
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void update_impl(const Geometry &geometry, double t, double dt);
  virtual void init_impl(const Geometry &geometry);

private:
  double m_next_update_time, m_last_update;
  int m_min_update_interval_years, m_patch_iter;

  unsigned int patch_lake_levels(const IceModelVec2S &bed,
                                 const IceModelVec2S &thk,
                                 const IceModelVec2S &sea_level);

};

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism

#endif /* _LAKE_LEVEL_PATCH_H_ */
