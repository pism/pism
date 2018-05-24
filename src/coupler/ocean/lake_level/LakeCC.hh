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
#include "pism/util/iceModelVec.hh"

namespace pism {
namespace ocean {
namespace lake_level {
class LakeCC : public LakeLevel {
public:
  LakeCC(IceGrid::ConstPtr g);
  virtual ~LakeCC();

protected:
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void update_impl(const Geometry &geometry, double my_t, double my_dt);
  virtual void init_impl(const Geometry &geometry);

protected:
  GeometryCalculator m_gc;
  IceModelVec2S m_target_lake_level;
  IceModelVec2S *m_ll_ptr;
  std::string m_option;
  double m_next_update_time;
  int m_min_update_interval_years, m_patch_iter;
  double m_lake_level_min, m_lake_level_max, m_lake_level_dh, m_drho, m_last_update, m_icefree_thickness, m_max_lake_fill_rate;
  bool m_update_patch, m_update_gradual, m_update_periodic, m_update_passive, m_update;

private:
  void process_options();
  unsigned int patch_lake_levels(const IceModelVec2S *bed, const IceModelVec2S *thk, const IceModelVec2S *sea_level);
  void gradually_fill(double dt, const IceModelVec2S *bed, const IceModelVec2S *thk, const IceModelVec2S *sea_level);
  void do_lake_update(const IceModelVec2S *bed, const IceModelVec2S *thk, const IceModelVec2S *sea_level);
  void prepare_mask_validity(const IceModelVec2S *thk, IceModelVec2Int &valid_mask);
};

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism
#endif /* _POLAKECC_H_ */
