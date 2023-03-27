/* Copyright (C) 2013, 2014, 2015, 2016, 2023 PISM Authors
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

#ifndef PISM_OCEAN_LAKE_CC_H
#define PISM_OCEAN_LAKE_CC_H

#include "pism/geometry/Geometry.hh"

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
  virtual DiagnosticList diagnostics_impl() const;

  IceModelVec2S m_target_level,
                m_fill_rate;

private:
  GeometryCalculator m_gc;
  std::string m_option;

  IceModelVec2S m_topg_overlay;

  double m_next_update_time,
         m_icefree_thickness,
         m_drho,
         m_lake_level_min,
         m_lake_level_max,
         m_lake_level_dh,
         m_max_lake_fill_rate;

  int m_max_update_interval_years,
      m_patch_iter,
      m_filter_size;

  bool m_keep_existing_lakes,
       m_check_sl_diagonal,
       m_use_topg_overlay,
       m_use_const_fill_rate;

  bool checkOceanBasinsVanished(const IceModelVec2S &bed,
                                const IceModelVec2S &old_sl,
                                const IceModelVec2S &new_sl);
  bool iterativelyPatchTargetLevel(const IceModelVec2S &bed,
                                   const IceModelVec2S &thk,
                                   const IceModelVec2S &sl,
                                   IceModelVec2S &target_level);
  unsigned int patch_lake_levels(const IceModelVec2S &bed,
                                 const IceModelVec2S &thk,
                                 const IceModelVec2S &sea_level,
                                 IceModelVec2S &lake_level);
  void updateLakeCC(const IceModelVec2S &bed,
                    const IceModelVec2S &thk,
                    const IceModelVec2S &sea_level,
                    const IceModelVec2S &eff_lake_level,
                    IceModelVec2S &lake_level);
  void compute_fill_rate(double dt,
                         const IceModelVec2S &lake_level,
                         const IceModelVec2S &bmb,
                         const IceModelVec2S &tc_calving,
                         const IceModelVec2S &tc_frontal_melt,
                         const IceModelVec2S &tc_forced_retreat,
                         IceModelVec2S &lake_fill_rate);
  void updateLakeLevelMinMax(const IceModelVec2S &lake_level,
                             const IceModelVec2S &target_level,
                             IceModelVec2S &min_level,
                             IceModelVec2S &max_level);
  bool prepareLakeLevel(const IceModelVec2S &target_level,
                        const IceModelVec2S &bed,
                        const IceModelVec2S &thk,
                        const IceModelVec2S &min_level,
                        const IceModelVec2S &old_ll,
                        const IceModelVec2S &old_sl,
                        IceModelVec2S &min_basin,
                        IceModelVec2S &lake_level);
  void gradually_fill(double dt,
                      double max_fill_rate,
                      const IceModelVec2S &target_level,
                      const IceModelVec2S &bed,
                      const IceModelVec2S &thk,
                      const IceModelVec2S &sea_level,
                      const IceModelVec2S &min_level,
                      const IceModelVec2S &max_level,
                      const IceModelVec2S &min_bed,
                      const IceModelVec2S &fill_rate,
                      IceModelVec2S &lake_level);

};

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism
#endif /* PISM_OCEAN_LAKE_CC_H */
