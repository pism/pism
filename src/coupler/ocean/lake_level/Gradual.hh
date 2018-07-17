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

#ifndef _LAKE_LEVEL_GRADUAL_H_
#define _LAKE_LEVEL_GRADUAL_H_

#include "pism/coupler/LakeLevel.hh"

namespace pism {
namespace ocean {
namespace lake_level {

class Gradual : public LakeLevel {
public:
  Gradual(IceGrid::ConstPtr g, std::shared_ptr<LakeLevel> in);

private:
  void update_impl(const Geometry &geometry, double t, double dt);
  void init_impl(const Geometry &geometry);
  IceModelVec2S m_target_level, m_min_level, m_max_level, m_min_basin,
                m_lake_area, m_lake_mass_input_discharge, m_lake_mass_input_basal,
                m_lake_mass_input_total, m_lake_fill_rate;
  IceModelVec2Int m_expansion_mask;
  bool m_use_const_fill_rate;
protected:
  double m_max_lake_fill_rate;
  void prepareLakeLevel(const IceModelVec2S &target_level,
                        const IceModelVec2S &bed,
                        IceModelVec2S &min_level,
                        IceModelVec2S &max_level,
                        IceModelVec2S &min_basin,
                        IceModelVec2Int &mask);
  void gradually_fill(double dt,
                      const IceModelVec2S &target_level,
                      const IceModelVec2S &bed,
                      const IceModelVec2S &thk,
                      const IceModelVec2S &sea_level,
                      const IceModelVec2S &min_level,
                      const IceModelVec2S &max_level,
                      const IceModelVec2S &min_bed);
  void compute_fill_rate(double dt);
  DiagnosticList diagnostics_impl() const;
};

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism

#endif /* _LAKE_LEVEL_GRADUAL_H_ */
