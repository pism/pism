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
  virtual void update_impl(const Geometry &geometry, double my_t, double my_dt);
  virtual void init_impl(const Geometry &geometry);
  virtual bool expandMargins_impl() const;
  virtual DiagnosticList diagnostics_impl() const;

private:
  GeometryCalculator m_gc;
  std::string m_option;
  double m_lake_level_min, m_lake_level_max, m_lake_level_dh, m_drho, m_icefree_thickness;
  bool m_filter_map, m_keep_existing_lakes, m_check_sl_diagonal;
  int m_n_filter;
  IceModelVec2Int m_valid_mask;

  void do_lake_update(const IceModelVec2S &bed, const IceModelVec2S &thk, const IceModelVec2S &sea_level, const IceModelVec2S &old_lake_level);
  void do_filter_map();
  void prepare_mask_validity(const IceModelVec2S &thk, const IceModelVec2S &bed, const IceModelVec2S &old_lake_level, IceModelVec2Int &valid_mask);
  void update_mask_sl_diagonal(const IceModelVec2S &sl, const IceModelVec2S &bed, const IceModelVec2S &thk, IceModelVec2Int &mask);
};

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism
#endif /* _POLAKECC_H_ */