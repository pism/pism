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

#ifndef _POSEALEVELGRADUAL_H_
#define _POSEALEVELGRADUAL_H_

#include "pism/coupler/SeaLevel.hh"
#include "pism/util/iceModelVec.hh"

namespace pism {
namespace ocean {
namespace sea_level {

class Gradual : public SeaLevel {
public:
  Gradual(IceGrid::ConstPtr g, std::shared_ptr<SeaLevel> in);
  virtual ~Gradual();

protected:
  virtual void update_impl(const Geometry &geometry, double my_t, double my_dt);
  virtual void init_impl(const Geometry &geometry);
  virtual DiagnosticList diagnostics_impl() const;

private:
  bool m_bootstrap;
  double m_max_fill_rate;

  IceModelVec2S m_target_level, m_min_basin, m_max_ll_basin;
  IceModelVec2Int m_expansion_mask;

  void prepareSeaLevel(const IceModelVec2S &target_level,
                       const IceModelVec2S &bed,
                       const IceModelVec2S &lake_level,
                       IceModelVec2Int &mask,
                       IceModelVec2S &min_basin,
                       IceModelVec2S &max_ll_basin,
                       IceModelVec2S &sea_level);
  void gradually_fill(const double dt,
                      const double max_fill_rate,
                      const IceModelVec2S &target_level,
                      const IceModelVec2S &bed,
                      const IceModelVec2S &thk,
                      const IceModelVec2S &min_bed,
                      IceModelVec2S &sea_level);

};

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism
#endif /* _POSEALEVELGRADUAL_H_ */
