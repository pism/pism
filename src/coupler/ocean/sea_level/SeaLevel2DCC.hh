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

#ifndef _POSEALEVEL2DCC_H_
#define _POSEALEVEL2DCC_H_

#include "pism/coupler/SeaLevel.hh"
#include "pism/util/iceModelVec.hh"

namespace pism {
namespace ocean {
namespace sea_level {

class SeaLevel2DCC : public SeaLevel {
public:
  SeaLevel2DCC(IceGrid::ConstPtr g, std::shared_ptr<SeaLevel> in);
  virtual ~SeaLevel2DCC();

protected:
  virtual MaxTimestep max_timestep_impl(double t) const;
  virtual void update_impl(const Geometry &geometry, double my_t, double my_dt);
  virtual void init_impl(const Geometry &geometry);
  virtual const bool expandMargins_impl() const;

private:
  std::string m_option_prefix;
  double m_next_update_time;
  double m_drho, m_offset;
  int m_update_interval_years;
  bool m_update_periodic, m_update_passive, m_update_startup, m_update;
  IceModelVec2Int m_mask;

  void process_options();
  void do_sl_mask_update(const IceModelVec2S &bed, const IceModelVec2S &thk);
};

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism
#endif /* _POSEALEVEL2DCC_H_ */
