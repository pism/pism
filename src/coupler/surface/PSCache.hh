/* Copyright (C) 2013, 2014 PISM Authors
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

#ifndef _PSCACHE_H_
#define _PSCACHE_H_

#include "PSModifier.hh"
#include "iceModelVec.hh"

namespace pism {

class PSCache : public PSModifier {
public:
  PSCache(IceGrid &g, SurfaceModel* in);
  virtual ~PSCache();

  virtual void init(Vars &vars);
  virtual void update(double my_t, double my_dt);
  virtual void ice_surface_mass_flux(IceModelVec2S &result);
  virtual void ice_surface_temperature(IceModelVec2S &result);
  virtual void ice_surface_liquid_water_fraction(IceModelVec2S &result);
  virtual void mass_held_in_surface_layer(IceModelVec2S &result);
  virtual void surface_layer_thickness(IceModelVec2S &result);

  virtual void define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype);
  virtual void write_variables(const std::set<std::string> &vars, const PIO &nc);

  virtual void max_timestep(double t, double &dt, bool &restrict);
protected:
  IceModelVec2S m_mass_flux, m_temperature, m_liquid_water_fraction,
    m_mass_held_in_surface_layer, m_surface_layer_thickness;
  double m_next_update_time;
  unsigned int m_update_interval_years;
};

} // end of namespace pism

#endif /* _PSCACHE_H_ */
