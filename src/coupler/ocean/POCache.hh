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

#ifndef _POCACHE_H_
#define _POCACHE_H_

#include "POModifier.hh"
#include "iceModelVec.hh"

namespace pism {

class POCache : public POModifier {
public:
  POCache(IceGrid &g, const Config &conf, OceanModel* in);
  virtual ~POCache();

  virtual PetscErrorCode init(Vars &vars);
  virtual PetscErrorCode update(double my_t, double my_dt);

  virtual PetscErrorCode sea_level_elevation(double &result);
  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);
  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);
  virtual PetscErrorCode melange_back_pressure_fraction(IceModelVec2S &result);

  virtual PetscErrorCode define_variables(std::set<std::string> vars, const PIO &nc,
                                          IO_Type nctype);
  virtual PetscErrorCode write_variables(std::set<std::string> vars, const PIO &nc);
  virtual PetscErrorCode max_timestep(double t, double &dt, bool &restrict);
protected:
  IceModelVec2S m_shelf_base_temperature, m_shelf_base_mass_flux,
    m_melange_back_pressure_fraction;
  double m_sea_level;
  double m_next_update_time;
  unsigned int m_update_interval_years;
  PetscErrorCode allocate_POCache();
};

} // end of namespace pism

#endif /* _POCACHE_H_ */
