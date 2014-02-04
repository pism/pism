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

class POCache : public POModifier {
public:
  POCache(IceGrid &g, const PISMConfig &conf, PISMOceanModel* in);
  virtual ~POCache();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode sea_level_elevation(PetscReal &result);
  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);
  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);

  virtual PetscErrorCode define_variables(std::set<std::string> vars, const PIO &nc,
                                          PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(std::set<std::string> vars, const PIO &nc);
protected:
  IceModelVec2S m_shelf_base_temperature, m_shelf_base_mass_flux;
  double m_sea_level;
  double m_next_update_time;
  unsigned int m_update_interval_years;
  PetscErrorCode allocate_POCache();
};

#endif /* _POCACHE_H_ */
