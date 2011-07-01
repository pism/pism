// Copyright (C) 2011 Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef _PODIRECTFORCING_H_
#define _PODIRECTFORCING_H_

#include "PASDirectForcing.hh"
#include "PISMOcean.hh"

class PODirectForcing : public PDirectForcing<PISMOceanModel>
{
public:
  PODirectForcing(IceGrid &g, const NCConfigVariable &conf)
    : PDirectForcing<PISMOceanModel>(g, conf)
  {
    temp_name       = "shelfbtemp";
    mass_flux_name  = "shelfbmassflux";
    bc_option_name = "-ocean_bc_file";
  }

  virtual ~PODirectForcing() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);

  virtual PetscErrorCode sea_level_elevation(PetscReal &result) {
    result = sea_level;
    return 0;
  }

  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);

  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);
};


#endif /* _PODIRECTFORCING_H_ */
