// Copyright (C) 2011, 2012 PISM Authors
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

#ifndef _POGIVENTH_H_
#define _POGIVENTH_H_

#include "PGivenClimate.hh"
#include "POModifier.hh"

class POGivenTH : public PGivenClimate<POModifier,PISMOceanModel>
{

public:
  POGivenTH(IceGrid &g, const NCConfigVariable &conf);
  virtual ~POGivenTH();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode sea_level_elevation(PetscReal &result);

  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);
  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);

  virtual PetscErrorCode calc_shelfbtemp_shelfbmassflux();

  virtual PetscErrorCode btemp_bmelt_3eqn( PetscReal rhow, PetscReal rhoi,
                                                        PetscReal sal_ocean, PetscReal temp_insitu, PetscReal zice,
                                                        PetscReal &temp_base, PetscReal &meltrate);

  virtual PetscErrorCode adiabatic_temperature_gradient(PetscReal salinity, PetscReal temp_insitu, PetscReal pressure, PetscReal &adlprt_out);
  virtual PetscErrorCode potential_temperature(PetscReal salinity,PetscReal temp_insitu,PetscReal pressure,
                                                    PetscReal reference_pressure, PetscReal& thetao);
  virtual PetscErrorCode insitu_temperature(PetscReal salinity, PetscReal thetao,
                                                 PetscReal pressure,PetscReal reference_pressure,
                                                 PetscReal &temp_insitu_out);

protected:
  IceModelVec2T *shelfbtemp, *shelfbmassflux;
  IceModelVec2S *ice_thickness;
  IceModelVec2T *theta_ocean, *salinity_ocean;

private:
  PetscErrorCode allocate_POGivenTH();
};

#endif /* _POGIVENTH_H_ */
