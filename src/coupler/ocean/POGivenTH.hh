// Copyright (C) 2011, 2012, 2014 PISM Authors
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
  POGivenTH(IceGrid &g, const PISMConfig &conf);
  virtual ~POGivenTH();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(double my_t, double my_dt);

  virtual PetscErrorCode sea_level_elevation(double &result);

  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);
  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);

private:
  IceModelVec2S shelfbtemp, shelfbmassflux;
  IceModelVec2S *ice_thickness;
  IceModelVec2T *theta_ocean, *salinity_ocean;

  PetscErrorCode calc_shelfbtemp_shelfbmassflux();

  PetscErrorCode btemp_bmelt_3eqn(double rhow, double rhoi,
                                  double sal_ocean, double temp_insitu, double zice,
                                  double &temp_base, double &meltrate);

  PetscErrorCode adiabatic_temperature_gradient(double salinity, double temp_insitu, double pressure, double &adlprt_out);
  PetscErrorCode potential_temperature(double salinity,double temp_insitu,double pressure,
                                       double reference_pressure, double& thetao);
  PetscErrorCode insitu_temperature(double salinity, double thetao,
                                    double pressure,double reference_pressure,
                                    double &temp_insitu_out);
  PetscErrorCode allocate_POGivenTH();
};

#endif /* _POGIVENTH_H_ */
