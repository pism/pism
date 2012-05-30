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

#ifndef _PSGIVEN_H_
#define _PSGIVEN_H_

#include "PGivenClimate.hh"
#include "PISMSurface.hh"
#include "PSModifier.hh"
#include "PISMAtmosphere.hh"

class PSGivenClimate : public PGivenClimate<PSModifier,PISMSurfaceModel>
{
public:
  PSGivenClimate(IceGrid &g, const NCConfigVariable &conf)
    : PGivenClimate<PSModifier,PISMSurfaceModel>(g, conf, NULL)
  {
    temp_name = "ice_surface_temp";
    mass_flux_name = "climatic_mass_balance";
    option_prefix = "-surface_given";
  }
  virtual ~PSGivenClimate() {}

  virtual void attach_atmosphere_model(PISMAtmosphereModel *input)
  { delete input; }

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode ice_surface_mass_flux(IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_temperature(IceModelVec2S &result);
};

#endif /* _PSGIVEN_H_ */
