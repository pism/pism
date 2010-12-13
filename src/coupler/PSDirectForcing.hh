// Copyright (C) 2010 Constantine Khroulev
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

#ifndef _PSDIRECTFORCING_H_
#define _PSDIRECTFORCING_H_

#include "PISMSurface.hh"
#include "iceModelVec2T.hh"

class PSDirectForcing : public PISMSurfaceModel
{
public:
  PSDirectForcing(IceGrid &g, const NCConfigVariable &conf)
    : PISMSurfaceModel(g, conf) {}

  virtual ~PSDirectForcing() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);

  virtual PetscErrorCode max_timestep(PetscReal t_years, PetscReal &dt_years);

  virtual void add_vars_to_output(string keyword, set<string> &result);

  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc, nc_type nctype);

  virtual PetscErrorCode write_variables(set<string> vars, string filename);

  virtual PetscErrorCode ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
					       IceModelVec2S &result);

  virtual PetscErrorCode ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
						 IceModelVec2S &result);

  virtual void attach_atmosphere_model(PISMAtmosphereModel *input)
  { delete input; }
protected:
  IceModelVec2T temperature, mass_flux;
  PetscReal bc_period;

  PetscReal my_mod(PetscReal input);
};

#endif /* _PSDIRECTFORCING_H_ */
