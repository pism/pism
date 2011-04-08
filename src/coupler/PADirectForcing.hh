// Copyright (C) 2010, 2011 Constantine Khroulev
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

#ifndef _PADIRECTFORCING_H_
#define _PADIRECTFORCING_H_

#include "PISMAtmosphere.hh"
#include "iceModelVec2T.hh"

class PADirectForcing : public PISMAtmosphereModel {
public:
  PADirectForcing(IceGrid &g, const NCConfigVariable &conf);
  virtual ~PADirectForcing();
  virtual PetscErrorCode max_timestep(PetscReal t_years, PetscReal &dt_years);
  virtual PetscErrorCode init(PISMVars &vars);
  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc, nc_type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
  virtual PetscErrorCode update(PetscReal t_years, PetscReal dt_years);
  virtual PetscErrorCode mean_precip(PetscReal t_years, PetscReal dt_years,
				     IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(PetscReal t_years, PetscReal dt_years,
					  IceModelVec2S &result); 
  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();  
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode temp_snapshot(PetscReal t_years, PetscReal dt_years,
				       IceModelVec2S &result);
protected:
  IceModelVec2T *temp, *precip;
};

#endif /* _PADIRECTFORCING_H_ */
