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

#ifndef _PADTFORCING_H_
#define _PADTFORCING_H_

#include "PScalarForcing.hh"
#include "PAModifier.hh"

class PA_delta_T : public PScalarForcing<PISMAtmosphereModel,PAModifier>
{
public:
  PA_delta_T(IceGrid &g, const NCConfigVariable &conf, PISMAtmosphereModel* in);
  virtual ~PA_delta_T() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result);

  virtual PetscErrorCode temp_time_series(int i, int j, int N,
                                          PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result);

  virtual void add_vars_to_output(string keyword,
                                  map<string,NCSpatialVariable> &result);

  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc,
                                          PISM_IO_Type nctype);

  virtual PetscErrorCode write_variables(set<string> vars, string filename);

protected:
  NCSpatialVariable air_temp, precipitation;
};


#endif /* _PADTFORCING_H_ */
