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

#ifndef _PAYEARLYCYCLE_H_
#define _PAYEARLYCYCLE_H_

#include "PISMAtmosphere.hh"
#include "iceModelVec.hh"

//! A class containing an incomplete implementation of an atmosphere model
//! based on a temperature parameterization using mean annual and mean July
//! (mean summer) temperatures and a cosine yearly cycle. Uses a stored
//! (constant in time) precipitation field.
class PAYearlyCycle : public PISMAtmosphereModel {
public:
  PAYearlyCycle(IceGrid &g, const NCConfigVariable &conf)
    : PISMAtmosphereModel(g, conf) {}
  virtual PetscErrorCode init(PISMVars &vars);
  virtual void add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
  //! This method implements the parameterization.
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt) = 0;
  virtual PetscErrorCode mean_precipitation(IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result);
  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result);
protected:
  PISMVars *variables;
  PetscScalar snow_temp_july_day;
  string reference, precip_filename;
  IceModelVec2S air_temp_mean_annual, air_temp_mean_july, precipitation;
  NCSpatialVariable air_temp_snapshot;
};

#endif /* _PAYEARLYCYCLE_H_ */
