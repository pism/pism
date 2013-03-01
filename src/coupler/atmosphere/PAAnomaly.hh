// Copyright (C) 2011, 2012, 2013 PISM Authors
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

#ifndef _PAANOMALY_H_
#define _PAANOMALY_H_

#include "PGivenClimate.hh"
#include "PAModifier.hh"

//! \brief Reads and uses air_temp and precipitation anomalies from a file.
class PAAnomaly : public PGivenClimate<PAModifier,PISMAtmosphereModel>
{
public:
  PAAnomaly(IceGrid &g, const NCConfigVariable &conf, PISMAtmosphereModel* in);
  virtual ~PAAnomaly();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode mean_precipitation(IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result); 
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result);

  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode precip_time_series(int i, int j, int N,
					    PetscReal *ts, PetscReal *values);

  virtual void add_vars_to_output(string keyword,
                                  map<string,NCSpatialVariable> &result);

  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc,
                                          PISM_IO_Type nctype);

  virtual PetscErrorCode write_variables(set<string> vars, string file);

protected:
  vector<PetscReal> ts_mod, ts_values;
  NCSpatialVariable air_temp, precipitation;
private:
  PetscErrorCode allocate_PAAnomaly();
};

#endif /* _PAANOMALY_H_ */
