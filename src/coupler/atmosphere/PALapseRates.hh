// Copyright (C) 2011 PISM Authors
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

#ifndef _PALAPSERATES_H_
#define _PALAPSERATES_H_

#include "PLapseRates.hh"
#include "PAModifier.hh"

class PALapseRates : public PLapseRates<PISMAtmosphereModel,PAModifier>
{
public:
  PALapseRates(IceGrid &g, const NCConfigVariable &conf, PISMAtmosphereModel* in)
    : PLapseRates<PISMAtmosphereModel,PAModifier>(g, conf, in)
  {
    precip_lapse_rate = 0;
    option_prefix = "-atmosphere";
  }

  virtual ~PALapseRates() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode mean_precip(IceModelVec2S &result);
  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result);

  virtual PetscErrorCode begin_pointwise_access();
  virtual PetscErrorCode end_pointwise_access();

  virtual PetscErrorCode temp_time_series(int i, int j, int N,
                                          PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result);
protected:
  PetscReal precip_lapse_rate;
};

#endif /* _PALAPSERATES_H_ */

