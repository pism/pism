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

#ifndef _PAMODIFIER_H_
#define _PAMODIFIER_H_

#include "PISMAtmosphere.hh"

class PAModifier : public Modifier<PISMAtmosphereModel>
{
public:
  PAModifier(IceGrid &g, const NCConfigVariable &conf, PISMAtmosphereModel* in)
    : Modifier<PISMAtmosphereModel>(g, conf, in) {}
  virtual ~PAModifier() {}

  virtual PetscErrorCode mean_precipitation(IceModelVec2S &result)
  {
    if (input_model != NULL) {
      PetscErrorCode ierr = input_model->mean_precipitation(result); CHKERRQ(ierr);
    }
    return 0;
  }

  virtual PetscErrorCode mean_annual_temp(IceModelVec2S &result)
  {
    if (input_model != NULL) {
      PetscErrorCode ierr = input_model->mean_annual_temp(result); CHKERRQ(ierr);
    }
    return 0;
  }

  virtual PetscErrorCode begin_pointwise_access()
  {
    if (input_model != NULL) {
      PetscErrorCode ierr = input_model->begin_pointwise_access(); CHKERRQ(ierr);
    }
    return 0;
  }

  virtual PetscErrorCode end_pointwise_access()
  {
    if (input_model != NULL) {
      PetscErrorCode ierr = input_model->end_pointwise_access(); CHKERRQ(ierr);
    }
    return 0;
  }

  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values)
  {
    if (input_model != NULL) {
      PetscErrorCode ierr = input_model->temp_time_series(i, j, N, ts, values); CHKERRQ(ierr);
    }
    return 0;
  }

  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result)
  {
    if (input_model != NULL) {
      PetscErrorCode ierr = input_model->temp_snapshot(result); CHKERRQ(ierr);
    }
    return 0;
  }
};

#endif /* _PAMODIFIER_H_ */
