// Copyright (C) 2012 PISM Authors
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

#ifndef _PAGENERICYEARLYCYCLE_H_
#define _PAGENERICYEARLYCYCLE_H_

#include "PAYearlyCycle.hh"

class Timeseries;

class PACosineYearlyCycle : public PAYearlyCycle {
public:
  PACosineYearlyCycle(IceGrid &g, const NCConfigVariable &conf)
    : PAYearlyCycle(g, conf), A(NULL) {}
  virtual ~PACosineYearlyCycle();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);
  virtual PetscErrorCode temp_time_series(int i, int j, int N,
					  PetscReal *ts, PetscReal *values);
  virtual PetscErrorCode temp_snapshot(IceModelVec2S &result);
protected:
  Timeseries *A;                 // amplitude scaling
};

#endif /* _PAGENERICYEARLYCYCLE_H_ */
