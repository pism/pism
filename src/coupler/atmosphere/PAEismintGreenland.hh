// Copyright (C) 2007-2011, 2013, 2014 Ed Bueler and Nathan Shemonski and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#ifndef __PA_EISMINT_Greenland
#define __PA_EISMINT_Greenland

#include "PAYearlyCycle.hh"

class PA_EISMINT_Greenland : public PAYearlyCycle {
public:
  PA_EISMINT_Greenland(IceGrid &g, const PISMConfig &conf);
  virtual ~PA_EISMINT_Greenland();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);
protected:
  virtual PetscReal greenhouse_shift(PetscReal my_t, PetscReal my_dt);
  bool do_greenhouse_warming;
  PetscReal greenhouse_warming_start_year;
  IceModelVec2S *lat, *surfelev;
};

#endif	// __PA_EISMINT_Greenland
