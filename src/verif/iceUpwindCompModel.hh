// Copyright (C) 2007 Ed Bueler
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

#ifndef __iceUpwindCompModel_hh
#define __iceUpwindCompModel_hh

#include <petsc.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "iceCompModel.hh"

class IceUpwindCompModel : public IceCompModel {

public:
  IceUpwindCompModel(IceGrid &g, ThermoGlenArrIce *i, const char mytest);
  virtual PetscErrorCode initFromOptions();
  virtual PetscErrorCode velocity(bool updateSIAVelocityAtDepth);

};

#endif /* __iceUpwindCompModel_hh */

