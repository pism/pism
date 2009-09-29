// Copyright (C) 2009 Constantine Khroulev
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

#include <petscda.h>
#include "../base/iceModel.hh"

#ifndef __iceUnitModel_hh
#define __iceUnitModel_hh

class IceUnitModel : public IceModel {
public:
  IceUnitModel(IceGrid &g) : IceModel(g) {}
  PetscErrorCode set_grid_defaults();
  PetscErrorCode set_vars_from_options();
  PetscErrorCode run();
  PetscErrorCode writeFiles(const char filename[]);

  PetscErrorCode testIceModelVec3();
  PetscErrorCode testIceModelVec3Bedrock();
};

#endif // 
