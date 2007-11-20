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

#ifndef __iceEISplModel_hh
#define __iceEISplModel_hh

#include "../base/grid.hh"
#include "../base/materials.hh"
#include "iceEISModel.hh"

//! This derived class does the plastic till and SSA modification of EISMINT II experiment I.  
class IceEISplModel : public IceEISModel {

public:
  IceEISplModel(IceGrid &g, IceType &i);
  virtual PetscErrorCode initFromOptions();
    
protected:
  PetscScalar*    phi_list;

  PetscErrorCode resetAccum();
  PetscErrorCode setTillProperties();

};

#endif /* __iceEISplModel_hh */
