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

#include <petsc.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "iceEISModel.hh"

//! This derived class does the plastic till and SSA modification of EISMINT II experiment I.  
class IceEISplModel : public IceEISModel {

public:
  IceEISplModel(IceGrid &g, IceType &i);
  virtual PetscErrorCode initFromOptions();
    
protected:
  PetscScalar              stream_width;
  static const PetscScalar DEFAULT_STREAM_WIDTH;

  static const PetscInt    phi_list_length; // = 5
  PetscScalar*             phi_list; // list of phi_list_length
  static const PetscScalar DEFAULT_TILL_PHI_LAKE;
  static const PetscScalar DEFAULT_TILL_PHI_STRONG;
  static const PetscScalar DEFAULT_TILL_PHI_UPSTREAM;
  static const PetscScalar DEFAULT_TILL_PHI_DOWNSTREAM;
  static const PetscScalar DEFAULT_TILL_PHI_OCEAN;

  PetscErrorCode resetAccum();
  PetscErrorCode setTillProperties();

};

#endif /* __iceEISplModel_hh */
