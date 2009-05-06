// Copyright (C) 2009 Andreas Aschwandend and Ed Bueler
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

#ifndef __iceEnthalpyModel_hh
#define __iceEnthalpyModel_hh

#include <petscvec.h>
#include "iceModel.hh"
#include "iceModelVec.hh"


//! Temporary class for development of enthalpy-based polythermal PISM.
/*!
Based for now on Bueler's reading of AB09 = A. Aschwandedn and H. Blatter,
"An enthalpy method for glaciers and ice sheets", submitted??
 */
class IceEnthalpyModel : public IceModel {

public:
  IceEnthalpyModel(IceGrid &g);
  virtual PetscErrorCode initFromFile(const char *);
  virtual PetscErrorCode write_extra_fields(const char filename[]);

protected:
  IceModelVec3  Enth3;

protected:
  virtual PetscErrorCode createVecs();
  
  PetscErrorCode setEnth3toCTSValue();
  PetscErrorCode setEnth3FromTemp_ColdIce();
  
  PetscScalar getAbsTemp(PetscScalar enth, PetscScalar p);

};

#endif

