// Copyright (C) 2008 Ed Bueler
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

#include <cstring>
#include <cstdlib>
#include <petscda.h>
#include <netcdf.h>
#include "nc_util.hh"

#ifndef __IceModelVec_hh
#define __IceModelVec_hh

//! Abstract class for reading, writing, allocating, and accessing a DA-based PETSc Vec from within IceModel.
class IceModelVec {
public:
  IceModelVec();
  ~IceModelVec();

  virtual PetscErrorCode  create(IceGrid* mygrid);
  virtual PetscErrorCode  destroy();

  virtual PetscErrorCode  setidNC(int my_ncid);
  virtual PetscErrorCode  setAttrsNC(
             const char my_long_name[], const char my_units[], const char my_pism_intent[]);
  virtual PetscErrorCode  setAttrsCFstandardNC(
             const char my_long_name[], const char my_units[], const char my_pism_intent[],
             const char my_standard_name[]);
  virtual PetscErrorCode  readFromNC();
  virtual PetscErrorCode  writeToNC();

  virtual PetscErrorCode  needAccessVals();
  virtual PetscErrorCode  doneAccessVals();
  virtual PetscErrorCode  beginComm();
  virtual PetscErrorCode  endComm();

  // FIXME:  make this protected once IceModelVec fully implemented!
  Vec          v;

protected:
  char         long_name[PETSC_MAX_PATH_LEN],
               standard_name[PETSC_MAX_PATH_LEN],
               units[PETSC_MAX_PATH_LEN],
               intent[PETSC_MAX_PATH_LEN];
  PetscTruth   allocated, has_standard_name;
  
  IceGrid      *grid;
  DA           da;
  int          ncid;

  virtual PetscErrorCode  checkAllocated();
  virtual PetscErrorCode  checkHaveArray();

  // FIXME: consider adding 
  //   virtual PetscErrorCode  checkSelfOwnsIt(const PetscInt i, const PetscInt j);
  //   virtual PetscErrorCode  checkSelfOwnsItGhosted(const PetscInt i, const PetscInt j);
};


//! Class for a 3d DA-based Vec for ice scalar quantities in IceModel.
class IceModelVec3 : public IceModelVec {
public:
  IceModelVec3();
  virtual PetscErrorCode  destroy();
  virtual PetscErrorCode  create(IceGrid* mygrid);

public:
  PetscScalar     getValZ(const PetscInt i, const PetscInt j, const PetscScalar z);
  PetscErrorCode  getValColumn(const PetscInt i, const PetscInt j, const PetscInt nlevels, 
                               PetscScalar *levels, PetscScalar *vals);
  PetscErrorCode  setValColumn(const PetscInt i, const PetscInt j, const PetscInt nlevels, 
                               PetscScalar *levels, PetscScalar *vals);

  PetscErrorCode  getHorSlice(Vec &gslice, const PetscScalar z);
  PetscErrorCode  getSurfaceValues(Vec &gsurf, PetscScalar **H);
 
protected:
  PetscScalar     ***array;
  PetscInt        Mz;
  PetscScalar     dz, *levels;
  
  PetscErrorCode  isLegalLevel(const PetscScalar z);
};


#endif /* __IceModelVec_hh */

