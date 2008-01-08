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
  virtual ~IceModelVec();

  virtual PetscErrorCode  create(IceGrid &mygrid, const char my_varname[]);
  virtual PetscErrorCode  destroy();

  virtual PetscErrorCode  setVaridNC(const int my_varid);
  virtual PetscErrorCode  setAttrsNC(const int my_varid,
             const char my_long_name[], const char my_units[], const char my_pism_intent[]);
  virtual PetscErrorCode  setAttrsCFstandardNC(const int my_varid,
             const char my_long_name[], const char my_units[], const char my_pism_intent[],
             const char my_standard_name[]);

  virtual PetscErrorCode  readVecNC(const int ncid, PetscTruth *exists);

  virtual PetscErrorCode  writeAttrsNC(const int ncid);
  virtual PetscErrorCode  putVecNC(const int ncid, Vec g, const int *s, const int *c, int dims, 
                                   void *a_mpi, int a_size);

  virtual PetscErrorCode  needAccessToVals();
  virtual PetscErrorCode  doneAccessToVals();
  virtual PetscErrorCode  beginGhostComm();
  virtual PetscErrorCode  endGhostComm();

  virtual PetscErrorCode  setToConstant(const PetscScalar c);
 
  // FIXME:  make this protected once IceModelVec fully implemented!
  Vec          v;

protected:
  char         varname[PETSC_MAX_PATH_LEN],       // usually the name of the NetCDF variable (unless
                                                  //   there is no corresponding NetCDF var)
               long_name[PETSC_MAX_PATH_LEN],     // NetCDF attribute
               units[PETSC_MAX_PATH_LEN],         // NetCDF attribute
               pism_intent[PETSC_MAX_PATH_LEN],   // NetCDF attribute
               standard_name[PETSC_MAX_PATH_LEN]; // NetCDF attribute; sometimes specified in CF convention
  PetscTruth   allocated, has_standard_name;
  
  IceGrid      *grid;
  DA           *da;
  
  int          varid_nc;

  virtual PetscErrorCode  checkAllocated();
  virtual PetscErrorCode  checkHaveArray();

  // FIXME: consider adding 
  //   virtual PetscErrorCode  checkSelfOwnsIt(const PetscInt i, const PetscInt j);
  //   virtual PetscErrorCode  checkSelfOwnsItGhosted(const PetscInt i, const PetscInt j);
};


struct planeStar {
  PetscScalar ij, ip1, im1, jp1, jm1;
};

struct planeBox {
  PetscScalar ij, ip1, im1, jp1, jm1, ip1jp1, im1jp1, ip1jm1, im1jm1;
};


//  NOTE: apply following abstraction to *age*, i.e. vtau, first

//! Class for a 3d DA-based Vec for ice scalar quantities in IceModel.
class IceModelVec3 : public IceModelVec {
public:
  IceModelVec3();
  virtual PetscErrorCode  destroy();
  virtual PetscErrorCode  create(IceGrid &mygrid, const char my_varname[]);

  virtual PetscErrorCode  needAccessToVals();
  virtual PetscErrorCode  doneAccessToVals();

public:
  // call needAccessToVals() before set???() or get???() *and* call doneAccessToVals() afterward
  PetscErrorCode  setValColumn(const PetscInt i, const PetscInt j, const PetscInt nlevels, 
                               PetscScalar *levelsIN, PetscScalar *valsIN);
  PetscErrorCode  setToConstantColumn(const PetscInt i, const PetscInt j,
                                      const PetscScalar c);

  PetscScalar     getValZ(const PetscInt i, const PetscInt j, const PetscScalar z);

  PetscErrorCode  getPlaneStarZ(const PetscInt i, const PetscInt j, const PetscScalar z,
                                planeStar *star);
  PetscErrorCode  getPlaneBoxZ(const PetscInt i, const PetscInt j, const PetscScalar z,
                               planeBox *box);

  PetscErrorCode  getValColumn(const PetscInt i, const PetscInt j, const PetscInt nlevels, 
                               PetscScalar *levelsIN, PetscScalar *valsOUT);

  PetscErrorCode  getHorSlice(Vec &gslice, const PetscScalar z);
  PetscErrorCode  getSurfaceValuesVec2d(Vec &gsurf, Vec myH);
  PetscErrorCode  getSurfaceValuesArray2d(Vec &gsurf, PetscScalar **H);

protected:
  PetscScalar     ***array;
  
  virtual PetscErrorCode checkHaveArray();
  PetscErrorCode         isLegalLevel(const PetscScalar z);

  PetscErrorCode  getInternalColumn(const PetscInt i, const PetscInt j, PetscInt *nlevels, 
                                    PetscScalar **levelsOUT, PetscScalar **valsOUT);
};


#endif /* __IceModelVec_hh */

