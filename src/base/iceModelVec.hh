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

  virtual PetscErrorCode  create(IceGrid &mygrid, const char my_varname[], bool local);
  virtual PetscErrorCode  destroy();

  virtual PetscErrorCode  setVaridNC(const int my_varid);
  virtual PetscErrorCode  setAttrsNC(const int my_varid,
             const char my_long_name[], const char my_units[], const char my_pism_intent[]);
  virtual PetscErrorCode  setAttrsCFstandardNC(const int my_varid,
             const char my_long_name[], const char my_units[], const char my_pism_intent[],
             const char my_standard_name[]);

  virtual PetscErrorCode  findVecNC(const int ncid, PetscTruth *exists);  //FIXME: not implemented

  virtual PetscErrorCode  getVecNC(const int ncid, const int *s, const int *c, int dims, 
                                   void *a_mpi, int a_size);

  virtual PetscErrorCode  writeAttrsNC(const int ncid);  //FIXME: not implemented

  virtual PetscErrorCode  putVecNC(const int ncid, const int *s, const int *c, int dims, 
                                   void *a_mpi, int a_size);

  virtual PetscErrorCode  regridVecNC(const char *vars, char c, int dim_flag, LocalInterpCtx &lic);

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
  PetscTruth   has_standard_name;
  
  IceGrid      *grid;
  DA           da;
  bool         localp, IOwnDA;
  
  int          varid_nc;

  void         *array;  // will be PetscScalar** or PetscScalar*** in derived classes

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


// Class for a 2d DA-based Vec, with STAR stencil, for ice scalar quantities in IceModel.
class IceModelVec2 : public IceModelVec {
public:
  IceModelVec2();
  virtual PetscErrorCode  create(IceGrid &my_grid, const char my_varname[], bool local);

  PetscScalar     getVal(const PetscInt i, const PetscInt j);
  PetscErrorCode  getPlaneStar(const PetscInt i, const PetscInt j, planeStar *star);

  PetscScalar**   arrayGet();  // FIXME: potentially dangerous

protected:
  PetscErrorCode  create(IceGrid &my_grid, const char my_varname[], bool local, DAStencilType my_sten);
};


// Class for a 2d DA-based Vec, with BOX stencil, for ice scalar quantities in IceModel.
class IceModelVec2Box : public IceModelVec2 {
public:
  IceModelVec2Box();
  virtual PetscErrorCode  create(IceGrid &mygrid, const char my_varname[], bool local);

  PetscErrorCode  getPlaneBox(const PetscInt i, const PetscInt j, planeBox *box);
};


//! Class for a 3d DA-based Vec for bedrock scalar quantities in IceModel.
class IceModelVec3Bedrock : public IceModelVec {
public:
  IceModelVec3Bedrock();
  virtual PetscErrorCode  create(IceGrid &mygrid, const char my_varname[], bool local);

  PetscErrorCode  setInternalColumn(const PetscInt i, const PetscInt j, PetscScalar *valsIN);
  PetscErrorCode  setToConstantColumn(const PetscInt i, const PetscInt j,
                                      const PetscScalar c);
  PetscErrorCode  getInternalColumn(const PetscInt i, const PetscInt j, PetscScalar **valsOUT);

  PetscErrorCode  setValColumn(const PetscInt i, const PetscInt j, const PetscInt nlevels, 
                               PetscScalar *levelsIN, PetscScalar *valsIN);
  PetscErrorCode  getValColumn(const PetscInt i, const PetscInt j, const PetscInt nlevels, 
                               PetscScalar *levelsIN, PetscScalar *valsOUT);

protected:  
  PetscErrorCode  isLegalLevel(const PetscScalar z);
};


// see iceModelVec3.cc for implementation of next class:

//! Class for a 3d DA-based Vec for ice scalar quantities in IceModel.
class IceModelVec3 : public IceModelVec {
public:
  IceModelVec3();
  virtual PetscErrorCode  create(IceGrid &mygrid, const char my_varname[], bool local);
  PetscErrorCode          createSameDA(IceModelVec3 imv3_dasource, 
                                       IceGrid &mygrid, const char my_varname[], bool local);

  // note the IceModelVec3 with this method must be *local* while imv3_source must be *global*
  virtual PetscErrorCode  beginGhostCommTransfer(IceModelVec3 imv3_source);
  virtual PetscErrorCode  endGhostCommTransfer(IceModelVec3 imv3_source);

  // call needAccessToVals() before set???() or get???() *and* call doneAccessToVals() afterward
  PetscErrorCode  setValColumn(const PetscInt i, const PetscInt j, const PetscInt nlevels, 
                               PetscScalar *levelsIN, PetscScalar *valsIN);
  PetscErrorCode  setToConstantColumn(const PetscInt i, const PetscInt j, const PetscScalar c);

  PetscScalar     getValZ(const PetscInt i, const PetscInt j, const PetscScalar z);

  PetscErrorCode  getPlaneStarZ(const PetscInt i, const PetscInt j, const PetscScalar z,
                                planeStar *star);

//Idea; see IceModel::horizontalVelocitySIARegular(), for example; would only work for 4 neighbors
//   if used with two IceModelVec3:
//  PetscErrorCode  averageStagNeighCol(const PetscInt i, const PetscInt j, 
//                                      const PetscInt nlevels, PetscScalar *levelsIN, 
//                                      PetscScalar *avOUT);

  PetscErrorCode  getValColumn(const PetscInt i, const PetscInt j, const PetscInt nlevels, 
                               PetscScalar *levelsIN, PetscScalar *valsOUT);

  PetscErrorCode  getHorSlice(Vec &gslice, const PetscScalar z);
  PetscErrorCode  getSurfaceValuesVec2d(Vec &gsurf, Vec myH);
  PetscErrorCode  getSurfaceValuesArray2d(Vec &gsurf, PetscScalar **H);

protected:  
  PetscErrorCode  isLegalLevel(const PetscScalar z);
  PetscErrorCode  getInternalColumn(const PetscInt i, const PetscInt j, PetscInt *nlevels, 
                                    PetscScalar **levelsOUT, PetscScalar **valsOUT);
};

#endif /* __IceModelVec_hh */

