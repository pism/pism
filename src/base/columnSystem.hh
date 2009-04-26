// Copyright (C) 2004-2009 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __columnSystem_hh
#define __columnSystem_hh

#include <petsc.h>
#include "iceModelVec.hh"

//! Virtual class for a tridiagonal system to solve in a column of ice.
/*!
Because both IceModel::ageStep() and IceModel::temperatureStep() set up
and solve a tridiagonal system of equations, this is worth abstracting.

This base class just holds the tridiagonal system and the ability to
solve it.  Derived classes will actually set up instances of the system.
 */
class columnSystemCtx {

public:
  columnSystemCtx(int my_nmax); /*! allocate a tridiagonal system
                                 of maximum size nmax */
  ~columnSystemCtx();           //! deallocate it

protected:
  PetscInt    nmax;
  PetscScalar *L,
              *Lp,
              *D,
              *U,
              *rhs,
              *work;

  // deliberately protected so only derived classes can use
  PetscErrorCode solveTridiagonalSystem(
    PetscInt n,    // size of instance; n <= nmax
    PetscScalar **x);
};


class ageSystemCtx : public columnSystemCtx {

public:
  ageSystemCtx(int my_Mz);
  PetscErrorCode ageSetConstants(
    PetscScalar dx, PetscScalar dy, PetscScalar dtTempAge, PetscScalar nuEQ, 
    PetscScalar *zlevEQ);
  PetscErrorCode ageColumnWiseSetUpAndSolve(
    PetscInt i, PetscInt j, PetscInt ks,
    PetscScalar *u, PetscScalar *v, PetscScalar *w,
    IceModelVec3 &tau3,
    PetscScalar **x);  

protected:
  PetscInt Mz;
  PetscScalar nuEQ, dx, dy, dtTempAge;
  PetscScalar *zlevEQ;
  IceModelVec3 *tau3; 
  
private:
  PetscInt callcount;
};

#endif	/* __columnSystem_hh */


