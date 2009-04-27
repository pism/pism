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
  PetscScalar *L, *Lp, *D, *U, *rhs, *work; // vectors for tridiagonal system

  // deliberately protected so only derived classes can use
  PetscErrorCode solveTridiagonalSystem(PetscInt n, PetscScalar **x);
};


//! Tridiagonal linear system for vertical column in solving age advection problem.
class ageSystemCtx : public columnSystemCtx {

public:
  ageSystemCtx(int my_Mz);
  // arguments do not depend on which column
  PetscErrorCode ageSetConstants(
    PetscScalar dx, PetscScalar dy, PetscScalar dtTempAge, PetscScalar dzEQ);
  // arguments depend on which column
  PetscErrorCode ageColumnSetUpAndSolve(
    PetscInt i, PetscInt j, PetscInt ks,
    PetscScalar *u, PetscScalar *v, PetscScalar *w, IceModelVec3 &tau3,
    PetscScalar **x);  

protected:
  PetscInt callcount;
  PetscScalar dx, dy, dtTempAge, dzEQ, nuEQ;
};


//! Tridiagonal linear system for vertical column in solving temperature-based conservation of energy problem.
class tempSystemCtx : public columnSystemCtx {

public:
  tempSystemCtx(int my_Mz, int my_Mbz);
  // arguments do not depend on which column
  PetscErrorCode tempSetConstants(
    PetscScalar dx, PetscScalar dy, PetscScalar dtTempAge,
    PetscScalar dzEQ, PetscScalar dzbEQ,
    PetscScalar ice_rho, PetscScalar ice_c_p, PetscScalar ice_k,
    PetscScalar bed_thermal_rho, PetscScalar bed_thermal_c_p, PetscScalar bed_thermal_k);
  // arguments depend on which column
  PetscErrorCode tempColumnSetUpAndSolve(
    PetscInt i, PetscInt j,
    PetscInt ks, bool isMarginal, PetscScalar lambda,
    PetscScalar *T, PetscScalar *Tb,
    PetscScalar *u, PetscScalar *v, PetscScalar *w, PetscScalar *Sigma,
    PetscScalar Ghf_ij, PetscScalar Ts_ij, PetscScalar mask_ij,
    PetscScalar Tshelfbase_ij, PetscScalar Rb_ij,
    IceModelVec3 &T3,
    PetscScalar **x);  

protected:
  PetscInt    callcount, Mz, Mbz, k0;
  PetscScalar dx, dy, dtTempAge,
              dzEQ, dzbEQ,
              ice_rho, ice_c_p, ice_k,
              bed_thermal_rho, bed_thermal_c_p, bed_thermal_k,
              nuEQ, rho_c_I, rho_c_br, rho_c_av,
              iceK, iceR, brK, brR,
              rho_c_ratio, dzav, iceReff, brReff;
};

#endif	/* __columnSystem_hh */

