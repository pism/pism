// Copyright (C) 2009 Ed Bueler
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

//! Virtual base class.  Describes a tridiagonal system to solve in a column of ice.
/*!
Because both the age evolution and conservation of energy equations require us to set up
and solve a tridiagonal system of equations, this is structure is worth abstracting.

This base class just holds the tridiagonal system and the ability to
solve it, but does not insert entries into the relevant matrix locations.
Derived classes will actually set up instances of the system.

The sequence requires setting the column-independent (public) data members,
calling the initAllColumns() routine, and then setting up and solving
the system in each column.
 */
class columnSystemCtx {

public:
  columnSystemCtx(int my_nmax); /*! allocate a tridiagonal system
                                 of maximum size nmax */
  ~columnSystemCtx();           //! deallocate it

  PetscErrorCode viewColumnValues(PetscViewer viewer, 
                                  PetscScalar *v, PetscInt m, const char* info) const;

protected:
  PetscInt    nmax;
  PetscScalar *L, *Lp, *D, *U, *rhs, *work; // vectors for tridiagonal system

  // deliberately protected so only derived classes can use
  PetscErrorCode solveTridiagonalSystem(PetscInt n, PetscScalar **x);
};


//! Tridiagonal linear system for vertical column of age (pure advection) problem.
/*!
Call sequence like this:
\code
  ageSystemCtx foo;
  foo.dx = ...  // set public constants
  foo.u = ...   // set public pointers
  foo.initAllColumns();
  for (i in ownership) {
    for (j in ownership) {
      ks = ...
      foo.setIndicesThisColumn(i,j,ks);
      foo.solveThisColumn(x);
    }  
  }
\endcode
 */
class ageSystemCtx : public columnSystemCtx {

public:
  ageSystemCtx(PetscInt my_Mz);
  PetscErrorCode initAllColumns();
  PetscErrorCode setIndicesThisColumn(PetscInt my_i, PetscInt my_j, PetscInt my_ks);  
  PetscErrorCode solveThisColumn(PetscScalar **x);  

public:
  // constants which should be set before calling initForAllColumns()
  PetscScalar  dx,
               dy,
               dtAge,
               dzEQ;
  // pointers which should be set before calling initForAllColumns()
  PetscScalar  *u,
               *v,
               *w;
  IceModelVec3 *tau3;

protected: // used internally
  PetscInt    i, j, ks;
  PetscScalar nuEQ;
  bool        initAllDone,
              indicesValid;
};


//! Tridiagonal linear system for vertical column of temperature-based conservation of energy problem.
/*!
Call sequence like this:
\code
  tempSystemCtx foo;
  foo.dx = ...  // set public constants
  foo.u = ...   // set public pointers
  foo.initAllColumns();
  for (i in ownership) {
    for (j in ownership) {
      ks = ...
      foo.setIndicesThisColumn(i,j,ks);
      [COMPUTE OTHER PARAMS]
      foo.setSchemeParamsThisColumn(mask,isMarginal,lambda);  
      foo.setSurfaceBoundaryValuesThisColumn(Ts);
      foo.setBasalBoundaryValuesThisColumn(Ghf,Tshelfbase,Rb);
      foo.solveThisColumn(x);
    }  
  }
\endcode
 */
class tempSystemCtx : public columnSystemCtx {

public:
  tempSystemCtx(PetscInt my_Mz, PetscInt my_Mbz);
  PetscErrorCode initAllColumns();
  PetscErrorCode setIndicesThisColumn(PetscInt i, PetscInt j, PetscInt ks);  
  PetscErrorCode setSchemeParamsThisColumn(
                     PetscScalar my_mask, bool my_isMarginal, PetscScalar my_lambda);  
  PetscErrorCode setSurfaceBoundaryValuesThisColumn(PetscScalar my_Ts);
  PetscErrorCode setBasalBoundaryValuesThisColumn(
                     PetscScalar my_Ghf, PetscScalar my_Tshelfbase, PetscScalar my_Rb);
  PetscErrorCode solveThisColumn(PetscScalar **x);  

public:
  // constants which should be set before calling initForAllColumns()
  PetscScalar  dx,
               dy,
               dtTemp,
               dzEQ,
               dzbEQ,
               ice_rho,
               ice_c_p,
               ice_k,
               bed_thermal_rho,
               bed_thermal_c_p,
               bed_thermal_k;
  // pointers which should be set before calling initForAllColumns()
  PetscScalar  *T,
               *Tb,
               *u,
               *v,
               *w,
               *Sigma;
  IceModelVec3 *T3;

protected: // used internally
  PetscInt    Mz, Mbz, k0;
  PetscInt    i, j, ks;
  PetscScalar mask, lambda, Ts, Ghf, Tshelfbase, Rb;
  bool        isMarginal;
  PetscScalar nuEQ,
              rho_c_I,
              rho_c_br,
              rho_c_av,
              iceK,
              iceR,
              brK,
              brR,
              rho_c_ratio,
              dzav,
              iceReff,
              brReff;
  bool        initAllDone,
              indicesValid,
              schemeParamsValid,
              surfBCsValid,
              basalBCsValid;
};

#endif	/* __columnSystem_hh */

