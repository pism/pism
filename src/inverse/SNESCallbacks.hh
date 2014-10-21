// Copyright (C) 2012, 2014  David Maxwell
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#ifndef SNESPROBLEM_HH_2EZQQ4UH
#define SNESPROBLEM_HH_2EZQQ4UH

#include <petscsnes.h>
#include <petscdmda.h>

namespace pism {

template<class Problem, class VecArrayType>
class SNESDMCallbacks {
public:
  SNESDMCallbacks() {
    m_cbData.dm = NULL;
    m_cbData.p = NULL;
  };
  
  PetscErrorCode connect(SNES snes, Problem &p, DM dm, Vec r, Mat J, Mat Jpc=NULL) {
    PetscErrorCode ierr;
    if (m_cbData.dm != NULL) {
      throw RuntimeError::formatted("SNESDMCallbacks already connected");
    }
    m_cbData.dm = dm;
    m_cbData.p = &p;
    
    ierr = DMDASetLocalFunction(dm,(DMDALocalFunction1) SNESDMCallbacks<Problem,VecArrayType>::formFunctionCallback); CHKERRQ(ierr);
    ierr = DMDASetLocalJacobian(dm,(DMDALocalFunction1) SNESDMCallbacks<Problem,VecArrayType>::formJacobianCallback); CHKERRQ(ierr);
    ierr = SNESSetDM(snes, dm); CHKERRQ(ierr);
    ierr = SNESSetFunction(snes, r, SNESDAFormFunction, &m_cbData); CHKERRQ(ierr);
    if (Jpc==NULL) Jpc = J;
    ierr = SNESSetJacobian(snes, J, Jpc, SNESDAComputeJacobian, &m_cbData); CHKERRQ(ierr);
    return 0;
  }

protected:
  //! Adaptor for gluing SNESDAFormFunction callbacks to a C++ class
  /* The callbacks from SNES are mediated via SNESDAFormFunction, which has the
     convention that its context argument is a pointer to a struct 
     having a DA as its first entry.  The SNESDMCallbackData fulfills 
     this requirement, and allows for passing the callback on to a
     class. */
  struct SNESDMCallbackData {
    DM           dm;
    Problem  *p;
  };
  static PetscErrorCode formFunctionCallback(DMDALocalInfo *info,
                                             VecArrayType xg, VecArrayType yg,
                                             SNESDMCallbackData *data) {
    PetscErrorCode ierr;
    ierr = data->p->assembleFunction(info,xg,yg); CHKERRQ(ierr);
    return 0;
  }

  static PetscErrorCode formJacobianCallback(DMDALocalInfo *info,
                                             VecArrayType xg, Mat J,
                                             SNESDMCallbackData *data) {
    PetscErrorCode ierr;
    ierr = data->p->assembleJacobian(info,xg,J); CHKERRQ(ierr);
    return 0;
  }

  SNESDMCallbackData m_cbData;
};

template<class Problem>
class SNESCallbacks {
public:
  
  static PetscErrorCode connect (SNES snes, Problem &p, Vec r, Mat J, Mat Jpc=NULL) {
    PetscErrorCode ierr;
    
    ierr = SNESSetFunction(snes,r,SNESCallbacks<Problem>::formFunctionCallback,&p); CHKERRQ(ierr);
    if (Jpc==NULL) Jpc = J;
    ierr = SNESSetJacobian(snes,J,Jpc,SNESCallbacks<Problem>::formJacobianCallback,&p); CHKERRQ(ierr);
    return 0;
  }

protected:
  static PetscErrorCode formFunctionCallback(SNES snes,Vec x,Vec f, void*ctx) {
    PetscErrorCode ierr;
    Problem *p = reinterpret_cast<Problem *>(ctx);
    ierr = p->assembleFunction(x,f); CHKERRQ(ierr);
    return 0;
  }

  static PetscErrorCode formJacobianCallback(SNES snes,
                                             Vec x, Mat *J, Mat*Jpc, MatStructure* structure,
                                             void *ctx) {
    PetscErrorCode ierr;
    Problem *p = reinterpret_cast<Problem *>(ctx);
    ierr = p->assembleJacobian(x,J); CHKERRQ(ierr);
    return 0;
  }
};

} // end of namespace pism

#endif /* end of include guard: SNESPROBLEM_HH_2EZQQ4UH */
