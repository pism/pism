// Copyright (C) 2012, 2014, 2015  David Maxwell
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

//! @note January 28, 2015: I don't think this code is used, which means
//! that it is almost certainly broken. -- CK

namespace pism {

template<class Problem, class VecArrayType>
class SNESDMCallbacks {
public:
  SNESDMCallbacks() {
    m_cbData.dm = NULL;
    m_cbData.p  = NULL;
  };

  void connect(SNES snes, Problem &p, DM dm, Vec r, Mat J, Mat Jpc=NULL) {
    PetscErrorCode ierr;
    if (m_cbData.dm != NULL) {
      throw RuntimeError::formatted("SNESDMCallbacks already connected");
    }
    m_cbData.dm = dm;
    m_cbData.p = &p;

    ierr = DMDASetLocalFunction(dm,
                                (DMDALocalFunction1)
                                SNESDMCallbacks<Problem,VecArrayType>::function_callback);
    PISM_CHK(ierr, "DMDASetLocalFunction");

    ierr = DMDASetLocalJacobian(dm,
                                (DMDALocalFunction1)
                                SNESDMCallbacks<Problem,VecArrayType>::jacobian_callback);
    PISM_CHK(ierr, "DMDASetLocalJacobian");

    ierr = SNESSetDM(snes, dm); PISM_CHK(ierr, "SNESSetDM");

    ierr = SNESSetFunction(snes, r, SNESDAFormFunction,
                           &m_cbData); PISM_CHK(ierr, "SNESSetFunction");

    if (Jpc == NULL) {
      Jpc = J;
    }

    ierr = SNESSetJacobian(snes, J, Jpc, SNESDAComputeJacobian,
                           &m_cbData); PISM_CHK(ierr, "SNESSetJacobian");
  }

protected:
  //! Adaptor for gluing SNESDAFormFunction callbacks to a C++ class
  /* The callbacks from SNES are mediated via SNESDAFormFunction, which has the
     convention that its context argument is a pointer to a struct
     having a DA as its first entry.  The SNESDMCallbackData fulfills
     this requirement, and allows for passing the callback on to a
     class. */
  struct CallbackData {
    DM dm;
    Problem  *p;
  };

  static PetscErrorCode function_callback(DMDALocalInfo *info,
                                          VecArrayType xg, VecArrayType yg,
                                          CallbackData *data) {
    PetscErrorCode ierr;
    ierr = data->p->assembleFunction(info,xg,yg); CHKERRQ(ierr);
    return 0;
  }

  static PetscErrorCode jacobian_callback(DMDALocalInfo *info,
                                          VecArrayType xg, Mat J,
                                          CallbackData *data) {
    PetscErrorCode ierr;
    ierr = data->p->assembleJacobian(info,xg,J); CHKERRQ(ierr);
    return 0;
  }

  CallbackData m_cbData;
};

template<class Problem>
class SNESCallbacks {
public:

  static void connect(SNES snes, Problem &p, Vec r, Mat J, Mat Jpc=NULL) {
    PetscErrorCode ierr;

    ierr = SNESSetFunction(snes, r, SNESCallbacks<Problem>::function_callback,
                           &p); PISM_CHK(ierr, "SNESSetFunction");

    if (Jpc == NULL) {
      Jpc = J;
    }

    ierr = SNESSetJacobian(snes, J, Jpc, SNESCallbacks<Problem>::jacobian_callback,
                           &p); PISM_CHK(ierr, "SNESSetJacobian");
  }

protected:
  static PetscErrorCode function_callback(SNES snes, Vec x, Vec f, void*ctx) {
    Problem *p = reinterpret_cast<Problem *>(ctx);
    PetscErrorCode ierr = p->assembleFunction(x,f); CHKERRQ(ierr);
    return 0;
  }

  static PetscErrorCode jacobian_callback(SNES snes,
                                          Vec x, Mat *J, Mat*Jpc, MatStructure* structure,
                                          void *ctx) {
    Problem *p = reinterpret_cast<Problem *>(ctx);
    PetscErrorCode ierr = p->assembleJacobian(x,J); CHKERRQ(ierr);
    return 0;
  }
};

} // end of namespace pism

#endif /* end of include guard: SNESPROBLEM_HH_2EZQQ4UH */
