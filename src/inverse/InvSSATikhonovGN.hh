// Copyright (C) 2012  David Maxwell
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

#ifndef INVSSATIKHONOVGN_HH_SIU7F33G
#define INVSSATIKHONOVGN_HH_SIU7F33G

#include "iceModelVec.hh"
#include "InvSSAForwardProblem.hh"
#include "Functional.hh"
#include "TerminationReason.hh"

template<class C,PetscErrorCode (C::*MultiplyCallback)(Vec,Vec) >
class MatrixMultiplyCallback {
public:
  static PetscErrorCode connect(Mat A) {
    PetscErrorCode ierr;
    ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void))MatrixMultiplyCallback::multiply); CHKERRQ(ierr); 
    return 0;
  }
protected:
  static PetscErrorCode multiply(Mat A, Vec x, Vec y) {
    PetscErrorCode ierr;
    C *ctx;
    ierr = MatShellGetContext(A,&ctx); CHKERRQ(ierr);
    ierr = (ctx->*MultiplyCallback)(x,y); CHKERRQ(ierr);
    return 0;
  }
};

class InvSSATikhonovGN {
public:
  typedef IceModelVec2S DesignVec;
  typedef IceModelVec2V StateVec;
  // typedef InvSSATikhonovGNListener Listener;

  InvSSATikhonovGN( InvSSAForwardProblem &ssaforward, DesignVec &d0, StateVec &u_obs, PetscReal eta, 
                    IPFunctional<DesignVec> &designFunctional, IPFunctional<StateVec> &stateFunctional);

  ~InvSSATikhonovGN();
  
  virtual StateVec &stateSolution() {
    return m_ssaforward.solution();
  }

  virtual DesignVec &designSolution() {
    return m_d;
  }

  virtual PetscErrorCode setInitialGuess( DesignVec &d) {
    PetscErrorCode ierr;
    ierr = m_d.copy_from(d); CHKERRQ(ierr);
    return 0;
  }

  virtual PetscErrorCode apply_GN(IceModelVec2S &h, IceModelVec2S &out);
  virtual PetscErrorCode apply_GN(Vec h, Vec out);

  virtual PetscErrorCode assemble_GN_rhs(DesignVec &out);

  virtual PetscErrorCode init(TerminationReason::Ptr &reason);
  virtual PetscErrorCode solve(TerminationReason::Ptr &reason);
  virtual PetscErrorCode solve_linearized(TerminationReason::Ptr &reason);

protected:

  PetscErrorCode construct();
  PetscErrorCode destruct();

  InvSSAForwardProblem &m_ssaforward;

  DesignVec m_x;
  DesignVec m_y;

  DesignVec m_tmp_D1;
  StateVec  m_tmp_S1;
  StateVec  m_tmp_S2;

  DesignVec  m_GN_rhs;

  DesignVec &m_d0;
  DesignVec m_dGlobal;
  DesignVec m_d;
  DesignVec m_d_diff;

  StateVec &m_u_obs;
  StateVec m_u_diff;

  KSP m_ksp;  
  Mat m_mat_GN;

  PetscReal m_eta;
  IPFunctional<DesignVec> &m_designFunctional;
  IPFunctional<StateVec> &m_stateFunctional;

  MPI_Comm m_comm;

};



#endif /* end of include guard: INVSSATIKHONOVGN_HH_SIU7F33G */
