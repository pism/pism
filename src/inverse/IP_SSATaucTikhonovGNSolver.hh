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

#ifndef IP_SSATAUCTIKHONOVGN_HH_SIU7F33G
#define IP_SSATAUCTIKHONOVGN_HH_SIU7F33G

#include "iceModelVec.hh"
#include "IP_SSATaucForwardProblem.hh"
#include "functional/IPFunctional.hh"
#include "TerminationReason.hh"

namespace pism {

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

class IP_SSATaucTikhonovGNSolver {
public:
  typedef IceModelVec2S DesignVec;
  typedef IceModelVec2V StateVec;
  // typedef IP_SSATaucTikhonovGNSolverListener Listener;

  IP_SSATaucTikhonovGNSolver( IP_SSATaucForwardProblem &ssaforward, DesignVec &d0, StateVec &u_obs, double eta, 
                              IPInnerProductFunctional<DesignVec> &designFunctional, IPInnerProductFunctional<StateVec> &stateFunctional);

  ~IP_SSATaucTikhonovGNSolver();
  
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

  //! Sets the desired target misfit (in units of \f$\sqrt{J_{\rm misfit}}\f$).
  virtual PetscErrorCode setTargetMisfit( double misfit) {
    m_target_misfit = misfit;
    return 0;
  }

  virtual PetscErrorCode evaluateGNFunctional(DesignVec &h, double *value);

  virtual PetscErrorCode apply_GN(IceModelVec2S &h, IceModelVec2S &out);
  virtual PetscErrorCode apply_GN(Vec h, Vec out);

  virtual PetscErrorCode init(TerminationReason::Ptr &reason);

  virtual PetscErrorCode check_convergence(TerminationReason::Ptr &reason); 
  
  virtual PetscErrorCode solve(TerminationReason::Ptr &reason);

  virtual PetscErrorCode evaluate_objective_and_gradient(TerminationReason::Ptr &reason);

protected:

  virtual PetscErrorCode assemble_GN_rhs(DesignVec &out);

  virtual PetscErrorCode solve_linearized(TerminationReason::Ptr &reason);

  virtual PetscErrorCode compute_dlogalpha(double *dalpha, TerminationReason::Ptr &reason);

  virtual PetscErrorCode linesearch(TerminationReason::Ptr &reason);

  PetscErrorCode construct();
  PetscErrorCode destruct();

  IP_SSATaucForwardProblem &m_ssaforward;

  DesignVec m_x;
  DesignVec m_y;

  DesignVec m_tmp_D1Global;
  DesignVec m_tmp_D2Global;
  DesignVec m_tmp_D1Local;
  DesignVec m_tmp_D2Local;
  StateVec  m_tmp_S1Global;
  StateVec  m_tmp_S2Global;
  StateVec  m_tmp_S1Local;
  StateVec  m_tmp_S2Local;

  DesignVec  m_GN_rhs;

  DesignVec &m_d0;
  DesignVec m_dGlobal;
  DesignVec m_d;
  DesignVec m_d_diff;
  DesignVec m_d_diff_lin;

  DesignVec m_h;
  DesignVec m_hGlobal;
  DesignVec m_dalpha_rhs;
  DesignVec m_dh_dalpha;
  DesignVec m_dh_dalphaGlobal;

  DesignVec m_grad_design;
  DesignVec m_grad_state;
  DesignVec m_gradient;
  
  double m_val_design, m_val_state, m_value;

  StateVec &m_u_obs;
  StateVec m_u_diff;

  KSP m_ksp;  
  Mat m_mat_GN;

  double m_eta;
  IPInnerProductFunctional<DesignVec> &m_designFunctional;
  IPInnerProductFunctional<StateVec> &m_stateFunctional;

  double m_alpha;
  double m_logalpha;
  double m_target_misfit;

  int m_iter, m_iter_max;
  bool m_tikhonov_adaptive;
  double m_vel_scale;
  double m_tikhonov_rtol, m_tikhonov_atol, m_tikhonov_ptol;

  MPI_Comm m_comm;

};

} // end of namespace pism

#endif /* end of include guard: IP_SSATAUCTIKHONOVGN_HH_SIU7F33G */
