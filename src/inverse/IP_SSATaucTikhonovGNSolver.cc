// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017  David Maxwell and Constantine Khroulev
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

#include "IP_SSATaucTikhonovGNSolver.hh"
#include "pism/util/TerminationReason.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"

namespace pism {
namespace inverse {

IP_SSATaucTikhonovGNSolver::IP_SSATaucTikhonovGNSolver(IP_SSATaucForwardProblem &ssaforward,
                                                       DesignVec &d0, StateVec &u_obs, double eta,
                                                       IPInnerProductFunctional<DesignVec> &designFunctional,
                                                       IPInnerProductFunctional<StateVec> &stateFunctional)
  : m_ssaforward(ssaforward), m_d0(d0), m_u_obs(u_obs), m_eta(eta),
    m_designFunctional(designFunctional), m_stateFunctional(stateFunctional),
    m_target_misfit(0.0) {
  PetscErrorCode ierr;
  IceGrid::ConstPtr grid = m_d0.grid();
  m_comm = grid->com;

  unsigned int design_stencil_width = m_d0.stencil_width();
  unsigned int state_stencil_width = m_u_obs.stencil_width();

  m_x.create(grid, "x", WITH_GHOSTS, design_stencil_width);

  m_tmp_D1Global.create(grid, "work vector", WITHOUT_GHOSTS, 0);
  m_tmp_D2Global.create(grid, "work vector", WITHOUT_GHOSTS, 0);
  m_tmp_S1Global.create(grid, "work vector", WITHOUT_GHOSTS, 0);
  m_tmp_S2Global.create(grid, "work vector", WITHOUT_GHOSTS, 0);

  m_tmp_D1Local.create(grid, "work vector", WITH_GHOSTS, design_stencil_width);
  m_tmp_D2Local.create(grid, "work vector", WITH_GHOSTS, design_stencil_width);
  m_tmp_S1Local.create(grid, "work vector", WITH_GHOSTS, state_stencil_width);
  m_tmp_S2Local.create(grid, "work vector", WITH_GHOSTS, state_stencil_width);

  m_GN_rhs.create(grid, "GN_rhs", WITHOUT_GHOSTS, 0);

  m_dGlobal.create(grid, "d (sans ghosts)", WITHOUT_GHOSTS, 0);

  m_d.reset(new DesignVec(grid, "d", WITH_GHOSTS, design_stencil_width));
  m_d_diff.create(grid, "d_diff", WITH_GHOSTS, design_stencil_width);
  m_d_diff_lin.create(grid, "d_diff linearized", WITH_GHOSTS, design_stencil_width);
  m_h.create(grid, "h", WITH_GHOSTS, design_stencil_width);
  m_hGlobal.create(grid, "h (sans ghosts)", WITHOUT_GHOSTS);

  m_dalpha_rhs.create(grid, "dalpha rhs", WITHOUT_GHOSTS);
  m_dh_dalpha.create(grid, "dh_dalpha", WITH_GHOSTS, design_stencil_width);
  m_dh_dalphaGlobal.create(grid, "dh_dalpha", WITHOUT_GHOSTS);
  m_u_diff.create(grid, "du", WITH_GHOSTS, state_stencil_width);

  m_grad_design.create(grid, "grad design", WITHOUT_GHOSTS);
  m_grad_state.create(grid, "grad design", WITHOUT_GHOSTS);
  m_gradient.create(grid, "grad design", WITHOUT_GHOSTS);

  ierr = KSPCreate(grid->com, m_ksp.rawptr());
  PISM_CHK(ierr, "KSPCreate");

  ierr = KSPSetOptionsPrefix(m_ksp, "inv_gn_");
  PISM_CHK(ierr, "KSPSetOptionsPrefix");

  double ksp_rtol = 1e-5; // Soft tolerance
  ierr = KSPSetTolerances(m_ksp, ksp_rtol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
  PISM_CHK(ierr, "KSPSetTolerances");

  ierr = KSPSetType(m_ksp, KSPCG);
  PISM_CHK(ierr, "KSPSetType");

  PC pc;
  ierr = KSPGetPC(m_ksp, &pc);
  PISM_CHK(ierr, "KSPGetPC");

  ierr = PCSetType(pc, PCNONE);
  PISM_CHK(ierr, "PCSetType");

  ierr = KSPSetFromOptions(m_ksp);
  PISM_CHK(ierr, "KSPSetFromOptions");

  int nLocalNodes  = grid->xm()*grid->ym();
  int nGlobalNodes = grid->Mx()*grid->My();
  ierr = MatCreateShell(grid->com, nLocalNodes, nLocalNodes,
                        nGlobalNodes, nGlobalNodes, this, m_mat_GN.rawptr());
  PISM_CHK(ierr, "MatCreateShell");

  typedef MatrixMultiplyCallback<IP_SSATaucTikhonovGNSolver,
                                 &IP_SSATaucTikhonovGNSolver::apply_GN> multCallback;
  multCallback::connect(m_mat_GN);

  m_alpha = 1./m_eta;
  m_logalpha = log(m_alpha);

  m_tikhonov_adaptive = options::Bool("-tikhonov_adaptive", "Tikhonov adaptive");

  m_iter_max = 1000;
  m_iter_max = options::Integer("-inv_gn_iter_max", "", m_iter_max);

  m_tikhonov_atol = grid->ctx()->config()->get_double("inverse.tikhonov.atol");
  m_tikhonov_rtol = grid->ctx()->config()->get_double("inverse.tikhonov.rtol");
  m_tikhonov_ptol = grid->ctx()->config()->get_double("inverse.tikhonov.ptol");

  m_log = d0.grid()->ctx()->log();
}

IP_SSATaucTikhonovGNSolver::~IP_SSATaucTikhonovGNSolver() {
  // empty
}

TerminationReason::Ptr IP_SSATaucTikhonovGNSolver::init() {
  return m_ssaforward.linearize_at(m_d0);
}

void IP_SSATaucTikhonovGNSolver::apply_GN(IceModelVec2S &x, IceModelVec2S &y) {
  this->apply_GN(x.vec(), y.vec());
}

//! @note This function has to return PetscErrorCode (it is used as a callback).
void  IP_SSATaucTikhonovGNSolver::apply_GN(Vec x, Vec y) {
  StateVec  &tmp_gS = m_tmp_S1Global;
  StateVec  &Tx     = m_tmp_S1Local;
  DesignVec &tmp_gD = m_tmp_D1Global;
  DesignVec &GNx    = m_tmp_D2Global;

  // FIXME: Needless copies for now.
  m_x.copy_from_vec(x);

  m_ssaforward.apply_linearization(m_x,Tx);
  Tx.update_ghosts();

  m_stateFunctional.interior_product(Tx,tmp_gS);

  m_ssaforward.apply_linearization_transpose(tmp_gS,GNx);

  m_designFunctional.interior_product(m_x,tmp_gD);
  GNx.add(m_alpha,tmp_gD);

  PetscErrorCode ierr = VecCopy(GNx.vec(), y); PISM_CHK(ierr, "VecCopy");
}

void IP_SSATaucTikhonovGNSolver::assemble_GN_rhs(DesignVec &rhs) {

  rhs.set(0);
  
  m_stateFunctional.interior_product(m_u_diff,m_tmp_S1Global);
  m_ssaforward.apply_linearization_transpose(m_tmp_S1Global,rhs);

  m_designFunctional.interior_product(m_d_diff,m_tmp_D1Global);
  rhs.add(m_alpha,m_tmp_D1Global);
  
  rhs.scale(-1);
}

TerminationReason::Ptr IP_SSATaucTikhonovGNSolver::solve_linearized() {
  PetscErrorCode ierr;

  this->assemble_GN_rhs(m_GN_rhs);

  ierr = KSPSetOperators(m_ksp,m_mat_GN,m_mat_GN);
  PISM_CHK(ierr, "KSPSetOperators");

  ierr = KSPSolve(m_ksp,m_GN_rhs.vec(),m_hGlobal.vec());
  PISM_CHK(ierr, "KSPSolve");

  KSPConvergedReason ksp_reason;
  ierr = KSPGetConvergedReason(m_ksp ,&ksp_reason);
  PISM_CHK(ierr, "KSPGetConvergedReason");
  
  m_h.copy_from(m_hGlobal);

  return TerminationReason::Ptr(new KSPTerminationReason(ksp_reason));
}

void IP_SSATaucTikhonovGNSolver::evaluateGNFunctional(DesignVec &h, double *value) {
  
  m_ssaforward.apply_linearization(h,m_tmp_S1Local);
  m_tmp_S1Local.update_ghosts();
  m_tmp_S1Local.add(1,m_u_diff);
  
  double sValue;
  m_stateFunctional.valueAt(m_tmp_S1Local,&sValue);
  
  m_tmp_D1Local.copy_from(m_d_diff);
  m_tmp_D1Local.add(1,h);
  
  double dValue;
  m_designFunctional.valueAt(m_tmp_D1Local,&dValue);
  
  *value = m_alpha*dValue + sValue;
}


TerminationReason::Ptr IP_SSATaucTikhonovGNSolver::check_convergence() {

  double designNorm, stateNorm, sumNorm;
  double dWeight, sWeight;
  dWeight = m_alpha;
  sWeight = 1;

  designNorm = m_grad_design.norm(NORM_2);
  stateNorm  = m_grad_state.norm(NORM_2);

  designNorm *= dWeight;
  stateNorm  *= sWeight;

  sumNorm = m_gradient.norm(NORM_2);

  m_log->message(2,
             "----------------------------------------------------------\n");
  m_log->message(2,
             "IP_SSATaucTikhonovGNSolver Iteration %d: misfit %g; functional %g \n",
             m_iter, sqrt(m_val_state)*m_vel_scale, m_value*m_vel_scale*m_vel_scale);
  if (m_tikhonov_adaptive) {
    m_log->message(2, "alpha %g; log(alpha) %g\n", m_alpha, m_logalpha);
  }
  double relsum = (sumNorm/std::max(designNorm,stateNorm));
  m_log->message(2,
             "design norm %g stateNorm %g sum %g; relative difference %g\n",
             designNorm, stateNorm, sumNorm, relsum);

  // If we have an adaptive tikhonov parameter, check if we have met
  // this constraint first.
  if (m_tikhonov_adaptive) {
    double disc_ratio = fabs((sqrt(m_val_state)/m_target_misfit) - 1.);
    if (disc_ratio > m_tikhonov_ptol) {
      return GenericTerminationReason::keep_iterating();
    }
  }
  
  if (sumNorm < m_tikhonov_atol) {
    return TerminationReason::Ptr(new GenericTerminationReason(1,"TIKHONOV_ATOL"));
  }

  if (sumNorm < m_tikhonov_rtol*std::max(designNorm,stateNorm)) {
    return TerminationReason::Ptr(new GenericTerminationReason(1,"TIKHONOV_RTOL"));
  }

  if (m_iter>m_iter_max) {
    return GenericTerminationReason::max_iter();
  } else {
    return GenericTerminationReason::keep_iterating();
  }
}

TerminationReason::Ptr IP_SSATaucTikhonovGNSolver::evaluate_objective_and_gradient() {

  TerminationReason::Ptr reason = m_ssaforward.linearize_at(*m_d);
  if (reason->failed()) {
    return reason;
  }

  m_d_diff.copy_from(*m_d);
  m_d_diff.add(-1,m_d0);

  m_u_diff.copy_from(*m_ssaforward.solution());
  m_u_diff.add(-1,m_u_obs);

  m_designFunctional.gradientAt(m_d_diff,m_grad_design);

  // The following computes the reduced gradient.
  StateVec &adjointRHS = m_tmp_S1Global;
  m_stateFunctional.gradientAt(m_u_diff,adjointRHS);  
  m_ssaforward.apply_linearization_transpose(adjointRHS,m_grad_state);

  m_gradient.copy_from(m_grad_design);
  m_gradient.scale(m_alpha);    
  m_gradient.add(1,m_grad_state);

  double valDesign, valState;
  m_designFunctional.valueAt(m_d_diff,&valDesign);
  m_stateFunctional.valueAt(m_u_diff,&valState);

  m_val_design = valDesign;
  m_val_state = valState;
  
  m_value = valDesign * m_alpha + valState;

  return reason;
}

TerminationReason::Ptr IP_SSATaucTikhonovGNSolver::linesearch() {
  PetscErrorCode ierr;

  TerminationReason::Ptr step_reason;

  double old_value = m_val_design * m_alpha + m_val_state;

  double descent_derivative;

  m_tmp_D1Global.copy_from(m_h);

  ierr = VecDot(m_gradient.vec(), m_tmp_D1Global.vec(), &descent_derivative);
  PISM_CHK(ierr, "VecDot");

  if (descent_derivative >=0) {
    printf("descent derivative: %g\n",descent_derivative);
    return TerminationReason::Ptr(new GenericTerminationReason(-1, "Not descent direction"));
  }

  double alpha = 1;
  m_tmp_D1Local.copy_from(*m_d);
  while(true) {
    m_d->add(alpha,m_h);  // Replace with line search.
    step_reason = this->evaluate_objective_and_gradient();
    if (step_reason->succeeded()) {
      if (m_value <= old_value + 1e-3*alpha*descent_derivative) {
        break;
      }
    }
    else {
      printf("forward solve failed in linsearch.  Shrinking.\n");
    }
    alpha *=.5;
    if (alpha<1e-20) {
      printf("alpha= %g; derivative = %g\n",alpha,descent_derivative);
      return TerminationReason::Ptr(new GenericTerminationReason(-1, "Too many step shrinks."));
    }
    m_d->copy_from(m_tmp_D1Local);
  }
  
  return GenericTerminationReason::success();
}

TerminationReason::Ptr IP_SSATaucTikhonovGNSolver::solve() {

  if (m_target_misfit == 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Call set target misfit prior to calling"
                                  " IP_SSATaucTikhonovGNSolver::solve.");
  }

  m_iter = 0;
  m_d->copy_from(m_d0);

  double dlogalpha = 0;

  TerminationReason::Ptr step_reason, reason;

  step_reason = this->evaluate_objective_and_gradient();
  if (step_reason->failed()) {
    reason.reset(new GenericTerminationReason(-1,"Forward solve"));
    reason->set_root_cause(step_reason);
    return reason;
  }

  while(true) {

    reason = this->check_convergence();
    if (reason->done()) {
      return reason;
    }

    if (m_tikhonov_adaptive) {
      m_logalpha += dlogalpha;
      m_alpha = exp(m_logalpha);
    }

    step_reason = this->solve_linearized();
    if (step_reason->failed()) {
      reason.reset(new GenericTerminationReason(-1,"Gauss Newton solve"));
      reason->set_root_cause(step_reason);
      return reason;
    }

    step_reason = this->linesearch();
    if (step_reason->failed()) {
      TerminationReason::Ptr cause = reason;
      reason.reset(new GenericTerminationReason(-1,"Linesearch"));
      reason->set_root_cause(step_reason);
      return reason;
    }

    if (m_tikhonov_adaptive) {
      step_reason = this->compute_dlogalpha(&dlogalpha);
      if (step_reason->failed()) {
        TerminationReason::Ptr cause = reason;
        reason.reset(new GenericTerminationReason(-1,"Tikhonov penalty update"));
        reason->set_root_cause(step_reason);
        return reason;
      }
    }

    m_iter++;
  }

  return reason;
}

TerminationReason::Ptr IP_SSATaucTikhonovGNSolver::compute_dlogalpha(double *dlogalpha) {

  PetscErrorCode ierr;

  // Compute the right-hand side for computing dh/dalpha.
  m_d_diff_lin.copy_from(m_d_diff);
  m_d_diff_lin.add(1,m_h);  
  m_designFunctional.interior_product(m_d_diff_lin,m_dalpha_rhs);
  m_dalpha_rhs.scale(-1);

  // Solve linear equation for dh/dalpha. 
  ierr = KSPSetOperators(m_ksp,m_mat_GN,m_mat_GN);
  PISM_CHK(ierr, "KSPSetOperators");

  ierr = KSPSolve(m_ksp,m_dalpha_rhs.vec(),m_dh_dalphaGlobal.vec());
  PISM_CHK(ierr, "KSPSolve");

  m_dh_dalpha.copy_from(m_dh_dalphaGlobal);

  KSPConvergedReason ksp_reason;
  ierr = KSPGetConvergedReason(m_ksp,&ksp_reason);
  PISM_CHK(ierr, "KSPGetConvergedReason");

  if (ksp_reason<0) {
    return TerminationReason::Ptr(new KSPTerminationReason(ksp_reason));
  }

  // S1Local contains T(h) + F(x) - u_obs, i.e. the linearized misfit field.
  m_ssaforward.apply_linearization(m_h,m_tmp_S1Local);
  m_tmp_S1Local.update_ghosts();
  m_tmp_S1Local.add(1,m_u_diff);

  // Compute linearized discrepancy.
  double disc_sq;
  m_stateFunctional.dot(m_tmp_S1Local,m_tmp_S1Local,&disc_sq);

  // There are a number of equivalent ways to compute the derivative of the 
  // linearized discrepancy with respect to alpha, some of which are cheaper
  // than others to compute.  This equivalency relies, however, on having an 
  // exact solution in the Gauss-Newton step.  Since we only solve this with 
  // a soft tolerance, we lose equivalency.  We attempt a cheap computation,
  // and then do a sanity check (namely that the derivative is positive).
  // If this fails, we compute by a harder way that inherently yields a 
  // positive number.

  double ddisc_sq_dalpha;
  m_designFunctional.dot(m_dh_dalpha,m_d_diff_lin,&ddisc_sq_dalpha);
  ddisc_sq_dalpha *= -2*m_alpha;

  if (ddisc_sq_dalpha <= 0) {
    // Try harder.
    
    m_log->message(3,
               "Adaptive Tikhonov sanity check failed (dh/dalpha= %g <= 0)."
               " Tighten inv_gn_ksp_rtol?\n",
               ddisc_sq_dalpha);
    
    // S2Local contains T(dh/dalpha)
    m_ssaforward.apply_linearization(m_dh_dalpha,m_tmp_S2Local);
    m_tmp_S2Local.update_ghosts();

    double ddisc_sq_dalpha_a;
    m_stateFunctional.dot(m_tmp_S2Local,m_tmp_S2Local,&ddisc_sq_dalpha_a);
    double ddisc_sq_dalpha_b;
    m_designFunctional.dot(m_dh_dalpha,m_dh_dalpha,&ddisc_sq_dalpha_b);
    ddisc_sq_dalpha = 2*m_alpha*(ddisc_sq_dalpha_a+m_alpha*ddisc_sq_dalpha_b);

    m_log->message(3,
               "Adaptive Tikhonov sanity check recovery attempt: dh/dalpha= %g. \n",
               ddisc_sq_dalpha);

    // This is yet another alternative formula.
    // m_stateFunctional.dot(m_tmp_S1Local,m_tmp_S2Local,&ddisc_sq_dalpha);
    // ddisc_sq_dalpha *= 2;
  }

  // Newton's method formula.
  *dlogalpha = (m_target_misfit*m_target_misfit-disc_sq)/(ddisc_sq_dalpha*m_alpha);

  // It's easy to take steps that are too big when we are far from the solution.
  // So we limit the step size.
  double stepmax = 3;
  if (fabs(*dlogalpha)> stepmax) {
    double sgn = *dlogalpha > 0 ? 1 : -1;
    *dlogalpha = stepmax*sgn;
  }
  
  if (*dlogalpha<0) {
    *dlogalpha*=.5;
  }

  return GenericTerminationReason::success();
}

} // end of namespace inverse
} // end of namespace pism
