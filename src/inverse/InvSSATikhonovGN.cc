// Copyright (C) 2012, 2013  David Maxwell
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

#include "InvSSATikhonovGN.hh"
#include <assert.h>
#include "TerminationReason.hh"
#include "pism_options.hh"

InvSSATikhonovGN::InvSSATikhonovGN( InvSSAForwardProblem &ssaforward,
DesignVec &d0, StateVec &u_obs, PetscReal eta,
IPFunctional<DesignVec> &designFunctional, IPFunctional<StateVec> &stateFunctional):
m_ssaforward(ssaforward), m_d0(d0), m_u_obs(u_obs), m_eta(eta),
m_designFunctional(designFunctional), m_stateFunctional(stateFunctional)
{
  PetscErrorCode ierr;
  ierr = this->construct();
  assert(ierr==0);
}

InvSSATikhonovGN::~InvSSATikhonovGN() {
  PetscErrorCode ierr;
  ierr = this->destruct(); CHKERRCONTINUE(ierr);
  assert(ierr==0);
}


PetscErrorCode InvSSATikhonovGN::construct() {
  PetscErrorCode ierr;
  IceGrid &grid = *m_d0.get_grid();
  m_comm = grid.com;

  PetscInt design_stencil_width = m_d0.get_stencil_width();
  PetscInt state_stencil_width = m_u_obs.get_stencil_width();

  ierr = m_x.create(grid,"x",kHasGhosts,design_stencil_width); CHKERRQ(ierr);

  ierr = m_tmp_D1Global.create(grid,"work vector",kNoGhosts,0); CHKERRQ(ierr);
  ierr = m_tmp_D2Global.create(grid,"work vector",kNoGhosts,0); CHKERRQ(ierr);
  ierr = m_tmp_S1Global.create(grid,"work vector",kNoGhosts,0); CHKERRQ(ierr);
  ierr = m_tmp_S2Global.create(grid,"work vector",kNoGhosts,0); CHKERRQ(ierr);

  ierr = m_tmp_D1Local.create(grid,"work vector",kHasGhosts,design_stencil_width); CHKERRQ(ierr);
  ierr = m_tmp_D2Local.create(grid,"work vector",kHasGhosts,design_stencil_width); CHKERRQ(ierr);
  ierr = m_tmp_S1Local.create(grid,"work vector",kHasGhosts,state_stencil_width); CHKERRQ(ierr);
  ierr = m_tmp_S2Local.create(grid,"work vector",kHasGhosts,state_stencil_width); CHKERRQ(ierr);

  ierr = m_GN_rhs.create(grid,"GN_rhs",kNoGhosts,0); CHKERRQ(ierr);

  ierr = m_dGlobal.create(grid,"d (sans ghosts)",kNoGhosts,0); CHKERRQ(ierr);
  ierr = m_d.create(grid,"d",kHasGhosts,design_stencil_width); CHKERRQ(ierr);
  ierr = m_d_diff.create(grid,"d_diff",kHasGhosts,design_stencil_width); CHKERRQ(ierr);
  ierr = m_d_diff_lin.create(grid,"d_diff linearized",kHasGhosts,design_stencil_width); CHKERRQ(ierr);
  ierr = m_h.create(grid,"h",kHasGhosts,design_stencil_width); CHKERRQ(ierr);
  ierr = m_hGlobal.create(grid,"h (sans ghosts)",kNoGhosts); CHKERRQ(ierr);
  
  ierr = m_dalpha_rhs.create(grid,"dalpha rhs",kNoGhosts); CHKERRQ(ierr);
  ierr = m_dh_dalpha.create(grid,"dh_dalpha",kHasGhosts,design_stencil_width); CHKERRQ(ierr);
  ierr = m_dh_dalphaGlobal.create(grid,"dh_dalpha",kNoGhosts); CHKERRQ(ierr);
  ierr = m_u_diff.create(grid,"du",kHasGhosts,state_stencil_width); CHKERRQ(ierr);

  ierr = m_grad_design.create(grid,"grad design",kNoGhosts); CHKERRQ(ierr);
  ierr = m_grad_state.create(grid,"grad design",kNoGhosts); CHKERRQ(ierr);
  ierr = m_gradient.create(grid,"grad design",kNoGhosts); CHKERRQ(ierr);

  ierr = KSPCreate(grid.com, &m_ksp); CHKERRQ(ierr);
  ierr = KSPSetOptionsPrefix(m_ksp,"inv_gn_"); CHKERRQ(ierr);
  PetscReal ksp_rtol = 1e-5; // Soft tolerance
  ierr = KSPSetTolerances(m_ksp,ksp_rtol,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
  ierr = KSPSetType(m_ksp,KSPCG); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(m_ksp,&pc); CHKERRQ(ierr);
  ierr = PCSetType(pc,PCNONE); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(m_ksp); CHKERRQ(ierr);  

  PetscInt nLocalNodes  = grid.xm*grid.ym;
  PetscInt nGlobalNodes = grid.Mx*grid.My;
  ierr = MatCreateShell(grid.com,nLocalNodes,nLocalNodes,nGlobalNodes,nGlobalNodes,this,&m_mat_GN); CHKERRQ(ierr);

  typedef MatrixMultiplyCallback<InvSSATikhonovGN,&InvSSATikhonovGN::apply_GN> multCallback;
  ierr = multCallback::connect(m_mat_GN);

  m_alpha = 1./m_eta;
  m_logalpha = log(m_alpha);
  m_vel_scale = grid.config.get("inv_ssa_velocity_scale");
  m_rms_error = grid.config.get("inv_ssa_target_rms_misfit")/m_vel_scale;

  ierr = PISMOptionsIsSet("-tikhonov_adaptive", m_tikhonov_adaptive); CHKERRQ(ierr);
  
  m_iter_max = 1000; bool flag;
  ierr = PISMOptionsInt("-inv_gn_iter_max", "",m_iter_max,flag); CHKERRQ(ierr);  

  m_tikhonov_atol = grid.config.get("tikhonov_atol");
  m_tikhonov_rtol = grid.config.get("tikhonov_rtol");
  m_tikhonov_ptol = grid.config.get("tikhonov_ptol");

  return 0;
}

PetscErrorCode InvSSATikhonovGN::destruct() {
  PetscErrorCode ierr;
  ierr = KSPDestroy(&m_ksp); CHKERRQ(ierr);
  ierr = MatDestroy(&m_mat_GN);
  return 0;
}

PetscErrorCode InvSSATikhonovGN::init(TerminationReason::Ptr &reason) {
  PetscErrorCode ierr;
  ierr = m_ssaforward.linearize_at(m_d0,reason); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSATikhonovGN::apply_GN(IceModelVec2S &x,IceModelVec2S &y) {
  PetscErrorCode ierr;
  ierr = this->apply_GN(x.get_vec(),y.get_vec()); CHKERRQ(ierr);
  return 0; 
}

PetscErrorCode InvSSATikhonovGN::apply_GN(Vec x, Vec y) {
  PetscErrorCode ierr;

  StateVec &tmp_gS    = m_tmp_S1Global;
  StateVec &Tx        = m_tmp_S1Local;
  DesignVec &tmp_gD   = m_tmp_D1Global;
  DesignVec  &GNx      = m_tmp_D2Global;
  
  // FIXME: Needless copies for now.
  ierr = m_x.copy_from(x); CHKERRQ(ierr);

  ierr = m_ssaforward.apply_linearization(m_x,Tx); CHKERRQ(ierr);
  ierr = Tx.update_ghosts(); CHKERRQ(ierr);
  
  ierr = m_stateFunctional.interior_product(Tx,tmp_gS); CHKERRQ(ierr);
  
  ierr = m_ssaforward.apply_linearization_transpose(tmp_gS,GNx); CHKERRQ(ierr);

  ierr = m_designFunctional.interior_product(m_x,tmp_gD); CHKERRQ(ierr);
  ierr = GNx.add(m_alpha,tmp_gD); CHKERRQ(ierr);

  ierr = GNx.copy_to(y); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSSATikhonovGN::assemble_GN_rhs(DesignVec &rhs) {
  PetscErrorCode ierr;

  ierr = rhs.set(0); CHKERRQ(ierr);
  
  ierr = m_stateFunctional.interior_product(m_u_diff,m_tmp_S1Global); CHKERRQ(ierr);
  ierr = m_ssaforward.apply_linearization_transpose(m_tmp_S1Global,rhs); CHKERRQ(ierr);

  ierr = m_designFunctional.interior_product(m_d_diff,m_tmp_D1Global); CHKERRQ(ierr);
  ierr = rhs.add(m_alpha,m_tmp_D1Global);
  
  ierr = rhs.scale(-1); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSSATikhonovGN::solve_linearized(TerminationReason::Ptr &reason) {
  PetscErrorCode ierr;

  ierr = this->assemble_GN_rhs(m_GN_rhs); CHKERRQ(ierr);

  ierr = KSPSetOperators(m_ksp,m_mat_GN,m_mat_GN,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_ksp,m_GN_rhs.get_vec(),m_hGlobal.get_vec()); CHKERRQ(ierr);

  KSPConvergedReason ksp_reason;
  ierr = KSPGetConvergedReason(m_ksp,&ksp_reason); CHKERRQ(ierr);
  
  ierr = m_h.copy_from(m_hGlobal); CHKERRQ(ierr);

  reason.reset( new KSPTerminationReason(ksp_reason) );

  return 0;
}

PetscErrorCode InvSSATikhonovGN::evaluateGNFunctional(DesignVec h, PetscReal *value) {
  PetscErrorCode ierr;
  
  ierr = m_ssaforward.apply_linearization(h,m_tmp_S1Local); CHKERRQ(ierr);
  ierr = m_tmp_S1Local.update_ghosts(); CHKERRQ(ierr);
  ierr = m_tmp_S1Local.add(1,m_u_diff);
  
  PetscReal sValue;
  ierr =  m_stateFunctional.valueAt(m_tmp_S1Local,&sValue); CHKERRQ(ierr);
  
  
  ierr = m_tmp_D1Local.copy_from(m_d_diff); CHKERRQ(ierr);
  ierr = m_tmp_D1Local.add(1,h); CHKERRQ(ierr);
  
  PetscReal dValue;
  ierr =  m_designFunctional.valueAt(m_tmp_D1Local,&dValue); CHKERRQ(ierr);
  
  *value = m_alpha*dValue + sValue;

  return 0;
}


PetscErrorCode InvSSATikhonovGN::check_convergence(TerminationReason::Ptr &reason) {
  PetscErrorCode ierr;

  PetscReal designNorm, stateNorm, sumNorm;
  PetscReal dWeight, sWeight;
  dWeight = m_alpha;
  sWeight = 1;

  ierr = m_grad_design.norm(NORM_2,designNorm); CHKERRQ(ierr);
  ierr = m_grad_state.norm(NORM_2,stateNorm); CHKERRQ(ierr);
  designNorm *= dWeight;
  stateNorm  *= sWeight;

  ierr = m_gradient.norm(NORM_2,sumNorm); CHKERRQ(ierr);

  ierr = verbPrintf(2,PETSC_COMM_WORLD,"----------------------------------------------------------\n",designNorm,stateNorm,sumNorm); CHKERRQ(ierr);
  ierr = verbPrintf(2,PETSC_COMM_WORLD,"InvSSATikhonovGN Iteration %d: misfit %g; functional %g \n",m_iter,sqrt(m_val_state)*m_vel_scale,m_value*m_vel_scale*m_vel_scale); CHKERRQ(ierr);
  if(m_tikhonov_adaptive) {
    ierr = verbPrintf(2,PETSC_COMM_WORLD,"alpha %g; log(alpha) %g\n",m_alpha,m_logalpha); CHKERRQ(ierr);
  }
  PetscReal relsum = (sumNorm/PetscMax(designNorm,stateNorm));
  ierr = verbPrintf(2,PETSC_COMM_WORLD,"design norm %g stateNorm %g sum %g; relative difference %g\n",designNorm,stateNorm,sumNorm,relsum); CHKERRQ(ierr);

  // If we have an adaptive tikhonov parameter, check if we have met
  // this constraint first.
  if(m_tikhonov_adaptive) {
    PetscReal disc_ratio = fabs( (sqrt(m_val_state)/m_rms_error) - 1.);
    if(disc_ratio > m_tikhonov_ptol) {
      reason = GenericTerminationReason::keep_iterating();
      return 0;
    }
  }
  
  if(sumNorm < m_tikhonov_atol) {
    reason.reset(new GenericTerminationReason(1,"TIKHONOV_ATOL"));
    return 0;
  }

  if( sumNorm < m_tikhonov_rtol*PetscMax(designNorm,stateNorm) ) {
    reason.reset(new GenericTerminationReason(1,"TIKHONOV_RTOL"));
    return 0;
  }

  if(m_iter>m_iter_max) {
    reason = GenericTerminationReason::max_iter();
  } else {
    reason = GenericTerminationReason::keep_iterating();
  }
  return 0;
}

PetscErrorCode InvSSATikhonovGN::evaluate_objective_and_gradient(TerminationReason::Ptr &reason) {
  PetscErrorCode ierr;

  ierr = m_ssaforward.linearize_at(m_d,reason); CHKERRQ(ierr);
  if(reason->failed()) {
    return 0;
  }

  ierr = m_d_diff.copy_from(m_d); CHKERRQ(ierr);
  ierr = m_d_diff.add(-1,m_d0); CHKERRQ(ierr);

  ierr = m_u_diff.copy_from(m_ssaforward.solution()); CHKERRQ(ierr);
  ierr = m_u_diff.add(-1,m_u_obs); CHKERRQ(ierr);

  ierr = m_designFunctional.gradientAt(m_d_diff,m_grad_design); CHKERRQ(ierr);

  // The following computes the reduced gradient.
  StateVec &adjointRHS = m_tmp_S1Global;
  ierr = m_stateFunctional.gradientAt(m_u_diff,adjointRHS); CHKERRQ(ierr);  
  ierr = m_ssaforward.apply_linearization_transpose(adjointRHS,m_grad_state); CHKERRQ(ierr);

  ierr = m_gradient.copy_from(m_grad_design); CHKERRQ(ierr);
  ierr = m_gradient.scale(m_alpha); CHKERRQ(ierr);    
  ierr = m_gradient.add(1,m_grad_state); CHKERRQ(ierr);

  PetscReal valDesign, valState;
  ierr = m_designFunctional.valueAt(m_d_diff,&valDesign); CHKERRQ(ierr);
  ierr = m_stateFunctional.valueAt(m_u_diff,&valState); CHKERRQ(ierr);

  m_val_design = valDesign;
  m_val_state = valState;
  
  m_value = valDesign * m_alpha + valState;

  return 0;
}

PetscErrorCode InvSSATikhonovGN::linesearch(TerminationReason::Ptr &reason) {
  PetscErrorCode ierr;

  TerminationReason::Ptr step_reason;

  PetscReal old_value = m_val_design * m_alpha + m_val_state;

  PetscReal descent_derivative;

  ierr = m_tmp_D1Global.copy_from(m_h); CHKERRQ(ierr);
  ierr = VecDot(m_gradient.get_vec(),m_tmp_D1Global.get_vec(),&descent_derivative);

  if(descent_derivative >=0 ) {
    printf("descent derivative: %g\n",descent_derivative);
    reason.reset(new GenericTerminationReason(-1,"Not descent direction"));
    return 0;
  }

  PetscReal alpha = 1;
  ierr = m_tmp_D1Local.copy_from(m_d); CHKERRQ(ierr);
  while(true) {
    ierr = m_d.add(alpha,m_h); CHKERRQ(ierr);  // Replace with line search.
    ierr = this->evaluate_objective_and_gradient(step_reason); CHKERRQ(ierr);
    if(step_reason->succeeded()) {
      if(m_value <= old_value + 1e-3*alpha*descent_derivative) {
        break;
      }
    }
    else {
      printf("forward solve failed in linsearch.  Shrinking.\n");
    }
    alpha *=.5;
    if(alpha<1e-20) {
      printf("alpha= %g; derivative = %g\n",alpha,descent_derivative);
      reason.reset(new GenericTerminationReason(-1,"Too many step shrinks."));
      return 0;
    }
    ierr = m_d.copy_from(m_tmp_D1Local); CHKERRQ(ierr);
  }
  
  reason = GenericTerminationReason::success();
  return 0;
}

PetscErrorCode InvSSATikhonovGN::solve(TerminationReason::Ptr &reason) {
  PetscErrorCode ierr;

  m_iter = 0;
  ierr = m_d.copy_from(m_d0); CHKERRQ(ierr);

  PetscReal dlogalpha = 0;

  TerminationReason::Ptr step_reason;

  this->evaluate_objective_and_gradient(step_reason);
  if(step_reason->failed()) {
    reason.reset(new GenericTerminationReason(-1,"Forward solve"));
    reason->set_root_cause(step_reason);
    return 0;
  }

  while(true) {

    ierr = this->check_convergence(reason); CHKERRQ(ierr);
    if(reason->done()) {
      return 0;
    }

    if(m_tikhonov_adaptive) {
      m_logalpha += dlogalpha;
      m_alpha = exp(m_logalpha);
    }

    ierr = this->solve_linearized(step_reason); CHKERRQ(ierr);
    if(step_reason->failed()) {
      reason.reset(new GenericTerminationReason(-1,"Gauss Newton solve"));
      reason->set_root_cause(step_reason);
      return 0;
    }

    ierr = this->linesearch(step_reason); CHKERRQ(ierr);
    if(step_reason->failed()) {
      TerminationReason::Ptr cause = reason;
      reason.reset(new GenericTerminationReason(-1,"Linesearch"));
      reason->set_root_cause(step_reason);
      return 0;
    }

    if(m_tikhonov_adaptive) {
      ierr = this->compute_dlogalpha(&dlogalpha,step_reason); CHKERRQ(ierr);
      if(step_reason->failed()) {
        TerminationReason::Ptr cause = reason;
        reason.reset(new GenericTerminationReason(-1,"Tikhonov penalty update"));
        reason->set_root_cause(step_reason);
        return 0;
      }
    }

    m_iter++;
  }
  return 0;
}

PetscErrorCode InvSSATikhonovGN::compute_dlogalpha(PetscReal *dlogalpha, TerminationReason::Ptr &reason) {

  PetscErrorCode ierr;

  // Compute the right-hand side for computing dh/dalpha.
  ierr = m_d_diff_lin.copy_from(m_d_diff); CHKERRQ(ierr);
  ierr = m_d_diff_lin.add(1,m_h); CHKERRQ(ierr);  
  ierr = m_designFunctional.interior_product(m_d_diff_lin,m_dalpha_rhs); CHKERRQ(ierr);
  ierr = m_dalpha_rhs.scale(-1);

  // Solve linear equation for dh/dalpha. 
  ierr = KSPSetOperators(m_ksp,m_mat_GN,m_mat_GN,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_ksp,m_dalpha_rhs.get_vec(),m_dh_dalphaGlobal.get_vec()); CHKERRQ(ierr);
  ierr = m_dh_dalpha.copy_from(m_dh_dalphaGlobal); CHKERRQ(ierr);

  KSPConvergedReason ksp_reason;
  ierr = KSPGetConvergedReason(m_ksp,&ksp_reason); CHKERRQ(ierr);
  if(ksp_reason<0) {
    reason.reset( new KSPTerminationReason(ksp_reason) );
    return 0;
  }

  // S1Local contains T(h) + F(x) - u_obs, i.e. the linearized misfit field.
  ierr = m_ssaforward.apply_linearization(m_h,m_tmp_S1Local); CHKERRQ(ierr);
  ierr = m_tmp_S1Local.update_ghosts(); CHKERRQ(ierr);
  ierr = m_tmp_S1Local.add(1,m_u_diff); CHKERRQ(ierr);

  // Compute linearized discrepancy.
  PetscReal disc_sq;
  ierr = m_stateFunctional.dot(m_tmp_S1Local,m_tmp_S1Local,&disc_sq); CHKERRQ(ierr);

  // There are a number of equivalent ways to compute the derivative of the 
  // linearized discrepancy with respect to alpha, some of which are cheaper
  // than others to compute.  This equivalency relies, however, on having an 
  // exact solution in the Gauss-Newton step.  Since we only solve this with 
  // a soft tolerance, we lose equivalency.  We attempt a cheap computation,
  // and then do a sanity check (namely that the derivative is positive).
  // If this fails, we compute by a harder way that inherently yields a 
  // positive number.

  PetscReal ddisc_sq_dalpha;
  ierr = m_designFunctional.dot(m_dh_dalpha,m_d_diff_lin,&ddisc_sq_dalpha);
  ddisc_sq_dalpha *= -2*m_alpha;

  if(ddisc_sq_dalpha <= 0) {
    // Try harder.
    
    ierr = verbPrintf(3,PETSC_COMM_WORLD,"Adaptive Tikhonov sanity check failed (dh/dalpha= %g <= 0).  Tighten inv_gn_ksp_rtol?\n",ddisc_sq_dalpha); CHKERRQ(ierr);
    
    // S2Local contains T(dh/dalpha)
    ierr = m_ssaforward.apply_linearization(m_dh_dalpha,m_tmp_S2Local); CHKERRQ(ierr);
    ierr = m_tmp_S2Local.update_ghosts(); CHKERRQ(ierr);

    PetscReal ddisc_sq_dalpha_a;
    ierr = m_stateFunctional.dot(m_tmp_S2Local,m_tmp_S2Local,&ddisc_sq_dalpha_a); CHKERRQ(ierr);
    PetscReal ddisc_sq_dalpha_b;
    ierr = m_designFunctional.dot(m_dh_dalpha,m_dh_dalpha,&ddisc_sq_dalpha_b); CHKERRQ(ierr);
    ddisc_sq_dalpha = 2*m_alpha*(ddisc_sq_dalpha_a+m_alpha*ddisc_sq_dalpha_b);

    ierr = verbPrintf(3,PETSC_COMM_WORLD,"Adaptive Tikhonov sanity check recovery attempt: dh/dalpha= %g. \n",ddisc_sq_dalpha); CHKERRQ(ierr);

    // This is yet another alternative formula.
    // ierr = m_stateFunctional.dot(m_tmp_S1Local,m_tmp_S2Local,&ddisc_sq_dalpha); CHKERRQ(ierr);
    // ddisc_sq_dalpha *= 2;
  }

  // Newton's method formula.
  *dlogalpha = (m_rms_error*m_rms_error-disc_sq)/(ddisc_sq_dalpha*m_alpha);

  // It's easy to take steps that are too big when we are far from the solution.
  // So we limit the step size.
  PetscReal stepmax = 3;
  if(fabs(*dlogalpha)> stepmax) {
    PetscReal sgn = *dlogalpha > 0 ? 1 : -1;
    *dlogalpha = stepmax*sgn;
  }
  
  if(*dlogalpha<0) {
    *dlogalpha*=.5;
  }

  reason = GenericTerminationReason::success();

  return 0;
}
