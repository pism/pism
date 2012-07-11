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
  ierr = m_h.create(grid,"h",kHasGhosts,design_stencil_width); CHKERRQ(ierr);
  ierr = m_hGlobal.create(grid,"h (sans ghosts)",kNoGhosts); CHKERRQ(ierr);
  
  ierr = m_dalpha_rhs.create(grid,"dalpha rhs",kNoGhosts); CHKERRQ(ierr);
  ierr = m_dh_dalpha.create(grid,"dh_dalpha",kHasGhosts,design_stencil_width); CHKERRQ(ierr);
  ierr = m_dh_dalphaGlobal.create(grid,"dh_dalpha",kNoGhosts); CHKERRQ(ierr);
  ierr = m_u_diff.create(grid,"du",kHasGhosts,state_stencil_width); CHKERRQ(ierr);

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
  m_vel_scale = grid.config.get("inv_ssa_velocity_scale");
  m_rms_error = grid.config.get("inv_ssa_target_rms_misfit")/m_vel_scale;

  ierr = PISMOptionsIsSet("-tikhonov_adaptive", m_tikhonov_adaptive); CHKERRQ(ierr);
  
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
  ierr = Tx.beginGhostComm(); CHKERRQ(ierr);
  ierr = Tx.endGhostComm(); CHKERRQ(ierr);
  
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
  ierr = m_tmp_S1Local.beginGhostComm(); CHKERRQ(ierr);
  ierr = m_tmp_S1Local.endGhostComm(); CHKERRQ(ierr);
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
  
  PetscReal dVal, sVal;
  ierr = m_designFunctional.valueAt(m_d_diff,&dVal); CHKERRQ(ierr);
  ierr = m_stateFunctional.valueAt(m_u_diff,&sVal); CHKERRQ(ierr);
  PetscReal F = m_alpha*dVal + sVal;
  
  printf("InvSSATikhonovGN Iteration %d: misfit %g; functional %g alpha %g\n",m_iter,sqrt(sVal)*m_vel_scale,F*m_vel_scale*m_vel_scale,m_alpha);
  
  PetscInt iter_max = 10; bool flag;
  ierr = PISMOptionsInt("-inv_gn_iter_max", "",iter_max,flag); CHKERRQ(ierr);

  if(m_iter>iter_max) {
    reason = GenericTerminationReason::max_iter();
  } else {
    reason = GenericTerminationReason::keep_iterating();
  }
  return 0;
}

PetscErrorCode InvSSATikhonovGN::solve(TerminationReason::Ptr &reason) {
  PetscErrorCode ierr;

  m_iter = 0;
  ierr = m_d.copy_from(m_d0); CHKERRQ(ierr);

  TerminationReason::Ptr step_reason;
  while(true) {
    ierr = m_ssaforward.linearize_at(m_d,step_reason); CHKERRQ(ierr);
    if(step_reason->failed()) {
      reason.reset(new GenericTerminationReason(-1,"Forward solve"));
      reason->set_root_cause(step_reason);
      return 0;
    }

    ierr = m_d_diff.copy_from(m_d); CHKERRQ(ierr);
    ierr = m_d_diff.add(-1,m_d0); CHKERRQ(ierr);

    ierr = m_u_diff.copy_from(m_ssaforward.solution()); CHKERRQ(ierr);
    ierr = m_u_diff.add(-1,m_u_obs); CHKERRQ(ierr);

    ierr = this->check_convergence(reason); CHKERRQ(ierr);
    if(reason->done()) {
      return 0;
    }

    ierr = this->solve_linearized(step_reason); CHKERRQ(ierr);
    if(step_reason->failed()) {
      reason.reset(new GenericTerminationReason(-1,"Gauss Newton solve"));
      reason->set_root_cause(step_reason);
      return 0;
    }

    PetscReal dalpha = 0;
    if(m_tikhonov_adaptive) {
      ierr = this->compute_dalpha(&dalpha,step_reason); CHKERRQ(ierr);
      if(step_reason->failed()) {
        TerminationReason::Ptr cause = reason;
        reason.reset(new GenericTerminationReason(-1,"Tikhonov penalty update"));
        reason->set_root_cause(step_reason);
        return 0;
      }
    }

    ierr = m_d.add(1,m_h); CHKERRQ(ierr);  // Replace with line search.

    m_alpha += dalpha;
    m_iter++;
  }
  return 0;
}

PetscErrorCode InvSSATikhonovGN::assemble_dalpha_rhs(DesignVec &rhs) {
  PetscErrorCode ierr;
  ierr = m_tmp_D1Local.copy_from(m_d_diff); CHKERRQ(ierr);
  ierr = m_tmp_D1Local.add(1,m_h); CHKERRQ(ierr);  
  ierr = m_designFunctional.interior_product(m_tmp_D1Local,rhs); CHKERRQ(ierr);
  ierr = rhs.scale(-1); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSATikhonovGN::compute_dalpha(PetscReal *dalpha, TerminationReason::Ptr &reason) {
  PetscErrorCode ierr;

  ierr = this->assemble_dalpha_rhs(m_dalpha_rhs); CHKERRQ(ierr);
 
  ierr = KSPSetOperators(m_ksp,m_mat_GN,m_mat_GN,SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_ksp,m_dalpha_rhs.get_vec(),m_dh_dalphaGlobal.get_vec()); CHKERRQ(ierr);
  ierr = m_dh_dalpha.copy_from(m_dh_dalphaGlobal); CHKERRQ(ierr);

  KSPConvergedReason ksp_reason;
  ierr = KSPGetConvergedReason(m_ksp,&ksp_reason); CHKERRQ(ierr);
  if(ksp_reason<0) {
    reason.reset( new KSPTerminationReason(ksp_reason) );
    return 0;
  }

  ierr = m_tmp_D1Local.copy_from(m_d_diff); CHKERRQ(ierr);
  ierr = m_tmp_D1Local.add(1,m_h); CHKERRQ(ierr);

  ierr = m_ssaforward.apply_linearization(m_h,m_tmp_S1Local); CHKERRQ(ierr);
  ierr = m_tmp_S1Local.beginGhostComm(); CHKERRQ(ierr);
  ierr = m_tmp_S1Local.endGhostComm(); CHKERRQ(ierr);
  ierr = m_tmp_S1Local.add(1,m_u_diff); CHKERRQ(ierr);

  ierr = m_designFunctional.interior_product(m_tmp_D1Local,m_tmp_D1Global); CHKERRQ(ierr);
  ierr = m_tmp_D1Global.scale(m_alpha); CHKERRQ(ierr);
  
  ierr = m_stateFunctional.interior_product(m_tmp_S1Local,m_tmp_S1Global); CHKERRQ(ierr);
  ierr = m_ssaforward.apply_linearization_transpose(m_tmp_S1Global,m_tmp_D2Global); CHKERRQ(ierr);
  
  PetscReal n1,n2,dn;
  ierr = m_tmp_D1Global.norm(NORM_2,n1); CHKERRQ(ierr);
  ierr = m_tmp_D2Global.norm(NORM_2,n2); CHKERRQ(ierr);
  ierr = m_tmp_D2Global.add(-1,m_tmp_D1Global); CHKERRQ(ierr);
  ierr = m_tmp_D2Global.norm(NORM_2,dn); CHKERRQ(ierr);
  
  printf("Norms Adx %g T^tB(Tx-y) %g diff %g\n",n1,n2,dn);

  PetscReal disc_sq;
  ierr = m_stateFunctional.dot(m_tmp_S1Local,m_tmp_S1Local,&disc_sq); CHKERRQ(ierr);

  ierr = m_ssaforward.apply_linearization(m_dh_dalpha,m_tmp_S2Local); CHKERRQ(ierr);
  ierr = m_tmp_S2Local.beginGhostComm(); CHKERRQ(ierr);
  ierr = m_tmp_S2Local.endGhostComm(); CHKERRQ(ierr);

  PetscReal ddisc_sq_dalpha;
  ierr = m_stateFunctional.dot(m_tmp_S1Local,m_tmp_S2Local,&ddisc_sq_dalpha); CHKERRQ(ierr);
  ddisc_sq_dalpha *= 2;

  PetscReal ddisc_sq_dalpha4_a;
  ierr = m_stateFunctional.dot(m_tmp_S2Local,m_tmp_S2Local,&ddisc_sq_dalpha4_a); CHKERRQ(ierr);
  PetscReal ddisc_sq_dalpha4_b;
  ierr = m_designFunctional.dot(m_dh_dalpha,m_dh_dalpha,&ddisc_sq_dalpha4_b); CHKERRQ(ierr);
  PetscReal ddisc_sq_dalpha4 = 2*m_alpha*(ddisc_sq_dalpha4_a+m_alpha*ddisc_sq_dalpha4_b);
  
  PetscReal ddisc_sq_dalpha2;
  ierr = m_designFunctional.dot(m_dh_dalpha,m_tmp_D1Local,&ddisc_sq_dalpha2);
  ddisc_sq_dalpha2 *= -2*m_alpha;
  
  PetscReal ddisc_sq_dalpha3;
  ierr = this->apply_GN(m_dh_dalphaGlobal,m_tmp_D2Global);
  ierr = VecDot(m_dh_dalphaGlobal.get_vec(),m_tmp_D2Global.get_vec(),&ddisc_sq_dalpha3);
  ddisc_sq_dalpha3 *= 2*m_alpha;

  printf("four derivatives: %g %.10g %.10g %.10g\n",ddisc_sq_dalpha,ddisc_sq_dalpha2,ddisc_sq_dalpha3,ddisc_sq_dalpha4);

  *dalpha = (m_rms_error*m_rms_error-disc_sq)/ddisc_sq_dalpha2;

  printf("disc_sq %g ddisc_sq_dalpha %g alt %g\n",disc_sq,ddisc_sq_dalpha,ddisc_sq_dalpha2);

  printf("disc %.10g desired %.10g alpha %.10g dalpha %.10g\n",sqrt(disc_sq)*m_vel_scale,m_rms_error*m_vel_scale,m_alpha,*dalpha);

  reason = GenericTerminationReason::success();
  return 0;
}
