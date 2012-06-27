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

#include "InvSSA_LCLTikhonov.hh"
#include "InvSSATikhonov.hh"
#include <assert.h>

typedef IceModelVec2S  DesignVec;
typedef IceModelVec2V  StateVec;

// typedef TikhonovProblemListener<InverseProblem> Listener;
// typedef typename Listener::Ptr ListenerPtr;

InvSSA_LCLTikhonov::InvSSA_LCLTikhonov( InvSSATikhonov &invProb,
InvSSA_LCLTikhonov::DesignVec &d0, InvSSA_LCLTikhonov::StateVec &u_obs, PetscReal eta):
m_d0(d0), m_u_obs(u_obs), m_eta(eta), m_invProblem(invProb) {
  PetscErrorCode ierr;
  ierr = this->construct();
  assert(ierr==0);
}

PetscErrorCode InvSSA_LCLTikhonov::construct() {
  PetscErrorCode ierr;

  IceGrid &grid = *m_d0.get_grid();

  PetscReal stressScale = grid.config.get("tauc_param_tauc_scale");
  m_constraintsScale = grid.Lx*grid.Ly*4*stressScale;

  m_velocityScale = grid.config.get("inv_ssa_velocity_scale")/secpera;


  PetscInt design_stencil_width = m_d0.get_stencil_width();
  PetscInt state_stencil_width = m_u_obs.get_stencil_width();
  ierr = m_d.create(grid, "design variable", kHasGhosts, design_stencil_width); CHKERRQ(ierr);
  ierr = m_d_Jdesign.create(grid, "Jdesign design variable", kHasGhosts, design_stencil_width); CHKERRQ(ierr);
  ierr = m_dGlobal.create(grid, "design variable (global)", kNoGhosts, design_stencil_width); CHKERRQ(ierr);
  ierr = m_dGlobal.copy_from(m_d0); CHKERRQ(ierr);

  ierr = m_uGlobal.create(grid, "state variable (global)", kNoGhosts, state_stencil_width); CHKERRQ(ierr);
  ierr = m_u.create(grid, "state variable", kHasGhosts, state_stencil_width); CHKERRQ(ierr);
  ierr = m_du.create(grid, "du", kHasGhosts, state_stencil_width); CHKERRQ(ierr);
  ierr = m_u_Jdesign.create(grid, "Jdesign state variable", kHasGhosts, state_stencil_width); CHKERRQ(ierr);
  
  ierr = m_u_diff.create( grid, "state residual", kHasGhosts, state_stencil_width); CHKERRQ(ierr);
  ierr = m_d_diff.create( grid, "design residual", kHasGhosts, design_stencil_width); CHKERRQ(ierr);
  ierr = m_dzeta.create(grid,"dzeta",kHasGhosts,design_stencil_width); CHKERRQ(ierr);

  ierr = m_grad_penalty.create( grid, "penalty gradient", kNoGhosts, state_stencil_width); CHKERRQ(ierr);
  ierr = m_grad_objective.create( grid, "objective gradient", kNoGhosts, design_stencil_width); CHKERRQ(ierr);

  ierr = m_invProblem.get_da(&m_da); CHKERRQ(ierr);
  
  ierr = m_constraints.create(grid,"PDE constraints",kNoGhosts,design_stencil_width); CHKERRQ(ierr);
  
  ierr = DMGetMatrix(m_da, "baij", &m_Jstate); CHKERRQ(ierr);

  PetscInt nLocalNodes  = grid.xm*grid.ym;
  PetscInt nGlobalNodes = grid.Mx*grid.My;
  ierr = MatCreateShell(grid.com,2*nLocalNodes,nLocalNodes,2*nGlobalNodes,nGlobalNodes,this,&m_Jdesign); CHKERRQ(ierr);
  ierr = MatShellSetOperation(m_Jdesign,MATOP_MULT,(void(*)(void))InvSSA_LCLTikhonov_applyJacobianDesign); CHKERRQ(ierr);
  ierr = MatShellSetOperation(m_Jdesign,MATOP_MULT_TRANSPOSE,(void(*)(void))InvSSA_LCLTikhonov_applyJacobianDesignTranspose); CHKERRQ(ierr);

  m_x.reset(new TwoBlockVec(m_dGlobal.get_vec(),m_uGlobal.get_vec()));
  return 0;
}

InvSSA_LCLTikhonov::~InvSSA_LCLTikhonov() 
{
  PetscErrorCode ierr;
  ierr = this->destruct();
  assert(ierr==0);
}

PetscErrorCode InvSSA_LCLTikhonov::destruct() {
  PetscErrorCode ierr;
  ierr = MatDestroy(&m_Jstate); CHKERRQ(ierr);
  ierr = MatDestroy(&m_Jdesign); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSSA_LCLTikhonov::setInitialGuess( DesignVec &d0) {
  PetscErrorCode ierr;
  ierr = m_dGlobal.copy_from(d0); CHKERRQ(ierr);
  return 0;
}

InvSSA_LCLTikhonov::StateVec &InvSSA_LCLTikhonov::stateSolution() {
  // PetscErrorCode ierr;
  
  // FIXME!
  m_x->scatterToB(m_uGlobal.get_vec()); //CHKERRQ(ierr);
  m_uGlobal.scale(m_velocityScale); //CHKERRQ(ierr);
  
  return m_uGlobal;
}

InvSSA_LCLTikhonov::DesignVec &InvSSA_LCLTikhonov::designSolution() {
  m_x->scatterToA(m_d.get_vec());
  return m_d;
}

PetscErrorCode InvSSA_LCLTikhonov::connect(TaoSolver tao) {
  PetscErrorCode ierr;
  ierr = TaoSetStateDesignIS(tao, m_x->blockBIndexSet() /*state*/ , m_x->blockAIndexSet() /*design*/); CHKERRQ(ierr);
  ierr = TaoObjGradCallback<InvSSA_LCLTikhonov>::connect(tao,*this); CHKERRQ(ierr);
  ierr = TaoLCLCallbacks<InvSSA_LCLTikhonov>::connect(tao,*this,m_constraints.get_vec(),m_Jstate,m_Jdesign); CHKERRQ(ierr);
  ierr = TaoMonitorCallback<InvSSA_LCLTikhonov>::connect(tao,*this); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode InvSSA_LCLTikhonov::monitorTao(TaoSolver tao) {
  PetscErrorCode ierr;
  
  PetscInt its;
  ierr =  TaoGetSolutionStatus(tao, &its, NULL, NULL, NULL, NULL, NULL ); CHKERRQ(ierr);
  
  int nListeners = m_listeners.size();
  for(int k=0; k<nListeners; k++) {
   ierr = m_listeners[k]->iteration(*this,m_eta,
                 its,m_valObjective,m_valPenalty,
                 m_d, m_d_diff, m_grad_objective,
                 m_invProblem.solution(), m_u_diff, m_grad_penalty,
                 m_constraints); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode InvSSA_LCLTikhonov::evaluateObjectiveAndGradient(TaoSolver /*tao*/, Vec x, PetscReal *value, Vec gradient) {
  PetscErrorCode ierr;

  ierr = m_x->scatter(x,m_dGlobal.get_vec(),m_uGlobal.get_vec()); CHKERRQ(ierr);
  ierr = m_uGlobal.scale(m_velocityScale); CHKERRQ(ierr);
  
  // Variable 'm_dGlobal' has no ghosts.  We need ghosts for computation with the design variable.
  ierr = m_d.copy_from(m_dGlobal); CHKERRQ(ierr);

  ierr = m_d_diff.copy_from(m_d); CHKERRQ(ierr);
  ierr = m_d_diff.add(-1,m_d0); CHKERRQ(ierr);
  ierr = m_invProblem.evalGradObjective(m_d_diff,m_grad_objective); CHKERRQ(ierr);
  m_grad_objective.scale(1/m_eta);

  ierr = m_u_diff.copy_from(m_uGlobal); CHKERRQ(ierr);
  ierr = m_u_diff.add(-1, m_u_obs); CHKERRQ(ierr);
  ierr = m_invProblem.evalGradPenalty(m_u_diff,m_grad_penalty); CHKERRQ(ierr);
  ierr = m_grad_penalty.scale(m_velocityScale); CHKERRQ(ierr);

  m_x->gather(m_grad_objective.get_vec(),m_grad_penalty.get_vec(),gradient);

  // PetscReal go,gp,x_norm, g_norm;
  // ierr = m_grad_objective.norm(NORM_2,go); CHKERRQ(ierr);
  // ierr = m_grad_penalty.norm(NORM_2,gp); CHKERRQ(ierr);

  // ierr = VecNorm(*m_x,NORM_2,&x_norm); CHKERRQ(ierr);
  // ierr = VecNorm(gradient,NORM_2,&g_norm); CHKERRQ(ierr);
  // // ierr = VecView(*m_x,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  // printf("grad_objective %g grad_penalty, %g x residual: %g same? %g\n",go,gp,x_norm,g_norm);

  PetscReal valObjective, valPenalty;
  ierr = m_invProblem.evalObjective(m_d_diff,&valObjective); CHKERRQ(ierr);
  ierr = m_invProblem.evalPenalty(m_u_diff,&valPenalty); CHKERRQ(ierr);

  m_valObjective = valObjective;
  m_valPenalty = valPenalty;

  *value = valObjective / m_eta + valPenalty;

  return 0;
}

PetscErrorCode InvSSA_LCLTikhonov::formInitialGuess(Vec *x) {
  PetscErrorCode ierr;
  bool success;
  ierr = m_d.copy_from(m_dGlobal); CHKERRQ(ierr);
  ierr = m_invProblem.linearizeAt(m_d,success); CHKERRQ(ierr);
  ierr = m_uGlobal.copy_from(m_invProblem.solution()); CHKERRQ(ierr);
  ierr = m_uGlobal.scale(1./m_velocityScale); CHKERRQ(ierr);  
  m_uGlobal.set(0);

  ierr = m_x->gather(m_dGlobal.get_vec(),m_uGlobal.get_vec()); CHKERRQ(ierr);

  // This is probably irrelevant.
  ierr = m_uGlobal.scale(m_velocityScale); CHKERRQ(ierr);  

  *x =  *m_x;
  return 0;
}

PetscErrorCode InvSSA_LCLTikhonov::evaluateConstraints(TaoSolver, Vec x, Vec r) {
  PetscErrorCode ierr;

  ierr = m_x->scatter(x,m_dGlobal.get_vec(),m_uGlobal.get_vec()); CHKERRQ(ierr);
  ierr = m_uGlobal.scale(m_velocityScale); CHKERRQ(ierr);

  ierr = m_d.copy_from(m_dGlobal); CHKERRQ(ierr);
  ierr = m_u.copy_from(m_uGlobal); CHKERRQ(ierr);

  ierr = m_invProblem.set_zeta(m_d); CHKERRQ(ierr);

  ierr = m_invProblem.assembleFunction(m_u, r); CHKERRQ(ierr);

  ierr = VecScale(r,1./m_constraintsScale);

  return 0;
}

PetscErrorCode InvSSA_LCLTikhonov::evaluateConstraintsJacobianState(TaoSolver, Vec x, Mat *Jstate, Mat * /*Jpc*/, Mat * /*Jinv*/, MatStructure *s) {
  PetscErrorCode ierr;

  ierr = m_x->scatter(x,m_dGlobal.get_vec(),m_uGlobal.get_vec()); CHKERRQ(ierr);
  ierr = m_uGlobal.scale(m_velocityScale); CHKERRQ(ierr);
  
  ierr = m_d.copy_from(m_dGlobal); CHKERRQ(ierr);
  ierr = m_u.copy_from(m_uGlobal); CHKERRQ(ierr);

  ierr = m_invProblem.set_zeta(m_d); CHKERRQ(ierr);
  ierr = m_invProblem.assembleJacobian(m_u,*Jstate); CHKERRQ(ierr);
  *s = SAME_NONZERO_PATTERN;

  ierr = MatScale(*Jstate,m_velocityScale/m_constraintsScale); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode  InvSSA_LCLTikhonov::evaluateConstraintsJacobianDesign(TaoSolver, Vec x, Mat* /*Jdesign*/) {
  PetscErrorCode ierr;
  // I'm not sure if the following are necessary (i.e. will the copies that happen
  // in evaluateObjectiveAndGradient be sufficient) but we'll do them here
  // just in case.
  ierr = m_x->scatter(x,m_dGlobal.get_vec(),m_uGlobal.get_vec()); CHKERRQ(ierr);
  ierr = m_uGlobal.scale(m_velocityScale); CHKERRQ(ierr);
  ierr = m_d_Jdesign.copy_from(m_dGlobal); CHKERRQ(ierr);
  ierr = m_u_Jdesign.copy_from(m_uGlobal); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSSA_LCLTikhonov::applyConstraintsJacobianDesign(Vec x, Vec y) {
  PetscErrorCode ierr;
  ierr = m_dzeta.copy_from(x); CHKERRQ(ierr);

  PetscReal **zeta_a;
  ierr = m_d_Jdesign.get_array(zeta_a); CHKERRQ(ierr);

  PISMVector2 **u_a;
  ierr = m_u_Jdesign.get_array(u_a); CHKERRQ(ierr);
  
  PetscReal **dzeta_a;
  ierr = m_dzeta.get_array(dzeta_a); CHKERRQ(ierr);
  
  PISMVector2 **y_a;
  ierr = DMDAVecGetArray(m_da,y,&y_a); CHKERRQ(ierr);

  ierr = m_invProblem.assemble_T_rhs(zeta_a, u_a, dzeta_a, y_a); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(m_da,y,&y_a); CHKERRQ(ierr);
  
  ierr = VecScale(y,-1./m_constraintsScale); CHKERRQ(ierr);
  
  ierr = m_d_Jdesign.end_access(); CHKERRQ(ierr);
  ierr = m_u_Jdesign.end_access(); CHKERRQ(ierr);
  ierr = m_dzeta.end_access(); CHKERRQ(ierr);

  // PetscReal normX, normY;
  // ierr = VecNorm(x,NORM_2,&normX); CHKERRQ(ierr);
  // ierr = VecNorm(y,NORM_2,&normY); CHKERRQ(ierr);  
  // printf("applyConstraintsJacobianDesign %g %g\n",normX,normY);
  
  
  return 0;
}

PetscErrorCode InvSSA_LCLTikhonov::applyConstraintsJacobianDesignTranspose(Vec x, Vec y) {
  PetscErrorCode ierr;
  
  ierr = m_du.copy_from(x); CHKERRQ(ierr);

  PetscReal **zeta_a;
  ierr = m_d_Jdesign.get_array(zeta_a); CHKERRQ(ierr);

  PISMVector2 **u_a;
  ierr = m_u_Jdesign.get_array(u_a); CHKERRQ(ierr);
  
  PISMVector2 **du_a;
  ierr = m_du.get_array(du_a); CHKERRQ(ierr);
  
  PetscReal **y_a;

  IceGrid &grid = *m_dGlobal.get_grid();
  ierr = DMDAVecGetArray(grid.da2,y,&y_a); CHKERRQ(ierr);
  
  ierr = m_invProblem.compute_Jdesign_transpose(zeta_a, u_a, du_a, y_a); CHKERRQ(ierr);


  ierr = DMDAVecRestoreArray(grid.da2,y,&y_a); CHKERRQ(ierr);

  ierr = VecScale(y,1./m_constraintsScale); CHKERRQ(ierr);
  
  ierr = m_d_Jdesign.end_access(); CHKERRQ(ierr);
  ierr = m_u_Jdesign.end_access(); CHKERRQ(ierr);
  ierr = m_du.end_access(); CHKERRQ(ierr);

  // PetscReal normX, normY;
  // ierr = VecNorm(x,NORM_2,&normX); CHKERRQ(ierr);
  // ierr = VecNorm(y,NORM_2,&normY); CHKERRQ(ierr);  
  // printf("applyConstraintsJacobianDesignTranspose %g %g\n",normX,normY);
  
  return 0;
}

PetscErrorCode InvSSA_LCLTikhonov_applyJacobianDesign(Mat A, Vec x, Vec y) {
  PetscErrorCode ierr;
  InvSSA_LCLTikhonov *ctx;
  ierr = MatShellGetContext(A,&ctx); CHKERRQ(ierr);
  ierr = ctx->applyConstraintsJacobianDesign(x,y); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode InvSSA_LCLTikhonov_applyJacobianDesignTranspose(Mat A, Vec x, Vec y) {
  PetscErrorCode ierr;
  InvSSA_LCLTikhonov *ctx;
  ierr = MatShellGetContext(A,&ctx); CHKERRQ(ierr);
  ierr = ctx->applyConstraintsJacobianDesignTranspose(x,y); CHKERRQ(ierr);

  return 0;
}

void InvSSALCLTikhonovAddListener(InvSSA_LCLTikhonov &problem, 
                 PythonLCLTikhonovSVListener::Ptr listener ) {
  std::tr1::shared_ptr<InvSSA_LCLTikhonovPythonListenerBridge> bridge(new InvSSA_LCLTikhonovPythonListenerBridge(listener) );
  problem.addListener(bridge);
}
