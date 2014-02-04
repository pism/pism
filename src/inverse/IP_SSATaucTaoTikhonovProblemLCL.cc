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

#include "IP_SSATaucTaoTikhonovProblemLCL.hh"
#include <assert.h>

typedef IceModelVec2S  DesignVec;
typedef IceModelVec2V  StateVec;

// typedef TikhonovProblemListener<InverseProblem> Listener;
// typedef typename Listener::Ptr ListenerPtr;

IP_SSATaucTaoTikhonovProblemLCL::IP_SSATaucTaoTikhonovProblemLCL( IP_SSATaucForwardProblem &ssaforward,
IP_SSATaucTaoTikhonovProblemLCL::DesignVec &d0, IP_SSATaucTaoTikhonovProblemLCL::StateVec &u_obs, double eta,
IPFunctional<DesignVec> &designFunctional, IPFunctional<StateVec> &stateFunctional):
m_ssaforward(ssaforward), m_d0(d0), m_u_obs(u_obs), m_eta(eta),
m_designFunctional(designFunctional), m_stateFunctional(stateFunctional)
{
  PetscErrorCode ierr;
  ierr = this->construct();
  assert(ierr==0);
}

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL::construct() {
  PetscErrorCode ierr;

  IceGrid &grid = *m_d0.get_grid();

  double stressScale = grid.config.get("design_param_tauc_scale");
  m_constraintsScale = grid.Lx*grid.Ly*4*stressScale;

  m_velocityScale = grid.config.get("inv_ssa_velocity_scale", "m/year", "m/second");


  int design_stencil_width = m_d0.get_stencil_width();
  int state_stencil_width = m_u_obs.get_stencil_width();
  ierr = m_d.create(grid, "design variable", WITH_GHOSTS, design_stencil_width); CHKERRQ(ierr);
  ierr = m_d_Jdesign.create(grid, "Jdesign design variable", WITH_GHOSTS, design_stencil_width); CHKERRQ(ierr);
  ierr = m_dGlobal.create(grid, "design variable (global)", WITHOUT_GHOSTS, design_stencil_width); CHKERRQ(ierr);
  ierr = m_dGlobal.copy_from(m_d0); CHKERRQ(ierr);

  ierr = m_uGlobal.create(grid, "state variable (global)", WITHOUT_GHOSTS, state_stencil_width); CHKERRQ(ierr);
  ierr = m_u.create(grid, "state variable", WITH_GHOSTS, state_stencil_width); CHKERRQ(ierr);
  ierr = m_du.create(grid, "du", WITH_GHOSTS, state_stencil_width); CHKERRQ(ierr);
  ierr = m_u_Jdesign.create(grid, "Jdesign state variable", WITH_GHOSTS, state_stencil_width); CHKERRQ(ierr);
  
  ierr = m_u_diff.create( grid, "state residual", WITH_GHOSTS, state_stencil_width); CHKERRQ(ierr);
  ierr = m_d_diff.create( grid, "design residual", WITH_GHOSTS, design_stencil_width); CHKERRQ(ierr);
  ierr = m_dzeta.create(grid,"dzeta",WITH_GHOSTS,design_stencil_width); CHKERRQ(ierr);

  ierr = m_grad_state.create( grid, "state gradient", WITHOUT_GHOSTS, state_stencil_width); CHKERRQ(ierr);
  ierr = m_grad_design.create( grid, "design gradient", WITHOUT_GHOSTS, design_stencil_width); CHKERRQ(ierr);

  ierr = m_constraints.create(grid,"PDE constraints",WITHOUT_GHOSTS,design_stencil_width); CHKERRQ(ierr);

  DM da;
  ierr = m_ssaforward.get_da(&da); CHKERRQ(ierr);
  ierr = DMCreateMatrix(da, "baij", &m_Jstate); CHKERRQ(ierr);

  int nLocalNodes  = grid.xm*grid.ym;
  int nGlobalNodes = grid.Mx*grid.My;
  ierr = MatCreateShell(grid.com,2*nLocalNodes,nLocalNodes,2*nGlobalNodes,nGlobalNodes,this,&m_Jdesign); CHKERRQ(ierr);
  ierr = MatShellSetOperation(m_Jdesign,MATOP_MULT,(void(*)(void))IP_SSATaucTaoTikhonovProblemLCL_applyJacobianDesign); CHKERRQ(ierr);
  ierr = MatShellSetOperation(m_Jdesign,MATOP_MULT_TRANSPOSE,(void(*)(void))IP_SSATaucTaoTikhonovProblemLCL_applyJacobianDesignTranspose); CHKERRQ(ierr);

  m_x.reset(new IPTwoBlockVec(m_dGlobal.get_vec(),m_uGlobal.get_vec()));
  return 0;
}

IP_SSATaucTaoTikhonovProblemLCL::~IP_SSATaucTaoTikhonovProblemLCL() 
{
  PetscErrorCode ierr;
  ierr = this->destruct();
  assert(ierr==0);
}

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL::destruct() {
  PetscErrorCode ierr;
  ierr = MatDestroy(&m_Jstate); CHKERRQ(ierr);
  ierr = MatDestroy(&m_Jdesign); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL::setInitialGuess( DesignVec &d0) {
  PetscErrorCode ierr;
  ierr = m_dGlobal.copy_from(d0); CHKERRQ(ierr);
  return 0;
}

IP_SSATaucTaoTikhonovProblemLCL::StateVec &IP_SSATaucTaoTikhonovProblemLCL::stateSolution() {
  // PetscErrorCode ierr;
  
  // FIXME!
  m_x->scatterToB(m_uGlobal.get_vec()); //CHKERRQ(ierr);
  m_uGlobal.scale(m_velocityScale); //CHKERRQ(ierr);
  
  return m_uGlobal;
}

IP_SSATaucTaoTikhonovProblemLCL::DesignVec &IP_SSATaucTaoTikhonovProblemLCL::designSolution() {
  m_x->scatterToA(m_d.get_vec()); //CHKERRQ(ierr);
  return m_d;
}

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL::connect(TaoSolver tao) {
  PetscErrorCode ierr;
  ierr = TaoSetStateDesignIS(tao, m_x->blockBIndexSet() /*state*/ , m_x->blockAIndexSet() /*design*/); CHKERRQ(ierr);
  ierr = TaoObjGradCallback<IP_SSATaucTaoTikhonovProblemLCL,&IP_SSATaucTaoTikhonovProblemLCL::evaluateObjectiveAndGradient>::connect(tao,*this); CHKERRQ(ierr);
  ierr = TaoLCLCallbacks<IP_SSATaucTaoTikhonovProblemLCL>::connect(tao,*this,m_constraints.get_vec(),m_Jstate,m_Jdesign); CHKERRQ(ierr);
  ierr = TaoMonitorCallback<IP_SSATaucTaoTikhonovProblemLCL>::connect(tao,*this); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL::monitorTao(TaoSolver tao) {
  PetscErrorCode ierr;
  
  int its;
  ierr =  TaoGetSolutionStatus(tao, &its, NULL, NULL, NULL, NULL, NULL ); CHKERRQ(ierr);
  
  int nListeners = m_listeners.size();
  for(int k=0; k<nListeners; k++) {
   ierr = m_listeners[k]->iteration(*this,m_eta,
                 its,m_val_design,m_val_state,
                 m_d, m_d_diff, m_grad_design,
                 m_ssaforward.solution(), m_u_diff, m_grad_state,
                 m_constraints); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL::evaluateObjectiveAndGradient(TaoSolver /*tao*/, Vec x, double *value, Vec gradient) {
  PetscErrorCode ierr;

  ierr = m_x->scatter(x,m_dGlobal.get_vec(),m_uGlobal.get_vec()); CHKERRQ(ierr);
  ierr = m_uGlobal.scale(m_velocityScale); CHKERRQ(ierr);
  
  // Variable 'm_dGlobal' has no ghosts.  We need ghosts for computation with the design variable.
  ierr = m_d.copy_from(m_dGlobal); CHKERRQ(ierr);

  ierr = m_d_diff.copy_from(m_d); CHKERRQ(ierr);
  ierr = m_d_diff.add(-1,m_d0); CHKERRQ(ierr);
  ierr = m_designFunctional.gradientAt(m_d_diff,m_grad_design); CHKERRQ(ierr);
  m_grad_design.scale(1/m_eta);

  ierr = m_u_diff.copy_from(m_uGlobal); CHKERRQ(ierr);
  ierr = m_u_diff.add(-1, m_u_obs); CHKERRQ(ierr);
  ierr = m_stateFunctional.gradientAt(m_u_diff,m_grad_state); CHKERRQ(ierr);
  ierr = m_grad_state.scale(m_velocityScale); CHKERRQ(ierr);

  m_x->gather(m_grad_design.get_vec(),m_grad_state.get_vec(),gradient);

  ierr = m_designFunctional.valueAt(m_d_diff,&m_val_design); CHKERRQ(ierr);
  ierr = m_stateFunctional.valueAt(m_u_diff,&m_val_state); CHKERRQ(ierr);

  *value = m_val_design / m_eta + m_val_state;

  return 0;
}

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL::formInitialGuess(Vec *x,TerminationReason::Ptr &reason) {
  PetscErrorCode ierr;
  ierr = m_d.copy_from(m_dGlobal); CHKERRQ(ierr);
  ierr = m_ssaforward.linearize_at(m_d,reason); CHKERRQ(ierr);
  if(reason->failed()) {
    return 0;
  }
  
  ierr = m_uGlobal.copy_from(m_ssaforward.solution()); CHKERRQ(ierr);
  ierr = m_uGlobal.scale(1./m_velocityScale); CHKERRQ(ierr);  

  ierr = m_x->gather(m_dGlobal.get_vec(),m_uGlobal.get_vec()); CHKERRQ(ierr);

  // This is probably irrelevant.
  ierr = m_uGlobal.scale(m_velocityScale); CHKERRQ(ierr);  

  *x =  *m_x;
  reason = GenericTerminationReason::success();
  return 0;
}

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL::evaluateConstraints(TaoSolver, Vec x, Vec r) {
  PetscErrorCode ierr;

  ierr = m_x->scatter(x,m_dGlobal.get_vec(),m_uGlobal.get_vec()); CHKERRQ(ierr);
  ierr = m_uGlobal.scale(m_velocityScale); CHKERRQ(ierr);

  ierr = m_d.copy_from(m_dGlobal); CHKERRQ(ierr);
  ierr = m_u.copy_from(m_uGlobal); CHKERRQ(ierr);

  ierr = m_ssaforward.set_design(m_d); CHKERRQ(ierr);

  ierr = m_ssaforward.assemble_residual(m_u, r); CHKERRQ(ierr);

  ierr = VecScale(r,1./m_constraintsScale);

  return 0;
}

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL::evaluateConstraintsJacobianState(TaoSolver, Vec x, Mat *Jstate, Mat * /*Jpc*/, Mat * /*Jinv*/, MatStructure *s) {
  PetscErrorCode ierr;

  ierr = m_x->scatter(x,m_dGlobal.get_vec(),m_uGlobal.get_vec()); CHKERRQ(ierr);
  ierr = m_uGlobal.scale(m_velocityScale); CHKERRQ(ierr);
  
  ierr = m_d.copy_from(m_dGlobal); CHKERRQ(ierr);
  ierr = m_u.copy_from(m_uGlobal); CHKERRQ(ierr);

  ierr = m_ssaforward.set_design(m_d); CHKERRQ(ierr);
  ierr = m_ssaforward.assemble_jacobian_state(m_u,*Jstate); CHKERRQ(ierr);
  *s = SAME_NONZERO_PATTERN;

  ierr = MatScale(*Jstate,m_velocityScale/m_constraintsScale); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode  IP_SSATaucTaoTikhonovProblemLCL::evaluateConstraintsJacobianDesign(TaoSolver, Vec x, Mat* /*Jdesign*/) {
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

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL::applyConstraintsJacobianDesign(Vec x, Vec y) {
  PetscErrorCode ierr;
  ierr = m_dzeta.copy_from(x); CHKERRQ(ierr);
  
  ierr = m_ssaforward.set_design(m_d_Jdesign); CHKERRQ(ierr);
  
  ierr = m_ssaforward.apply_jacobian_design(m_u_Jdesign, m_dzeta, y); CHKERRQ(ierr);
  
  ierr = VecScale(y,1./m_constraintsScale); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL::applyConstraintsJacobianDesignTranspose(Vec x, Vec y) {
  PetscErrorCode ierr;

  ierr = m_du.copy_from(x); CHKERRQ(ierr);

  ierr = m_ssaforward.set_design(m_d_Jdesign); CHKERRQ(ierr);

  ierr = m_ssaforward.apply_jacobian_design_transpose(m_u_Jdesign, m_du, y); CHKERRQ(ierr);

  ierr = VecScale(y,1./m_constraintsScale); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL_applyJacobianDesign(Mat A, Vec x, Vec y) {
  PetscErrorCode ierr;
  IP_SSATaucTaoTikhonovProblemLCL *ctx;
  ierr = MatShellGetContext(A,&ctx); CHKERRQ(ierr);
  ierr = ctx->applyConstraintsJacobianDesign(x,y); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL_applyJacobianDesignTranspose(Mat A, Vec x, Vec y) {
  PetscErrorCode ierr;
  IP_SSATaucTaoTikhonovProblemLCL *ctx;
  ierr = MatShellGetContext(A,&ctx); CHKERRQ(ierr);
  ierr = ctx->applyConstraintsJacobianDesignTranspose(x,y); CHKERRQ(ierr);

  return 0;
}
