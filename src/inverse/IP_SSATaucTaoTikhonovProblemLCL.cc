// Copyright (C) 2012, 2014, 2015, 2016, 2017  David Maxwell and Constantine Khroulev
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
#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"

namespace pism {
namespace inverse {

typedef IceModelVec2S  DesignVec;
typedef IceModelVec2V  StateVec;

// typedef TikhonovProblemListener<InverseProblem> Listener;
// typedef typename Listener::Ptr ListenerPtr;

IP_SSATaucTaoTikhonovProblemLCL::IP_SSATaucTaoTikhonovProblemLCL(IP_SSATaucForwardProblem &ssaforward,
                                                                 IP_SSATaucTaoTikhonovProblemLCL::DesignVec &d0,
                                                                 IP_SSATaucTaoTikhonovProblemLCL::StateVec &u_obs,
                                                                 double eta,
                                                                 IPFunctional<DesignVec> &designFunctional,
                                                                 IPFunctional<StateVec> &stateFunctional)
: m_ssaforward(ssaforward), m_d0(d0), m_u_obs(u_obs), m_eta(eta),
  m_designFunctional(designFunctional), m_stateFunctional(stateFunctional) {

  PetscErrorCode ierr;
  IceGrid::ConstPtr grid = m_d0.grid();

  double stressScale = grid->ctx()->config()->get_double("inverse.design.param_tauc_scale");
  m_constraintsScale = grid->Lx()*grid->Ly()*4*stressScale;

  m_velocityScale = grid->ctx()->config()->get_double("inverse.ssa.velocity_scale", "m second-1");


  int design_stencil_width = m_d0.stencil_width();
  int state_stencil_width = m_u_obs.stencil_width();
  m_d.reset(new DesignVec(grid, "design variable", WITH_GHOSTS, design_stencil_width));

  m_d_Jdesign.create(grid, "Jdesign design variable", WITH_GHOSTS, design_stencil_width);
  m_dGlobal.create(grid, "design variable (global)", WITHOUT_GHOSTS, design_stencil_width);
  m_dGlobal.copy_from(m_d0);

  m_uGlobal.reset(new StateVec(grid, "state variable (global)",
                               WITHOUT_GHOSTS, state_stencil_width));

  m_u.create(grid, "state variable", WITH_GHOSTS, state_stencil_width);
  m_du.create(grid, "du", WITH_GHOSTS, state_stencil_width);
  m_u_Jdesign.create(grid, "Jdesign state variable", WITH_GHOSTS, state_stencil_width);

  m_u_diff.reset(new StateVec(grid, "state residual", WITH_GHOSTS, state_stencil_width));

  m_d_diff.reset(new DesignVec(grid, "design residual", WITH_GHOSTS, design_stencil_width));

  m_dzeta.create(grid,"dzeta",WITH_GHOSTS,design_stencil_width);

  m_grad_state.reset(new StateVec(grid, "state gradient", WITHOUT_GHOSTS, state_stencil_width));

  m_grad_design.reset(new DesignVec(grid, "design gradient", WITHOUT_GHOSTS, design_stencil_width));

  m_constraints.reset(new StateVec(grid,"PDE constraints",WITHOUT_GHOSTS,design_stencil_width));

  DM da;
  m_ssaforward.get_da(&da);

  ierr = DMSetMatType(da, MATBAIJ);
  PISM_CHK(ierr, "DMSetMatType");

  ierr = DMCreateMatrix(da, m_Jstate.rawptr());
  PISM_CHK(ierr, "DMCreateMatrix");

  int nLocalNodes  = grid->xm()*grid->ym();
  int nGlobalNodes = grid->Mx()*grid->My();
  ierr = MatCreateShell(grid->com, 2*nLocalNodes, nLocalNodes, 2*nGlobalNodes, nGlobalNodes,
                        this, m_Jdesign.rawptr());
  PISM_CHK(ierr, "MatCreateShell");

  ierr = MatShellSetOperation(m_Jdesign, MATOP_MULT,
                              (void(*)(void))jacobian_design_callback);
  PISM_CHK(ierr, "MatShellSetOperation");

  ierr = MatShellSetOperation(m_Jdesign, MATOP_MULT_TRANSPOSE,
                              (void(*)(void))jacobian_design_transpose_callback);
  PISM_CHK(ierr, "MatShellSetOperation");

  m_x.reset(new IPTwoBlockVec(m_dGlobal.vec(),m_uGlobal->vec()));
}

IP_SSATaucTaoTikhonovProblemLCL::~IP_SSATaucTaoTikhonovProblemLCL()
{
  // empty
}

void IP_SSATaucTaoTikhonovProblemLCL::setInitialGuess(DesignVec &d0) {
  m_dGlobal.copy_from(d0);
}

IP_SSATaucTaoTikhonovProblemLCL::StateVec::Ptr IP_SSATaucTaoTikhonovProblemLCL::stateSolution() {

  m_x->scatterToB(m_uGlobal->vec());
  m_uGlobal->scale(m_velocityScale);

  return m_uGlobal;
}

IP_SSATaucTaoTikhonovProblemLCL::DesignVec::Ptr IP_SSATaucTaoTikhonovProblemLCL::designSolution() {
  m_x->scatterToA(m_d->vec()); //CHKERRQ(ierr);
  return m_d;
}

void IP_SSATaucTaoTikhonovProblemLCL::connect(Tao tao) {
  PetscErrorCode ierr;
  ierr = TaoSetStateDesignIS(tao,
                             m_x->blockBIndexSet() /*state*/,
                             m_x->blockAIndexSet() /*design*/);
  PISM_CHK(ierr, "TaoSetStateDesignIS");

  taoutil::TaoObjGradCallback<IP_SSATaucTaoTikhonovProblemLCL,
                              &IP_SSATaucTaoTikhonovProblemLCL::evaluateObjectiveAndGradient>::connect(tao, *this);

  taoutil::TaoLCLCallbacks<IP_SSATaucTaoTikhonovProblemLCL>::connect(tao, *this,
                                                            m_constraints->vec(),
                                                            m_Jstate, m_Jdesign);

  taoutil::TaoMonitorCallback<IP_SSATaucTaoTikhonovProblemLCL>::connect(tao,*this);
}

void IP_SSATaucTaoTikhonovProblemLCL::monitorTao(Tao tao) {
  PetscErrorCode ierr;

  // Has to be a PetscInt because of the TaoGetSolutionStatus call.
  PetscInt its;
  ierr = TaoGetSolutionStatus(tao, &its, NULL, NULL, NULL, NULL, NULL);
  PISM_CHK(ierr, "TaoGetSolutionStatus");

  int nListeners = m_listeners.size();
  for (int k = 0; k < nListeners; k++) {
    m_listeners[k]->iteration(*this, m_eta,
                              its, m_val_design, m_val_state,
                              m_d, m_d_diff, m_grad_design,
                              m_ssaforward.solution(),
                              m_u_diff,
                              m_grad_state,
                              m_constraints);
  }
}

void IP_SSATaucTaoTikhonovProblemLCL::evaluateObjectiveAndGradient(Tao /*tao*/, Vec x,
                                                                   double *value, Vec gradient) {

  m_x->scatter(x,m_dGlobal.vec(),m_uGlobal->vec());
  m_uGlobal->scale(m_velocityScale);

  // Variable 'm_dGlobal' has no ghosts.  We need ghosts for computation with the design variable.
  m_d->copy_from(m_dGlobal);

  m_d_diff->copy_from(*m_d);
  m_d_diff->add(-1,m_d0);
  m_designFunctional.gradientAt(*m_d_diff, *m_grad_design);
  m_grad_design->scale(1/m_eta);

  m_u_diff->copy_from(*m_uGlobal);
  m_u_diff->add(-1, m_u_obs);
  m_stateFunctional.gradientAt(*m_u_diff, *m_grad_state);
  m_grad_state->scale(m_velocityScale);

  m_x->gather(m_grad_design->vec(), m_grad_state->vec(), gradient);

  m_designFunctional.valueAt(*m_d_diff, &m_val_design);
  m_stateFunctional.valueAt(*m_u_diff, &m_val_state);

  *value = m_val_design / m_eta + m_val_state;
}

TerminationReason::Ptr IP_SSATaucTaoTikhonovProblemLCL::formInitialGuess(Vec *x) {
  m_d->copy_from(m_dGlobal);
  TerminationReason::Ptr reason = m_ssaforward.linearize_at(*m_d);
  if (reason->failed()) {
    return reason;
  }

  m_uGlobal->copy_from(*m_ssaforward.solution());
  m_uGlobal->scale(1.0 / m_velocityScale);

  m_x->gather(m_dGlobal.vec(), m_uGlobal->vec());

  // This is probably irrelevant.
  m_uGlobal->scale(m_velocityScale);

  *x =  *m_x;
  return GenericTerminationReason::success();
}

void IP_SSATaucTaoTikhonovProblemLCL::evaluateConstraints(Tao, Vec x, Vec r) {
  PetscErrorCode ierr;

  m_x->scatter(x,m_dGlobal.vec(),m_uGlobal->vec());
  m_uGlobal->scale(m_velocityScale);

  m_d->copy_from(m_dGlobal);
  m_u.copy_from(*m_uGlobal);

  m_ssaforward.set_design(*m_d);

  m_ssaforward.assemble_residual(m_u, r);

  ierr = VecScale(r,1./m_constraintsScale);
  PISM_CHK(ierr, "VecScale");
}

void IP_SSATaucTaoTikhonovProblemLCL::evaluateConstraintsJacobianState(Tao, Vec x,
                                                                       Mat Jstate,
                                                                       Mat /*Jpc*/,
                                                                       Mat /*Jinv*/,
                                                                       MatStructure *s) {
  PetscErrorCode ierr;

  m_x->scatter(x, m_dGlobal.vec(), m_uGlobal->vec());
  m_uGlobal->scale(m_velocityScale);

  m_d->copy_from(m_dGlobal);
  m_u.copy_from(*m_uGlobal);

  m_ssaforward.set_design(*m_d);

  m_ssaforward.assemble_jacobian_state(m_u, Jstate);

  *s = SAME_NONZERO_PATTERN;

  ierr = MatScale(Jstate, m_velocityScale / m_constraintsScale);
  PISM_CHK(ierr, "MatScale");
}

void IP_SSATaucTaoTikhonovProblemLCL::evaluateConstraintsJacobianDesign(Tao, Vec x, Mat /*Jdesign*/) {
  // I'm not sure if the following are necessary (i.e. will the copies that happen
  // in evaluateObjectiveAndGradient be sufficient) but we'll do them here
  // just in case.
  m_x->scatter(x,m_dGlobal.vec(),m_uGlobal->vec());
  m_uGlobal->scale(m_velocityScale);
  m_d_Jdesign.copy_from(m_dGlobal);
  m_u_Jdesign.copy_from(*m_uGlobal);
}

void IP_SSATaucTaoTikhonovProblemLCL::applyConstraintsJacobianDesign(Vec x, Vec y) {
  m_dzeta.copy_from_vec(x);

  m_ssaforward.set_design(m_d_Jdesign);

  m_ssaforward.apply_jacobian_design(m_u_Jdesign, m_dzeta, y);

  PetscErrorCode ierr = VecScale(y,1./m_constraintsScale);
  PISM_CHK(ierr, "VecScale");
}

void IP_SSATaucTaoTikhonovProblemLCL::applyConstraintsJacobianDesignTranspose(Vec x, Vec y) {
  m_du.copy_from_vec(x);

  m_ssaforward.set_design(m_d_Jdesign);

  m_ssaforward.apply_jacobian_design_transpose(m_u_Jdesign, m_du, y);

  PetscErrorCode ierr = VecScale(y, 1.0 / m_constraintsScale);
  PISM_CHK(ierr, "VecScale");
}

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL::jacobian_design_callback(Mat A, Vec x, Vec y) {
  try {
    IP_SSATaucTaoTikhonovProblemLCL *ctx;
    PetscErrorCode ierr = MatShellGetContext(A,&ctx);
    PISM_CHK(ierr, "MatShellGetContext");

    ctx->applyConstraintsJacobianDesign(x,y);
  } catch (...) {
    MPI_Comm com = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)A, &com); CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL::jacobian_design_transpose_callback(Mat A, Vec x,
                                                                                   Vec y) {
  try {
    IP_SSATaucTaoTikhonovProblemLCL *ctx;
    PetscErrorCode ierr = MatShellGetContext(A,&ctx);
    PISM_CHK(ierr, "MatShellGetContext");

    ctx->applyConstraintsJacobianDesignTranspose(x,y);
  } catch (...) {
    MPI_Comm com = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)A, &com); CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

} // end of namespace inverse
} // end of namespace pism
