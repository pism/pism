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

#ifndef TIKHONOVPROBLEM_HH_70YYTOXB
#define TIKHONOVPROBLEM_HH_70YYTOXB

#include "TaoUtil.hh"
#include <petsc.h>

#include "TikhonovProblemListener.hh"

template<class InverseProblem>
class TikhonovProblem {
public:
  
  typedef typename InverseProblem::DesignVec DesignVec;
  typedef typename InverseProblem::StateVec StateVec;
  typedef TikhonovProblemListener<InverseProblem> Listener;
  typedef typename Listener::Ptr ListenerPtr;
  
  TikhonovProblem( InverseProblem &inv_problem,
                   DesignVec &d0, StateVec &u_obs, PetscReal eta):
  m_invProblem(inv_problem), m_d0(d0), m_u_obs(u_obs), m_eta(eta) {
    IceGrid &grid = *m_d0.get_grid();
    m_comm = grid.com;

    m_tikhonov_atol = grid.config.get("tikhonov_atol");
    m_tikhonov_rtol = grid.config.get("tikhonov_rtol");

    PetscInt design_stencil_width = m_d0.get_stencil_width();
    PetscInt state_stencil_width = m_u_obs.get_stencil_width();
    m_d.create(grid, "design variable", kHasGhosts, design_stencil_width);
    m_dGlobal.create(grid, "design variable (global)", kNoGhosts, design_stencil_width);
    m_dGlobal.copy_from(m_d0);

    m_u_diff.create( grid, "state residual", kHasGhosts, state_stencil_width);
    m_d_diff.create( grid, "design residual", kHasGhosts, design_stencil_width);

    m_grad_penalty.create( grid, "penalty gradient", kNoGhosts, design_stencil_width);
    m_grad_objective.create( grid, "objective gradient", kNoGhosts, design_stencil_width);
    m_grad.create( grid, "gradient", kNoGhosts, design_stencil_width);

  }
  virtual ~TikhonovProblem() {};

  virtual InverseProblem &invProblem() { return m_invProblem;}

  virtual void addListener(ListenerPtr listener) {
    m_listeners.push_back(listener);
  }

  virtual StateVec &stateSolution() {
    m_invProblem.solution();
  }

  virtual DesignVec &designSolution() {
    return m_d;
  }

  PetscErrorCode connect(TaoSolver tao) {
    PetscErrorCode ierr;
    ierr = TaoObjGradCallback<TikhonovProblem>::connect(tao,*this); CHKERRQ(ierr);
    ierr = TaoMonitorCallback<TikhonovProblem>::connect(tao,*this); CHKERRQ(ierr);
    ierr = TaoConvergenceCallback<TikhonovProblem>::connect(tao,*this); CHKERRQ(ierr);
    ierr = TaoGetVariableBoundsCallback<TikhonovProblem>::connect(tao,*this); CHKERRQ(ierr);
    return 0;
  }

  PetscErrorCode monitorTao(TaoSolver tao) {
    PetscErrorCode ierr;
    
    PetscInt its;
    ierr =  TaoGetSolutionStatus(tao, &its, NULL, NULL, NULL, NULL, NULL ); CHKERRQ(ierr);
    
    int nListeners = m_listeners.size();
    for(int k=0; k<nListeners; k++) {
     ierr = m_listeners[k]->iteration(*this,m_eta,
                   its,m_valObjective,m_valPenalty,
                   m_d, m_d_diff, m_grad_objective,
                   m_invProblem.solution(), m_u_diff, m_grad_penalty,
                   m_grad ); CHKERRQ(ierr);
    }
    return 0;
  }

  virtual PetscErrorCode convergenceTest(TaoSolver tao) {
    PetscErrorCode ierr;
    
    PetscReal designNorm, stateNorm, sumNorm;
    PetscReal dWeight, sWeight;
    if(m_eta>1) {
      dWeight = 1/m_eta;
      sWeight = 1;
    } else {
      dWeight = 1;
      sWeight = m_eta;      
    }
    
    ierr = m_grad_objective.norm(NORM_2,designNorm); CHKERRQ(ierr);
    ierr = m_grad_penalty.norm(NORM_2,stateNorm); CHKERRQ(ierr);
    ierr = m_grad.norm(NORM_2,sumNorm); CHKERRQ(ierr);
    designNorm *= dWeight;    
    stateNorm  *= sWeight;
    
    if( sumNorm < m_tikhonov_atol && sumNorm < m_tikhonov_rtol*PetscMax(designNorm,stateNorm) ) {
      ierr = TaoSetTerminationReason(tao,TAO_CONVERGED_USER); CHKERRQ(ierr);
    } else {
      ierr = TaoDefaultConvergenceTest(tao,NULL); CHKERRQ(ierr);
    }

    return 0;
  }

  virtual PetscErrorCode getVariableBounds(TaoSolver /*tao*/, Vec lo, Vec hi) {
    PetscErrorCode ierr;
    ierr = m_invProblem.getVariableBounds(lo,hi); CHKERRQ(ierr);
    return 0;
  }

  virtual PetscErrorCode evaluateObjective(TaoSolver tao, Vec x, PetscReal *value) {
    PetscErrorCode ierr;

    // Variable 'x' has no ghosts.  We need ghosts for computation with the design variable.
    ierr = m_d.copy_from(x); CHKERRQ(ierr);

    bool success;
    ierr = m_invProblem.linearizeAt(m_d, success); CHKERRQ(ierr);
    if(!success) {
      SETERRQ1(m_comm,1,"Failure in TikhonovProblem forward solve. %s",m_invProblem.reasonDescription().c_str());
    }

    ierr = m_d_diff.copy_from(m_d); CHKERRQ(ierr);
    ierr = m_d_diff.add(-1,m_d0); CHKERRQ(ierr);

    ierr = m_u_diff.copy_from(m_invProblem.solution()); CHKERRQ(ierr);
    ierr = m_u_diff.add(-1, m_u_obs); CHKERRQ(ierr);

    PetscReal valObjective, valPenalty;
    ierr = m_invProblem.evalObjective(m_d_diff,&valObjective); CHKERRQ(ierr);
    ierr = m_invProblem.evalPenalty(m_u_diff,&valPenalty); CHKERRQ(ierr);

    m_valObjective = valObjective;
    m_valPenalty = valPenalty;
    
    *value = valObjective / m_eta + valPenalty;
    return 0;
  }

  virtual PetscErrorCode evaluateObjectiveAndGradient(TaoSolver tao, Vec x, PetscReal *value, Vec gradient) {
    PetscErrorCode ierr;

    // Variable 'x' has no ghosts.  We need ghosts for computation with the design variable.
    ierr = m_d.copy_from(x); CHKERRQ(ierr);

    bool success;
    ierr = m_invProblem.linearizeAt(m_d, success); CHKERRQ(ierr);
    if(!success) {
      SETERRQ1(m_comm,1,"Failure in TikhonovProblem forward solve. %s",m_invProblem.reasonDescription().c_str());
    }

    ierr = m_d_diff.copy_from(m_d); CHKERRQ(ierr);
    ierr = m_d_diff.add(-1,m_d0); CHKERRQ(ierr);
    ierr = m_invProblem.evalGradObjective(m_d_diff,m_grad_objective); CHKERRQ(ierr);

    ierr = m_u_diff.copy_from(m_invProblem.solution()); CHKERRQ(ierr);
    ierr = m_u_diff.add(-1, m_u_obs); CHKERRQ(ierr);
    ierr = m_invProblem.evalGradPenaltyReduced(m_u_diff,m_grad_penalty); CHKERRQ(ierr);
    // ierr = evaluateGradientReducedFD2(m_invProblem,m_d,m_u_obs, m_grad_penalty); CHKERRQ(ierr);


    if(m_eta>1) {
      ierr = m_grad.copy_from(m_grad_objective); CHKERRQ(ierr);
      ierr = m_grad.scale(1./m_eta); CHKERRQ(ierr);    
      ierr = m_grad.add(1,m_grad_penalty); CHKERRQ(ierr);
    } else{
      ierr = m_grad.copy_from(m_grad_penalty); CHKERRQ(ierr);
      ierr = m_grad.scale(m_eta); CHKERRQ(ierr);    
      ierr = m_grad.add(1,m_grad_objective); CHKERRQ(ierr);
    }
    ierr = m_grad.copy_to(gradient); CHKERRQ(ierr);      

    PetscReal valObjective, valPenalty;
    ierr = m_invProblem.evalObjective(m_d_diff,&valObjective); CHKERRQ(ierr);
    ierr = m_invProblem.evalPenalty(m_u_diff,&valPenalty); CHKERRQ(ierr);

    m_valObjective = valObjective;
    m_valPenalty = valPenalty;
    
    if(m_eta > 1) {
      *value = valObjective / m_eta + valPenalty;
    } else {
      *value = valObjective + valPenalty*m_eta;
    }

    return 0;
  }
  
  virtual Vec formInitialGuess() {
    return m_dGlobal.get_vec();
  }
protected:

  DesignVec m_dGlobal;
  DesignVec m_d;
  DesignVec &m_d0;
  DesignVec m_d_diff;

  StateVec &m_u_obs;
  StateVec m_u_diff;

  DesignVec m_grad_objective;
  DesignVec m_grad_penalty;
  DesignVec m_grad;

  PetscReal m_eta;

  PetscReal m_valObjective;
  PetscReal m_valPenalty;

  InverseProblem &m_invProblem;

  std::vector<ListenerPtr> m_listeners;

  PetscReal m_tikhonov_atol;
  PetscReal m_tikhonov_rtol;

  MPI_Comm m_comm;

private:
  // Hide copy/assignment operations
  TikhonovProblem(TikhonovProblem const &);
  TikhonovProblem & operator=(TikhonovProblem const &);

};

template<class InverseProblem>
PetscErrorCode evaluateGradientReducedFD(InverseProblem &prob,IceModelVec2S &zeta,IceModelVec2V &u0, IceModelVec2S &gradient) {
  PetscErrorCode ierr;
  bool success;
  PetscReal h = PETSC_SQRT_MACHINE_EPSILON;

  IceGrid &grid = *zeta.get_grid();

  IceModelVec2V du;
  ierr = du.create(grid,"du",kHasGhosts,1); CHKERRQ(ierr);

  PetscReal f0;
  ierr = prob.linearizeAt(zeta,success); CHKERRQ(ierr);
  ierr = du.copy_from(prob.solution()); CHKERRQ(ierr);
  ierr = du.add(-1,u0); CHKERRQ(ierr);
  ierr = prob.evalPenalty(du, &f0); CHKERRQ(ierr);
  
  ierr = gradient.begin_access(); CHKERRQ(ierr);
  for(PetscInt i=grid.xs; i< grid.xs+grid.xm; i++) {
    for(PetscInt j=grid.ys; j< grid.ys+grid.ym; j++) {
      ierr = zeta.begin_access(); CHKERRQ(ierr);
      zeta(i,j) += h;
      ierr = zeta.end_access(); CHKERRQ(ierr);
      ierr = zeta.beginGhostComm(); CHKERRQ(ierr);
      ierr = zeta.endGhostComm(); CHKERRQ(ierr);
      
      ierr = prob.linearizeAt(zeta,success); CHKERRQ(ierr);

      ierr = zeta.begin_access(); CHKERRQ(ierr);
      zeta(i,j) -= h;
      ierr = zeta.end_access(); CHKERRQ(ierr);
      ierr = zeta.beginGhostComm(); CHKERRQ(ierr);
      ierr = zeta.endGhostComm(); CHKERRQ(ierr);

      ierr = du.copy_from(prob.solution()); CHKERRQ(ierr);
      ierr = du.add(-1,u0); CHKERRQ(ierr);
      PetscReal fh;
      ierr = prob.evalPenalty(du, &fh); CHKERRQ(ierr);
      gradient(i,j) = (fh-f0)/h;
    }
  }
  ierr = gradient.end_access(); CHKERRQ(ierr);
  ierr = prob.linearizeAt(zeta, success); CHKERRQ(ierr);
  return 0;
}

template<class InverseProblem>
PetscErrorCode evaluateGradientReducedFD2(InverseProblem &prob,IceModelVec2S &zeta,IceModelVec2V u0, IceModelVec2S gradient) {
  PetscErrorCode ierr;
  bool success;
  PetscReal h = PETSC_SQRT_MACHINE_EPSILON;

  IceGrid &grid = *zeta.get_grid();

  IceModelVec2V u;
  ierr = u.create(grid,"u",kHasGhosts,1); CHKERRQ(ierr);
  ierr = u.copy_from(prob.solution()); CHKERRQ(ierr);

  IceModelVec2V du;
  ierr = du.create(grid,"du",kHasGhosts,1); CHKERRQ(ierr);
  ierr = du.copy_from(u); CHKERRQ(ierr);
  ierr = du.add(-1,u0); CHKERRQ(ierr);

  IceModelVec2V gradPenalty;
  ierr = gradPenalty.create(grid,"gradPenalty",kNoGhosts,0); CHKERRQ(ierr);
  ierr = prob.evalGradPenalty(du,gradPenalty); CHKERRQ(ierr);
  
  ierr = gradient.begin_access(); CHKERRQ(ierr);
  for(PetscInt i=grid.xs; i< grid.xs+grid.xm; i++) {
    for(PetscInt j=grid.ys; j< grid.ys+grid.ym; j++) {
      ierr = zeta.begin_access(); CHKERRQ(ierr);
      zeta(i,j) += h;
      ierr = zeta.end_access(); CHKERRQ(ierr);
      ierr = zeta.beginGhostComm(); CHKERRQ(ierr);
      ierr = zeta.endGhostComm(); CHKERRQ(ierr);
      
      ierr = prob.linearizeAt(zeta,success); CHKERRQ(ierr);

      ierr = zeta.begin_access(); CHKERRQ(ierr);
      zeta(i,j) -= h;
      ierr = zeta.end_access(); CHKERRQ(ierr);
      ierr = zeta.beginGhostComm(); CHKERRQ(ierr);
      ierr = zeta.endGhostComm(); CHKERRQ(ierr);

      ierr = du.copy_from(prob.solution()); CHKERRQ(ierr);
      ierr = du.add(-1,u); CHKERRQ(ierr);
      ierr = du.scale(1/h); CHKERRQ(ierr);

      PetscReal g=0;
      ierr = du.begin_access();
      ierr = gradPenalty.begin_access();
      for(PetscInt k=grid.xs; k< grid.xs+grid.xm; k++) {
        for(PetscInt l=grid.ys; l< grid.ys+grid.ym; l++) {
          g += gradPenalty(k,l).u*du(k,l).u+gradPenalty(k,l).v*du(k,l).v;
        }
      }
      gradient(i,j) = g;
    }
  }
  ierr = gradient.end_access(); CHKERRQ(ierr);
  ierr = prob.linearizeAt(zeta, success); CHKERRQ(ierr);
  return 0;
}

#endif /* end of include guard: TIKHONOVPROBLEM_HH_70YYTOXB */
