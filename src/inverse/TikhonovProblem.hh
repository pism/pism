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

    // ierr = m_grad.set(0.); CHKERRQ(ierr);
    ierr = m_grad.copy_from(m_grad_objective); CHKERRQ(ierr);
    ierr = m_grad.scale(1./m_eta); CHKERRQ(ierr);    
    ierr = m_grad.add(1,m_grad_penalty); CHKERRQ(ierr);
    ierr = m_grad.copy_to(gradient); CHKERRQ(ierr);

    PetscReal valObjective, valPenalty;
    ierr = m_invProblem.evalObjective(m_d_diff,&valObjective); CHKERRQ(ierr);
    ierr = m_invProblem.evalPenalty(m_u_diff,&valPenalty); CHKERRQ(ierr);

    m_valObjective = valObjective;
    m_valPenalty = valPenalty;
    
    *value = valObjective / m_eta + valPenalty;

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

  MPI_Comm m_comm;
};

#endif /* end of include guard: TIKHONOVPROBLEM_HH_70YYTOXB */
