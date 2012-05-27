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

template<class InvProblem>
class TikhonovProblem {
public:
  typedef InvProblem::DomainVec DesignVec;
  typedef InvProblem::RangeVec  StateVec;
  
  TikhonovProblem( InvProblem &invProblem,
                   DomainVec &d0, StateVec &u_obs, PetscReal eta):
  m_d0(d0), m_u_obs(u_obs), m_eta(eta) {
    IceGrid &grid = *m_d0.get_grid();

    PetscInt design_stencil_width = m_d0.get_stencil_width();
    m_d.create(grid, "design variable", kHasGhosts, design_stencil_width);
    m_dGlobal.create(grid, "design variable (global)", kNoGhosts, design_stencil_width);
    m_dGlobal.copy(m_d);

    m_du.create( grid, "state residual", kHasGhosts, state_stencil_width);
    m_dd.create( grid, "design residual", kHasGhosts, design_stencil_width);

    m_grad_penalty.create( grid, "penalty gradient", kNoGhosts, 0);
    m_grad_objective.create( grid, "objective gradient", kNoGhosts, 0);
  }
  virtual ~TikhonovProblem() {};

  virtual InvProblem &invProblem() { return m_invProblem;}

  std::pair<DesignVec &,StateVec &> solution() {
    std::pair<DesignVec &,StateVec &> rv;
    rv.first = m_d;
    // FIXME: copy?
    rv.second = m_invProblem.solution();
    return rv;
  }

protected:
  virtual PetscErrorCode evaluateObjectiveAndGradient(TaoSolver tao, Vec x, PetscReal *value, Vec gradient) {
    PetscErrorCode ierr;
    
    // Variable 'x' has no ghosts.  We need ghosts for computation with the design variable.
    ierr = m_d.copy_from(x); CHKERRQ(ierr);

    bool success =  m_invProblem.linearizeAt(m_d);
    if(!success) {
      SETERRQ("Failure in TikhonovProblem forward solve: %s",m_invProblem.reason().c_str());
    }

    ierr = m_d_diff.copy_from(m_d); CHKERRQ(ierr);
    ierr = m_d_diff.add(-1,m_d0); CHKERRQ(ierr);
    ierr = m_invProblem.evalGradObjective(m_d_diff,m_grad_objective); CHKERRQ(ierr);

    ierr = m_u_diff.copy_from(m_invProblem.solution()); CHKERRQ(ierr);
    ierr = m_u_diff.add(-1, m_u_obs); CHKERRQ(ierr);
    ierr = m_invProblem.evalGradPenaltyReduced(m_u_diff,m_grad_penalty); CHKERRQ(ierr);
    ierr = m_grad_penalty.scale(m_eta); CHKERRQ(ierr);

    ierr = m_grad.copy_from(m_grad_penalty); CHKERRQ(ierr);
    ierr = m_grad.add(1,m_grad_objective); CHKERRQ(ierr);
    ierr = m_grad.copyTo(gradient); CHKERRQ(ierr);

    PetscReal valObjective, valPenalty;
    ierr = m_invProblem.evalObjective(m_d_diff,&valObjective); CHKERRQ(ierr);
    ierr = m_invProblem.evalPenalty(m_u_diff,&valPenalty); CHKERRQ(ierr);

    *value = valObjective + m_eta * valPenalty;
    
    return 0;
  }
  
  virtual Vec formInitialGuess() {
    return m_dGlobal.vec();
  }

  DesignVec m_dGlobal;
  DesignVec m_d;
  DesignVec &m_d0;
  DesignVec m_d_dff;

  StateVec &m_u_obs;
  StateVec m_u_diff;

  DesignVec m_grad_objective;
  DesignVec m_grad_penalty;
  DesignVec m_grad;

  PetscReal m_eta;

  InvProblem &m_invProblem;
};

#endif /* end of include guard: TIKHONOVPROBLEM_HH_70YYTOXB */
