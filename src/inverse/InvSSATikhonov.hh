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

#ifndef INVSSATIKHONOV_HH_GHK50ADU
#define INVSSATIKHONOV_HH_GHK50ADU

#include <tr1/memory>

#include "SSAForwardProblem.hh"
#include "TaoUtil.hh"
#include "Functional.hh"

class InvSSATikhonov;

class InvSSATikhonovListener {
public:
  typedef std::tr1::shared_ptr<InvSSATikhonovListener> Ptr;
  typedef IceModelVec2S DesignVec;
  typedef IceModelVec2V StateVec;
  
  InvSSATikhonovListener() {}
  virtual ~InvSSATikhonovListener() {}
  
  virtual PetscErrorCode 
  iteration( InvSSATikhonov &problem,
             PetscReal eta, PetscInt iter,
             PetscReal objectiveValue, PetscReal designValue,
             DesignVec &d, DesignVec &diff_d, DesignVec &grad_d,
             StateVec &u,  StateVec &diff_u,  DesignVec &grad_u,
             DesignVec &gradient) = 0;
};

class InvSSATikhonov
{
public:

  typedef IceModelVec2S DesignVec;
  typedef IceModelVec2V StateVec;
  typedef InvSSATikhonovListener Listener;

  InvSSATikhonov( SSAForwardProblem &ssaforward, DesignVec &d0, StateVec &u_obs, PetscReal eta, 
                  Functional<DesignVec> &designFunctional, Functional<StateVec> &stateFunctional);

  virtual ~InvSSATikhonov();

  virtual PetscErrorCode setInitialGuess( DesignVec &d) {
    PetscErrorCode ierr;
    ierr = m_dGlobal.copy_from(d); CHKERRQ(ierr);
    return 0;
  }

  virtual PetscErrorCode evaluateObjectiveAndGradient(TaoSolver tao, Vec x, PetscReal *value, Vec gradient);

  virtual void addListener(Listener::Ptr listener) {
    m_listeners.push_back(listener);
  }

  virtual StateVec &stateSolution() {
    return m_ssaforward.solution();
  }
  
  virtual DesignVec &designSolution() {
    return m_d;
  }

  virtual PetscErrorCode connect(TaoSolver tao);

  virtual PetscErrorCode monitorTao(TaoSolver tao);

  virtual PetscErrorCode convergenceTest(TaoSolver tao); 

  virtual PetscErrorCode getVariableBounds(TaoSolver /*tao*/, Vec lo, Vec hi); 

  virtual PetscErrorCode formInitialGuess(Vec *v) {
    *v = m_dGlobal.get_vec();
    return 0;
  }

protected:

  virtual PetscErrorCode construct();
  
  IceGrid *m_grid;
  
  SSAForwardProblem &m_ssaforward;

  DesignVec m_dGlobal;
  DesignVec m_d;
  DesignVec &m_d0;
  DesignVec m_d_diff;

  StateVec &m_u_obs;
  StateVec m_u_diff;

  StateVec m_adjointRHS;

  DesignVec m_grad_design;
  DesignVec m_grad_state;
  DesignVec m_grad;

  PetscReal m_eta;

  PetscReal m_val_design;
  PetscReal m_val_state;

  Functional<IceModelVec2S> &m_designFunctional;
  Functional<IceModelVec2V> &m_stateFunctional;

  std::vector<Listener::Ptr> m_listeners;

  PetscReal m_tikhonov_atol;
  PetscReal m_tikhonov_rtol;

};

#endif /* end of include guard: INVSSATIKHONOV_HH_GHK50ADU */
