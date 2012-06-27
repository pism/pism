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

#ifndef INVSSATIKHONOVLCL_HH_39UGM4S2
#define INVSSATIKHONOVLCL_HH_39UGM4S2


#include "TaoUtil.hh"
#include "TwoBlockVec.hh"
#include <petsc.h>
#include <tr1/memory>
#include "iceModelVec.hh"
#include "InvSSAForwardProblem.hh"
#include "Functional.hh"

class InvSSATikhonovLCL;
class InvSSATikhonovLCLListener {
public:
  typedef std::tr1::shared_ptr<InvSSATikhonovLCLListener> Ptr;
  typedef IceModelVec2S DesignVec;
  typedef IceModelVec2V StateVec;
  
  InvSSATikhonovLCLListener() {}
  virtual ~InvSSATikhonovLCLListener() {}
  
  virtual PetscErrorCode 
  iteration( InvSSATikhonovLCL &problem,
             PetscReal eta, PetscInt iter,
             PetscReal objectiveValue, PetscReal designValue,
             DesignVec &d, DesignVec &diff_d, DesignVec &grad_d,
             StateVec &u,   StateVec &diff_u,  StateVec &grad_u,
           StateVec &constraints) = 0;
};

PetscErrorCode InvSSATikhonovLCL_applyJacobianDesign(Mat A, Vec x, Vec y);
PetscErrorCode InvSSATikhonovLCL_applyJacobianDesignTranspose(Mat A, Vec x, Vec y);

class InvSSATikhonovLCL {
public:
  
  typedef IceModelVec2S  DesignVec;
  typedef IceModelVec2V  StateVec;

  typedef InvSSATikhonovLCLListener Listener;
  
  InvSSATikhonovLCL( InvSSAForwardProblem &ssaforward, DesignVec &d0, StateVec &u_obs, PetscReal eta,
                      Functional<DesignVec> &designFunctional, Functional<StateVec> &stateFunctional);

  virtual ~InvSSATikhonovLCL();

  virtual void addListener(Listener::Ptr listener) {
    m_listeners.push_back(listener);
  }

  virtual StateVec &stateSolution();
  virtual DesignVec &designSolution();

  virtual PetscErrorCode setInitialGuess( DesignVec &d0);

  PetscErrorCode connect(TaoSolver tao);

  PetscErrorCode monitorTao(TaoSolver tao);

  virtual PetscErrorCode evaluateObjectiveAndGradient(TaoSolver tao, Vec x, PetscReal *value, Vec gradient);
  
  virtual PetscErrorCode formInitialGuess(Vec *x);

  virtual PetscErrorCode evaluateConstraints(TaoSolver, Vec x, Vec r);

  virtual PetscErrorCode evaluateConstraintsJacobianState(TaoSolver, Vec x, Mat *Jstate, Mat *Jpc, Mat *Jinv, MatStructure *s);
  
  virtual PetscErrorCode evaluateConstraintsJacobianDesign(TaoSolver, Vec x, Mat* /*Jdesign*/);

  virtual PetscErrorCode applyConstraintsJacobianDesign(Vec x, Vec y);
  virtual PetscErrorCode applyConstraintsJacobianDesignTranspose(Vec x, Vec y);

protected:

  virtual PetscErrorCode construct();
  virtual PetscErrorCode destruct();

  InvSSAForwardProblem &m_ssaforward;

  std::auto_ptr<TwoBlockVec> m_x;

  DesignVec m_dGlobal;
  DesignVec m_d;
  DesignVec &m_d0;
  DesignVec m_d_diff;
  DesignVec m_dzeta;

  StateVec m_uGlobal;
  StateVec m_u;
  StateVec m_du;
  StateVec &m_u_obs;
  StateVec m_u_diff;

  DesignVec m_grad_design;
  StateVec  m_grad_state;

  PetscReal m_eta;

  PetscReal m_val_design;
  PetscReal m_val_state;

  StateVec m_constraints;
  Mat m_Jstate;
  Mat m_Jdesign;

  IceModelVec2S m_d_Jdesign;
  IceModelVec2V m_u_Jdesign;

  PetscReal m_constraintsScale;
  PetscReal m_velocityScale;

  Functional<IceModelVec2S> &m_designFunctional;
  Functional<IceModelVec2V> &m_stateFunctional;

  std::vector<Listener::Ptr> m_listeners;
};


#endif /* end of include guard: INVSSATIKHONOVLCL_HH_39UGM4S2 */
