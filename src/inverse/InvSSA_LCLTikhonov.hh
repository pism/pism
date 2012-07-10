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

#ifndef INVSSA_LCLTIKHONOV_HH_9T38Z13E
#define INVSSA_LCLTIKHONOV_HH_9T38Z13E


#include "TaoUtil.hh"
#include "TwoBlockVec.hh"
#include <petsc.h>
#include <tr1/memory>
#include "iceModelVec.hh"
#include "PythonTikhonovSVListener.hh"

class InvSSATikhonov;
// #include "TikhonovProblemListener.hh"


class InvSSA_LCLTikhonov;
class InvSSA_LCLTikhonovProblemListener {
public:
  typedef std::tr1::shared_ptr<InvSSA_LCLTikhonovProblemListener> Ptr;
  typedef IceModelVec2S DesignVec;
  typedef IceModelVec2V StateVec;
  
  InvSSA_LCLTikhonovProblemListener() {}
  virtual ~InvSSA_LCLTikhonovProblemListener() {}
  
  virtual PetscErrorCode 
  iteration( InvSSA_LCLTikhonov &problem,
             PetscReal eta, PetscInt iter,
             PetscReal objectiveValue, PetscReal designValue,
             DesignVec &d, DesignVec &diff_d, DesignVec &grad_d,
             StateVec &u,   StateVec &diff_u,  StateVec &grad_u,
             StateVec &constraints) {
               (void) problem; (void) eta; (void) d; 
               (void) diff_d; (void) grad_d; (void) u;
               (void) diff_u; (void) grad_u; (void) constraints;
               printf("Iteration %d: objValue %g designValue %g\n",iter,objectiveValue,designValue);
               return 0;};
};


class InvSSA_LCLTikhonovPythonListenerBridge: public InvSSA_LCLTikhonovProblemListener {
public:
  InvSSA_LCLTikhonovPythonListenerBridge(PythonLCLTikhonovSVListener::Ptr core) : m_core(core) { }
  PetscErrorCode iteration( InvSSA_LCLTikhonov &problem,
             PetscReal eta, PetscInt iter,
             PetscReal objectiveValue, PetscReal designValue,
             DesignVec &d, DesignVec &diff_d, DesignVec &grad_d,
             StateVec &u,   StateVec &diff_u,  StateVec &grad_u,
             StateVec &constraints) { 
    (void) problem;
    m_core->iteration(iter,eta,objectiveValue,designValue,d,diff_d,grad_d,
      u,diff_u,grad_u,constraints);
    return 0;
  }
protected:
  PythonLCLTikhonovSVListener::Ptr m_core;
};

void InvSSALCLTikhonovAddListener(InvSSA_LCLTikhonov &problem, 
PythonLCLTikhonovSVListener::Ptr listener );

PetscErrorCode InvSSA_LCLTikhonov_applyJacobianDesign(Mat A, Vec x, Vec y);
PetscErrorCode InvSSA_LCLTikhonov_applyJacobianDesignTranspose(Mat A, Vec x, Vec y);

class InvSSA_LCLTikhonov {
public:
  
  typedef IceModelVec2S  DesignVec;
  typedef IceModelVec2V  StateVec;

  typedef InvSSA_LCLTikhonovProblemListener Listener;
  
  InvSSA_LCLTikhonov( InvSSATikhonov &invProblem, DesignVec &d0, StateVec &u_obs, PetscReal eta);

  virtual ~InvSSA_LCLTikhonov();

  virtual void addListener(Listener::Ptr listener) {
    m_listeners.push_back(listener);
  }

  virtual StateVec &stateSolution();
  virtual DesignVec &designSolution();

  virtual InvSSATikhonov &invProblem() { return m_invProblem;}

  virtual PetscErrorCode setInitialGuess( DesignVec &d0);

  PetscErrorCode connect(TaoSolver tao);

  PetscErrorCode monitorTao(TaoSolver tao);

  virtual PetscErrorCode evaluateObjectiveAndGradient(TaoSolver tao, Vec x, PetscReal *value, Vec gradient);
  
  virtual Vec formInitialGuess();

  virtual PetscErrorCode evaluateConstraints(TaoSolver, Vec x, Vec r);

  virtual PetscErrorCode evaluateConstraintsJacobianState(TaoSolver, Vec x, Mat *Jstate, Mat *Jpc, Mat *Jinv, MatStructure *s);
  
  virtual PetscErrorCode evaluateConstraintsJacobianDesign(TaoSolver, Vec x, Mat* /*Jdesign*/);

  virtual PetscErrorCode applyConstraintsJacobianDesign(Vec x, Vec y);
  virtual PetscErrorCode applyConstraintsJacobianDesignTranspose(Vec x, Vec y);

protected:

  virtual PetscErrorCode construct();
  virtual PetscErrorCode destruct();

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

  DesignVec m_grad_objective;
  StateVec  m_grad_penalty;

  PetscReal m_eta;

  PetscReal m_valObjective;
  PetscReal m_valPenalty;

  InvSSATikhonov &m_invProblem;

  DM m_da;

  StateVec m_constraints;
  Mat m_Jstate;
  Mat m_Jdesign;

  IceModelVec2S m_d_Jdesign;
  IceModelVec2V m_u_Jdesign;

  PetscReal m_constraintsScale;
  PetscReal m_velocityScale;

  std::vector<Listener::Ptr> m_listeners;
};


#endif /* end of include guard: INVSSA_LCLTIKHONOV_HH_9T38Z13E */
