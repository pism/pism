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
#include <memory>
#include "iceModelVec.hh"

class InvSSATikhonov;
// #include "TikhonovProblemListener.hh"


PetscErrorCode InvSSA_LCLTikhonov_applyJacobianDesign(Mat A, Vec x, Vec y);
PetscErrorCode InvSSA_LCLTikhonov_applyJacobianDesignTranspose(Mat A, Vec x, Vec y);

class InvSSA_LCLTikhonov {
public:
  
  typedef IceModelVec2S  DesignVec;
  typedef IceModelVec2V  StateVec;

  // typedef TikhonovProblemListener<InverseProblem> Listener;
  // typedef typename Listener::Ptr ListenerPtr;
  
  InvSSA_LCLTikhonov( InvSSATikhonov &invProblem, DesignVec &d0, StateVec &u_obs, PetscReal eta);

  virtual ~InvSSA_LCLTikhonov();

  // virtual void addListener(ListenerPtr listener) {
  //   m_listeners.push_back(listener);
  // }

  virtual StateVec &stateSolution();
  virtual DesignVec &designSolution();

  virtual InvSSATikhonov &invProblem() { return m_invProblem;}

  PetscErrorCode connect(TaoSolver tao);

  PetscErrorCode monitorTao(TaoSolver tao);
  //  {
  // //   PetscErrorCode ierr;
  // //   
  // //   PetscInt its;
  // //   ierr =  TaoGetSolutionStatus(tao, &its, NULL, NULL, NULL, NULL, NULL ); CHKERRQ(ierr);
  // //   
  // //   int nListeners = m_listeners.size();
  // //   for(int k=0; k<nListeners; k++) {
  // //    ierr = m_listeners[k]->iteration(*this,m_eta,
  // //                  its,m_valObjective,m_valPenalty,
  // //                  m_d, m_d_diff, m_grad_objective,
  // //                  m_invProblem.solution(), m_u_diff, m_grad_penalty,
  // //                  m_grad ); CHKERRQ(ierr);
  // //   }
  // //   return 0;
  // // }

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

  Vec m_constraints;
  Mat m_Jstate;
  Mat m_Jdesign;

  IceModelVec2S m_d_Jdesign;
  IceModelVec2V m_u_Jdesign;

  PetscReal m_constraintsScale;

  // std::vector<ListenerPtr> m_listeners;
};


#endif /* end of include guard: INVSSA_LCLTIKHONOV_HH_9T38Z13E */

