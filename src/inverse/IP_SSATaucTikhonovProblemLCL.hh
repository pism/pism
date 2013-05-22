// Copyright (C) 2012  David Maxwell
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

#ifndef IP_SSATAUCTIKHONOVLCL_HH_39UGM4S2
#define IP_SSATAUCTIKHONOVLCL_HH_39UGM4S2


#include "TaoUtil.hh"
#include "IPTwoBlockVec.hh"
#include <petsc.h>
#include <tr1/memory>
#include "iceModelVec.hh"
#include "IP_SSATaucForwardProblem.hh"
#include "functional/IPFunctional.hh"

class IP_SSATaucTikhonovProblemLCL;
class IP_SSATaucTikhonovProblemLCLListener {
public:
  typedef std::tr1::shared_ptr<IP_SSATaucTikhonovProblemLCLListener> Ptr;
  typedef IceModelVec2S DesignVec;
  typedef IceModelVec2V StateVec;
  
  IP_SSATaucTikhonovProblemLCLListener() {}
  virtual ~IP_SSATaucTikhonovProblemLCLListener() {}
  
  virtual PetscErrorCode 
  iteration( IP_SSATaucTikhonovProblemLCL &problem,
             PetscReal eta, PetscInt iter,
             PetscReal objectiveValue, PetscReal designValue,
             DesignVec &d, DesignVec &diff_d, DesignVec &grad_d,
             StateVec &u,   StateVec &diff_u,  StateVec &grad_u,
           StateVec &constraints) = 0;
};

PetscErrorCode IP_SSATaucTikhonovProblemLCL_applyJacobianDesign(Mat A, Vec x, Vec y);
PetscErrorCode IP_SSATaucTikhonovProblemLCL_applyJacobianDesignTranspose(Mat A, Vec x, Vec y);

class IP_SSATaucTikhonovProblemLCL {
public:
  
  typedef IceModelVec2S  DesignVec;
  typedef IceModelVec2V  StateVec;

  typedef IP_SSATaucTikhonovProblemLCLListener Listener;
  
  IP_SSATaucTikhonovProblemLCL( IP_SSATaucForwardProblem &ssaforward, DesignVec &d0, StateVec &u_obs, PetscReal eta,
                      IPFunctional<DesignVec> &designFunctional, IPFunctional<StateVec> &stateFunctional);

  virtual ~IP_SSATaucTikhonovProblemLCL();

  virtual void addListener(Listener::Ptr listener) {
    m_listeners.push_back(listener);
  }

  virtual StateVec &stateSolution();
  virtual DesignVec &designSolution();

  virtual PetscErrorCode setInitialGuess( DesignVec &d0);

  PetscErrorCode connect(TaoSolver tao);

  PetscErrorCode monitorTao(TaoSolver tao);

  virtual PetscErrorCode evaluateObjectiveAndGradient(TaoSolver tao, Vec x, PetscReal *value, Vec gradient);
  
  virtual PetscErrorCode formInitialGuess(Vec *x,TerminationReason::Ptr &reason);

  virtual PetscErrorCode evaluateConstraints(TaoSolver, Vec x, Vec r);

  virtual PetscErrorCode evaluateConstraintsJacobianState(TaoSolver, Vec x, Mat *Jstate, Mat *Jpc, Mat *Jinv, MatStructure *s);
  
  virtual PetscErrorCode evaluateConstraintsJacobianDesign(TaoSolver, Vec x, Mat* /*Jdesign*/);

  virtual PetscErrorCode applyConstraintsJacobianDesign(Vec x, Vec y);
  virtual PetscErrorCode applyConstraintsJacobianDesignTranspose(Vec x, Vec y);

protected:

  virtual PetscErrorCode construct();
  virtual PetscErrorCode destruct();

  IP_SSATaucForwardProblem &m_ssaforward;

  std::auto_ptr<IPTwoBlockVec> m_x;

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

  IPFunctional<IceModelVec2S> &m_designFunctional;
  IPFunctional<IceModelVec2V> &m_stateFunctional;

  std::vector<Listener::Ptr> m_listeners;
};


#endif /* end of include guard: IP_SSATAUCTIKHONOVLCL_HH_39UGM4S2 */
