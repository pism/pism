// Copyright (C) 2012, 2014  David Maxwell and Constantine Khroulev
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
#include "iceModelVec.hh"
#include "IP_SSATaucForwardProblem.hh"
#include "functional/IPFunctional.hh"

namespace pism {

class IP_SSATaucTaoTikhonovProblemLCL;

//! Iteration callback class for IP_SSATaucTaoTikhonovProblemLCL
/*! A class for objects receiving iteration callbacks from a IP_SSATaucTaoTikhonovProblemLCL.  These 
  callbacks can be used to monitor the solution, plot iterations, print diagnostic messages, etc. 
  IP_SSATaucTaoTikhonovProblemLCLListeners are ususally used via a reference counted pointer 
  IP_SSATaucTaoTikhonovProblemLCLListeners::Ptr to allow for good memory management when Listeners are 
  created as subclasses of python classes.*/
class IP_SSATaucTaoTikhonovProblemLCLListener {
public:

#ifdef PISM_USE_TR1
  typedef std::tr1::shared_ptr<IP_SSATaucTaoTikhonovProblemLCLListener> Ptr;
#else
  typedef std::shared_ptr<IP_SSATaucTaoTikhonovProblemLCLListener> Ptr;
#endif


  typedef IceModelVec2S DesignVec;
  typedef IceModelVec2V StateVec;
  
  IP_SSATaucTaoTikhonovProblemLCLListener() {}
  virtual ~IP_SSATaucTaoTikhonovProblemLCLListener() {}
  
  //!Callback called after each iteration.
  virtual PetscErrorCode 
  iteration(IP_SSATaucTaoTikhonovProblemLCL &problem,  ///< The class calling the callback.
             double eta,                             ///< Tikhonov penalty parameter.
             int iter,                             ///< Current iteration count.
             double objectiveValue,                  ///< Value of the state functional.
             double designValue,                     ///< Value of the design functiona.
             DesignVec &d,                              ///< Value of the design variable.
             DesignVec &diff_d,                         ///< Diference between design variable and a-prior estimate.
             DesignVec &grad_d,                         ///< Gradient of design functional
             StateVec &u,                               ///< Value of state variable
             StateVec &diff_u,                          ///< Difference between state variable and desired value.
             StateVec &grad_u,                          ///< Gradient of state functional
             StateVec &constraints                      ///< Residual for state variable being a solution of the %SSA
             ) = 0;
};

PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL_applyJacobianDesign(Mat A, Vec x, Vec y);
PetscErrorCode IP_SSATaucTaoTikhonovProblemLCL_applyJacobianDesignTranspose(Mat A, Vec x, Vec y);

//! \brief Defines a Tikhonov minimization problem of determining \f$\tau_c\f$ from %SSA velocities to be solved with a TaoBasicSolver using the tao_lcl algorithm.
/*! Experimental and not particularly functional. */
class IP_SSATaucTaoTikhonovProblemLCL {
public:
  
  typedef IceModelVec2S  DesignVec;
  typedef IceModelVec2V  StateVec;

  typedef IP_SSATaucTaoTikhonovProblemLCLListener Listener;
  
  IP_SSATaucTaoTikhonovProblemLCL(IP_SSATaucForwardProblem &ssaforward, DesignVec &d0, StateVec &u_obs, double eta,
                                   IPFunctional<DesignVec> &designFunctional, IPFunctional<StateVec> &stateFunctional);

  virtual ~IP_SSATaucTaoTikhonovProblemLCL();

  virtual void addListener(Listener::Ptr listener) {
    m_listeners.push_back(listener);
  }

  virtual StateVec &stateSolution();
  virtual DesignVec &designSolution();

  virtual PetscErrorCode setInitialGuess(DesignVec &d0);

  PetscErrorCode connect(Tao tao);

  PetscErrorCode monitorTao(Tao tao);

  virtual PetscErrorCode evaluateObjectiveAndGradient(Tao tao, Vec x, double *value, Vec gradient);
  
  virtual PetscErrorCode formInitialGuess(Vec *x,TerminationReason::Ptr &reason);

  virtual PetscErrorCode evaluateConstraints(Tao, Vec x, Vec r);

  virtual PetscErrorCode evaluateConstraintsJacobianState(Tao, Vec x, Mat Jstate, Mat Jpc, Mat Jinv, MatStructure *s);
  
  virtual PetscErrorCode evaluateConstraintsJacobianDesign(Tao, Vec x, Mat Jdesign);

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

  double m_eta;

  double m_val_design;
  double m_val_state;

  StateVec m_constraints;
  Mat m_Jstate;
  Mat m_Jdesign;

  IceModelVec2S m_d_Jdesign;
  IceModelVec2V m_u_Jdesign;

  double m_constraintsScale;
  double m_velocityScale;

  IPFunctional<IceModelVec2S> &m_designFunctional;
  IPFunctional<IceModelVec2V> &m_stateFunctional;

  std::vector<Listener::Ptr> m_listeners;
};

} // end of namespace pism

#endif /* end of include guard: IP_SSATAUCTIKHONOVLCL_HH_39UGM4S2 */
