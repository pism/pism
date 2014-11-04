// Copyright (C) 2012,2013,2014  David Maxwell
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


#ifndef IPTAOTIKHONOVPROBLEM_HH_4NMM724B
#define IPTAOTIKHONOVPROBLEM_HH_4NMM724B

#ifdef PISM_USE_TR1
#include <tr1/memory>
#else
#include <memory>
#endif

#include "TaoUtil.hh"
#include "functional/IPFunctional.hh"
#include <assert.h>
#include "PISMConfig.hh"

template<class ForwardProblem> class IPTaoTikhonovProblem;

//! Iteration callback class for IPTaoTikhonovProblem
/*! A class for objects receiving iteration callbacks from a IPTaoTikhonovProblem.  These 
    callbacks can be used to monitor the solution, plot iterations, print diagnostic messages, etc. 
    IPTaoTikhonovProblemListeners are ususally used via a reference counted pointer 
    IPTaoTikhonovProblemListener::Ptr to allow for good memory management when Listeners are 
    created as subclasses of python classes. It would have been better to nest this inside of 
    IPTaoTikhonovProblem, but SWIG has a hard time with nested classes, so it's outer instead.*/
template<class ForwardProblem> class IPTaoTikhonovProblemListener {
public:
#ifdef PISM_USE_TR1
  typedef std::tr1::shared_ptr<IPTaoTikhonovProblemListener> Ptr;
#else
  typedef std::shared_ptr<IPTaoTikhonovProblemListener> Ptr;
#endif

  typedef typename ForwardProblem::DesignVec DesignVec;
  typedef typename ForwardProblem::StateVec StateVec;

  IPTaoTikhonovProblemListener() {}
  virtual ~IPTaoTikhonovProblemListener() {}

  //! The method called after each minimization iteration.
  virtual PetscErrorCode 
  iteration( IPTaoTikhonovProblem<ForwardProblem> &problem,
             double eta, int iter,
             double objectiveValue, double designValue,
             DesignVec &d, DesignVec &diff_d, DesignVec &grad_d,
             StateVec &u,  StateVec &diff_u,  DesignVec &grad_u,
             DesignVec &gradient) = 0;
};


//! \brief Defines a Tikhonov minimization problem to be solved with a TaoBasicSolver.
/*! Suppose \f$F\f$ is a map from a space \f$D\f$ of design variables to a space \f$S\f$ of
state variables and we wish to solve a possibly ill-posed problem of the form
\f[ F(d) = u \f]
where \f$u\f$ is know and \f$d\f$ is unknown.  Approximate solutions can be obtained by 
finding minimizers of an associated Tikhonov functional
\f[
J(d) = J_{S}(F(d)-u) + \frac{1}{\eta}J_{D}(d-d_0)
\f]
where \$J_{D}\$ and \$J_{S}\$ are functionals on the spaces \f$D\f$ and \f$S\f$ respectively,
\f$\eta\f$ is a penalty parameter, and \f$d_0\f$ is a best a-priori guess for the the solution.
The IPTaoTikhonovProblem class encapuslates all of the data required to formulate the minimization
problem as a Problem tha can be solved using a TaoBasicSolver. It is templated on the
the class ForwardProblem which defines the class of the forward map \f$F\f$ as well as the
spaces \f$D\f$ and \f$S\f$. An instance of ForwardProblem, along with 
specific functionals \f$J_D\f$ and \f$J_S\f$, the parameter \f$\eta\f$, and the data 
\f$y\f$ and \f$x_0\f$ are provided on constructing a IPTaoTikhonovProblem.

For example, if the SSATaucForwardProblem class defines the map taking yield stresses \f$\tau_c\f$
to the corresponding surface velocity field solving the SSA, a schematic setup of solving
the associated Tikhonov problem goes as follows.

\code
SSATaucForwardProblem forwardProblem(grid);
L2NormFunctional2S designFunctional(grid); //J_X
L2NormFunctional2V stateFunctional(grid);  //J_Y
IceModelVec2V u_obs;     // Set this to the surface velocity observations.
IceModelVec2S tauc_0;    // Set this to the initial guess for tauc.
double eta;           // Set this to the desired penalty parameter.

typedef InvSSATauc IPTaoTikhonovProblem<SSATaucForwardProblem>;
InvSSATauc tikhonovProblem(forwardProblem,tauc_0,u_obs,eta,designFunctional,stateFunctional);

TaoBasicSolver<InvSSATauc> solver(com,"tao_cg",tikhonovProblem);

TerminationReason::Ptr reason;
solver.solve(reason);

if(reason->succeeded()) {
  printf("Success: %s\n",reason->description().c_str());
} else {
  printf("Failure: %s\n",reason->description().c_str());
}
\endcode

The class ForwardProblem that defines the forward problem
must have the following characteristics:

<ol>
<li> Contains typedefs for DesignVec and StateVec that effectively
   define the function spaces \f$D\f$ and \f$S\f$.  E.g. 

\code
   typedef IceModelVec2S DesignVec;
   typedef IceModelVec2V StateVec;
\endcode
   would be appropriate for a map from basal yeild stress to surface velocities.

<li> A method
\code
  PetscErrorCode linearize_at( DesignVec &d, TerminationReason::Ptr &reason);
\endcode
that instructs the class to compute the value of F and 
anything needed to compute its linearization at \a d.   This is the first method
called when working with a new iterate of \a d.

<li> A method
\code
  StateVec &solution()
\endcode
that returns the most recently computed value of \f$F(d)\f$ 
as computed by a call to linearize_at.

<li> A method
\code
PetscErrorCode apply_linearization_transpose(StateVec &du, DesignVec &dzeta);
\endcode
that computes the action of \f$(F')^t\f$,
where \f$F'\f$ is the linearization of \f$F\f$ at the current iterate, and the transpose
is computed in the standard sense (i.e. thinking of \f$F'\f$ as a matrix with respect
to the bases implied by the DesignVec and StateVec spaces).  The need for a transpose arises
because
\f[
\frac{d}{dt} J_{S}(F(d+t\delta d)-u) = [DJ_S]_{k}\; F'_{kj} \; \delta d
\f]
and hence the gradient of the term \f$J_{S}(F(d)-u)\f$ with respect to \f$d\f$ is given
by
\f[
(F')^t (\nabla J_S)^t.
\f]
</ol>

*/
template<class ForwardProblem> class IPTaoTikhonovProblem
{
public:

  typedef typename ForwardProblem::DesignVec DesignVec;
  typedef typename ForwardProblem::StateVec StateVec;

  /*! Constructs a Tikhonov problem:
  
           Minimize \f$J(d) = J_S(F(d)-u_obs) + \frac{1}{\eta} J_D(d-d0)  \f$

      that can be solved with a TaoBasicSolver.
      
      @param forward Class defining the map F.  See class-level documentation for requirements of F.
      @param      d0 Best a-priori guess for the design parameter.
      @param   u_obs State parameter to match (i.e. approximately solve F(d)=u_obs)
      @param     eta Penalty parameter/Lagrange multiplier.  Take eta to zero to impose more regularization to an ill posed problem.
      @param   designFunctional The functional \f$J_D\f$
      @param    stateFunctional The functional \f$J_S\f$
  */

  IPTaoTikhonovProblem( ForwardProblem &forward, DesignVec &d0, StateVec &u_obs, double eta, 
                  IPFunctional<DesignVec>&designFunctional, IPFunctional<StateVec>&stateFunctional);

  virtual ~IPTaoTikhonovProblem();


  //! Sets the initial guess for minimization iterations. If this isn't set explicitly,
  //  the parameter \f$d0\f$ appearing the in the Tikhonov functional will be used.
  virtual PetscErrorCode setInitialGuess( DesignVec &d) {
    PetscErrorCode ierr;
    ierr = m_dGlobal.copy_from(d); CHKERRQ(ierr);
    return 0;
  }

  //! Callback provided to TAO for objective evaluation.
  virtual PetscErrorCode evaluateObjectiveAndGradient(Tao tao, Vec x, double *value, Vec gradient);

  //! Add an object to the list of objects to be called after each iteration.
  virtual void addListener( typename IPTaoTikhonovProblemListener<ForwardProblem>::Ptr listener) {
    m_listeners.push_back(listener);
  }

  //! Final value of \f$F(d)\f$, where \f$d\f$ is the solution of the minimization.
  virtual StateVec &stateSolution() {
    return m_forward.solution();
  }

  //! Value of \f$d\f$, the solution of the minimization problem.
  virtual DesignVec &designSolution() {
    return m_d;
  }

  //! Callback from TaoBasicSolver, used to wire the connections between a Tao and
  //  the current class.
  virtual PetscErrorCode connect(Tao tao);

  //! Callback from TAO after each iteration.  The call is forwarded to each element of our list of listeners.
  virtual PetscErrorCode monitorTao(Tao tao);

  //! Callback from TAO to detect convergence.  Allows us to implement a custom convergence check.
  virtual PetscErrorCode convergenceTest(Tao tao);

  //! Callback from TaoBasicSolver to form the starting iterate for the minimization.  See also
  //  setInitialGuess.
  virtual PetscErrorCode formInitialGuess(Vec *v, TerminationReason::Ptr & reason) {
    *v = m_dGlobal.get_vec();
    reason = GenericTerminationReason::success();
    return 0;
  }

protected:

  virtual PetscErrorCode construct();
  
  IceGrid *m_grid;
  
  ForwardProblem &m_forward;

  DesignVec m_d;           ///< Current iterate of design parameter
  DesignVec m_dGlobal;     ///< Initial iterate of design parameter, stored without ghosts for the benefit of TAO.
  DesignVec &m_d0;         ///< A-priori estimate of design parameter
  DesignVec m_d_diff;      ///< Storage for (m_d-m_d0)

  StateVec &m_u_obs;       ///< State parameter to match via F(d)=u_obs 
  StateVec m_u_diff;       ///< Storage for F(d)-u_obs

  StateVec m_adjointRHS;   ///< Temporary storage used in gradient computation.

  DesignVec m_grad_design; ///< Gradient of \f$J_D\f$ at the current iterate.
  DesignVec m_grad_state;  ///< Gradient of \f$J_S\f$ at the current iterate.
  DesignVec m_grad;        /**< Weighted sum of the design and state gradients
                                corresponding to the gradient of the Tikhonov functional \f$J\f$. */

  double m_eta;         ///<  Penalty parameter/Lagrange multiplier.

  double m_val_design;  ///<  Value of \f$J_D\f$ at the current iterate.
  double m_val_state;   ///<  Value of \f$J_S\f$ at the current iterate.

  IPFunctional<IceModelVec2S> &m_designFunctional;  //<! Implementation of \f$J_D\f$.
  IPFunctional<IceModelVec2V> &m_stateFunctional;   //<! Implementation of \f$J_S\f$.

  std::vector<typename IPTaoTikhonovProblemListener<ForwardProblem>::Ptr> m_listeners; ///< List of iteration callbacks.

  double m_tikhonov_atol;  ///< Convergence parameter: convergence stops when \f$||J_D||_2 <\f$ m_tikhonov_rtol.
  double m_tikhonov_rtol;  /**< Convergence parameter: convergence stops when \f$||J_D||_2 \f$ is 
                                  less than m_tikhonov_rtol times the maximum of the gradient of \f$J_S\f$ and
                                  \f$(1/\eta)J_D\f$.  This occurs when the two terms forming the sum of the gradient
                                  of \f$J\f$ point in roughly opposite directions with the same magnitude. */

};

template<class ForwardProblem> IPTaoTikhonovProblem<ForwardProblem>::IPTaoTikhonovProblem( ForwardProblem &forward,
                 DesignVec &d0, StateVec &u_obs, double eta,
                 IPFunctional<DesignVec> &designFunctional, IPFunctional<StateVec> &stateFunctional ):
                  m_forward(forward), m_d0(d0), m_u_obs(u_obs), m_eta(eta),
                  m_designFunctional(designFunctional), m_stateFunctional(stateFunctional)
{
  PetscErrorCode ierr = this->construct();
  CHKERRCONTINUE(ierr);
  assert(ierr == 0);
}


template<class ForwardProblem> PetscErrorCode IPTaoTikhonovProblem<ForwardProblem>::construct() {
  PetscErrorCode ierr;

  m_grid = m_d0.get_grid();

  m_tikhonov_atol = m_grid->config.get("tikhonov_atol");
  m_tikhonov_rtol = m_grid->config.get("tikhonov_rtol");

  int design_stencil_width = m_d0.get_stencil_width();
  int state_stencil_width = m_u_obs.get_stencil_width();
  ierr = m_d.create(*m_grid, "design variable", WITH_GHOSTS, design_stencil_width); CHKERRQ(ierr);
  ierr = m_dGlobal.create(*m_grid, "design variable (global)", WITHOUT_GHOSTS, design_stencil_width); CHKERRQ(ierr);
  ierr = m_dGlobal.copy_from(m_d0); CHKERRQ(ierr);

  ierr = m_u_diff.create( *m_grid, "state residual", WITH_GHOSTS, state_stencil_width); CHKERRQ(ierr);
  ierr = m_d_diff.create( *m_grid, "design residual", WITH_GHOSTS, design_stencil_width); CHKERRQ(ierr);

  ierr = m_grad_state.create( *m_grid, "state gradient", WITHOUT_GHOSTS, design_stencil_width); CHKERRQ(ierr);
  ierr = m_grad_design.create( *m_grid, "design gradient", WITHOUT_GHOSTS, design_stencil_width); CHKERRQ(ierr);
  ierr = m_grad.create( *m_grid, "gradient", WITHOUT_GHOSTS, design_stencil_width); CHKERRQ(ierr);

  ierr = m_adjointRHS.create(*m_grid,"work vector", WITHOUT_GHOSTS, design_stencil_width); CHKERRQ(ierr);

  return 0;
}

template<class ForwardProblem> IPTaoTikhonovProblem<ForwardProblem>::~IPTaoTikhonovProblem() {}

template<class ForwardProblem> PetscErrorCode IPTaoTikhonovProblem<ForwardProblem>::connect(Tao tao) {
  PetscErrorCode ierr;
  typedef TaoObjGradCallback<IPTaoTikhonovProblem<ForwardProblem>,&IPTaoTikhonovProblem<ForwardProblem>::evaluateObjectiveAndGradient> ObjGradCallback; 
  ierr = ObjGradCallback::connect(tao,*this); CHKERRQ(ierr);
  ierr = TaoMonitorCallback< IPTaoTikhonovProblem<ForwardProblem> >::connect(tao,*this); CHKERRQ(ierr);
  ierr = TaoConvergenceCallback< IPTaoTikhonovProblem<ForwardProblem> >::connect(tao,*this); CHKERRQ(ierr);

  double fatol = 1e-10, frtol = 1e-20;
  double gatol = PETSC_DEFAULT, grtol = PETSC_DEFAULT, gttol = PETSC_DEFAULT;
  ierr = TaoSetTolerances(tao, fatol, frtol, gatol, grtol, gttol); CHKERRQ(ierr);

  return 0;
}

template<class ForwardProblem> PetscErrorCode IPTaoTikhonovProblem<ForwardProblem>::monitorTao(Tao tao) {
  PetscErrorCode ierr;
  
  int its;
  ierr =  TaoGetSolutionStatus(tao, &its, NULL, NULL, NULL, NULL, NULL ); CHKERRQ(ierr);
  
  int nListeners = m_listeners.size();
  for(int k=0; k<nListeners; k++) {
   ierr = m_listeners[k]->iteration(*this,m_eta,
                 its,m_val_design,m_val_state,
                 m_d, m_d_diff, m_grad_design,
                 m_forward.solution(), m_u_diff, m_grad_state,
                 m_grad );
   CHKERRQ(ierr);
  }
  return 0;
}

template<class ForwardProblem> PetscErrorCode IPTaoTikhonovProblem<ForwardProblem>::convergenceTest(Tao tao) {
  PetscErrorCode ierr;
  double designNorm, stateNorm, sumNorm;
  double dWeight, sWeight;
  dWeight = 1/m_eta;
  sWeight = 1;
  
  ierr = m_grad_design.norm(NORM_2,designNorm); CHKERRQ(ierr);
  ierr = m_grad_state.norm(NORM_2,stateNorm); CHKERRQ(ierr);
  ierr = m_grad.norm(NORM_2,sumNorm); CHKERRQ(ierr);
  designNorm *= dWeight;    
  stateNorm  *= sWeight;
  
#if PETSC_VERSION_LT(3,5,0)
  if (sumNorm < m_tikhonov_atol) {
    ierr = TaoSetTerminationReason(tao, TAO_CONVERGED_GATOL); CHKERRQ(ierr);    
  } else if( sumNorm < m_tikhonov_rtol*PetscMax(designNorm,stateNorm) ) {
    ierr = TaoSetTerminationReason(tao,TAO_CONVERGED_USER); CHKERRQ(ierr);
  } else {
    ierr = TaoDefaultConvergenceTest(tao,NULL); CHKERRQ(ierr);
  }
#else
  if (sumNorm < m_tikhonov_atol) {
    ierr = TaoSetConvergedReason(tao, TAO_CONVERGED_GATOL); CHKERRQ(ierr);
  } else if (sumNorm < m_tikhonov_rtol*PetscMax(designNorm,stateNorm)) {
    ierr = TaoSetConvergedReason(tao,TAO_CONVERGED_USER); CHKERRQ(ierr);
  } else {
    ierr = TaoDefaultConvergenceTest(tao,NULL); CHKERRQ(ierr);
  }
#endif
  return 0;
}

template<class ForwardProblem> PetscErrorCode IPTaoTikhonovProblem<ForwardProblem>::evaluateObjectiveAndGradient(Tao tao, Vec x, double *value, Vec gradient) {
  PetscErrorCode ierr;

  // Variable 'x' has no ghosts.  We need ghosts for computation with the design variable.
  ierr = m_d.copy_from_vec(x); CHKERRQ(ierr);

  TerminationReason::Ptr reason;
  ierr = m_forward.linearize_at(m_d, reason); CHKERRQ(ierr);
  if(reason->failed()) {
    ierr = verbPrintf(2,m_grid->com,"IPTaoTikhonovProblem::evaluateObjectiveAndGradient failure in forward solve\n%s\n",reason->description().c_str()); CHKERRQ(ierr);
#if PETSC_VERSION_LT(3,5,0)
    ierr = TaoSetTerminationReason(tao,TAO_DIVERGED_USER); CHKERRQ(ierr);
#else
    ierr = TaoSetConvergedReason(tao,TAO_DIVERGED_USER); CHKERRQ(ierr);
#endif
    return 0;
  }

  ierr = m_d_diff.copy_from(m_d); CHKERRQ(ierr);
  ierr = m_d_diff.add(-1,m_d0); CHKERRQ(ierr);
  ierr = m_designFunctional.gradientAt(m_d_diff,m_grad_design); CHKERRQ(ierr);

  ierr = m_u_diff.copy_from(m_forward.solution()); CHKERRQ(ierr);
  ierr = m_u_diff.add(-1, m_u_obs); CHKERRQ(ierr);

  // The following computes the reduced gradient.
  ierr = m_stateFunctional.gradientAt(m_u_diff,m_adjointRHS); CHKERRQ(ierr);  
  ierr = m_forward.apply_linearization_transpose(m_adjointRHS,m_grad_state); CHKERRQ(ierr);

  ierr = m_grad.copy_from(m_grad_design); CHKERRQ(ierr);
  ierr = m_grad.scale(1./m_eta); CHKERRQ(ierr);    
  ierr = m_grad.add(1,m_grad_state); CHKERRQ(ierr);

  ierr = VecCopy(m_grad.get_vec(), gradient); CHKERRQ(ierr);

  double valDesign, valState;
  ierr = m_designFunctional.valueAt(m_d_diff,&valDesign); CHKERRQ(ierr);
  ierr = m_stateFunctional.valueAt(m_u_diff,&valState); CHKERRQ(ierr);

  m_val_design = valDesign;
  m_val_state = valState;
  
  *value = valDesign / m_eta + valState;

  return 0;
}



#endif /* end of include guard: IPTAOTIKHONOVPROBLEM_HH_4NMM724B */
