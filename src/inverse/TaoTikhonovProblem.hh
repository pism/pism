// Copyright (C) 2012,2013  David Maxwell
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


#ifndef TAOTIKHONOVPROBLEM_HH_4NMM724B
#define TAOTIKHONOVPROBLEM_HH_4NMM724B

#include <tr1/memory>

#include "TaoUtil.hh"
#include "Functional.hh"
#include <assert.h>

template<class ForwardProblem> class TaoTikhonovProblem;

/*! A class for objects receiving iteration callbacks from a TaoTikhonovProblem.  These 
    callbacks can be used to monitor the solution, plot iterations, print diagnostic messages, etc. 
    TaoTikhonovProblemListeners are ususally used via a reference counted pointer 
    TaoTikhonovProblemListener::Ptr to allow for good memory management when Listeners are 
    created as subclasses of python classes. It would have been better to nest this inside of 
    TaoTikhonovProblem, but SWIG has a hard time with nested classes, so it's outer instead.*/
template<class ForwardProblem> class TaoTikhonovProblemListener {
public:
  typedef std::tr1::shared_ptr<TaoTikhonovProblemListener> Ptr;

  typedef typename ForwardProblem::DesignVec DesignVec;
  typedef typename ForwardProblem::StateVec StateVec;

  TaoTikhonovProblemListener() {}
  virtual ~TaoTikhonovProblemListener() {}

  //! The method called after each minimization iteration.
  virtual PetscErrorCode 
  iteration( TaoTikhonovProblem<ForwardProblem> &problem,
             PetscReal eta, PetscInt iter,
             PetscReal objectiveValue, PetscReal designValue,
             DesignVec &d, DesignVec &diff_d, DesignVec &grad_d,
             StateVec &u,  StateVec &diff_u,  DesignVec &grad_u,
             DesignVec &gradient) = 0;
};


//! \brief Defines a Tikhonov minimization problem to be solved with a TaoBasicSolver.
/*! Suppose \f$F\f$ is a map from a space \f$D\f$ of design variables to a space \f$S\f$ of
state variables and we wish to solve a possibly ill-posed problem of the form
\f[ F(d) = u \f]
where \f$u\f$ is know and \f$d\f$ is unknown.  Approximate solutions can be obtained by 
finding minimizers of an associated Tikhonov problem
\[
J(d) = J_{S}(F(d)-u) + \frac{1}{\eta}J_{D}(d-d_0)
\]
where \$J_{D}\$ and \$J_{S}\$ are functionals on the spaces \f$D\f$ and \f$S\f$ respectively,
$\eta$ is a penalty paramter, and \f$d_0\f$ is a best initial guess for the the solution.

The TaoTikhonovProblemTAO class encapuslates all of the data required to formulate the minimization
problem as a Problem tha can be solved using a TaoBasicSolver. It is templated on the
the class ForwardProblem which defines the class of the forward map \f$F\f$ as well as the
spaces \f$D\f$ and \f$S\f$. An instance of ForwardProblem, along with 
specific functionals \f$J_D\f$ and f$J_S\f$, the parameter \f$\eta\f$, and the data 
\f$y\f$ and $x_0$ are provided on constructing a TaoTikhonovProblem.

For example, if the SSATaucForwardProblem class defines the map taking yield stresses \f$\tau_c\f$
to the corresponding surface velocity field solving the SSA, a schematic setup of solving
the associated Tikhonov problem goes as follows.

\code
SSATaucForwardProblem forwardProblem(grid);
L2NormFunctional2S designFunctional(grid); //J_X
L2NormFunctional2V stateFunctional(grid);  //J_Y
IceModelVec2V u_obs;     // Set this to the surface velocity observations.
IceModelVec2S tauc_0;    // Set this to the initial guess for tauc.
PetscReal eta;           // Set this to the desired penalty parameter.

typedef InvSSATauc TaoTikhonovProblem<SSATaucForwardProblem>;
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

*/
template<class ForwardProblem> class TaoTikhonovProblem
{
public:

  typedef typename ForwardProblem::DesignVec DesignVec;
  typedef typename ForwardProblem::StateVec StateVec;

  /*! Constructs a Tikhonov problem:
  Minimize \f$J(d) = J_Y(F(d)-u) + J_X(d-d0)
  
  */

  TaoTikhonovProblem( ForwardProblem &forward, DesignVec &d0, StateVec &u_obs, PetscReal eta, 
                  Functional<DesignVec>&designFunctional, Functional<StateVec>&stateFunctional);

  virtual ~TaoTikhonovProblem();

  virtual PetscErrorCode setInitialGuess( DesignVec &d) {
    PetscErrorCode ierr;
    ierr = m_dGlobal.copy_from(d); CHKERRQ(ierr);
    return 0;
  }

  virtual PetscErrorCode evaluateObjectiveAndGradient(TaoSolver tao, Vec x, PetscReal *value, Vec gradient);

  virtual void addListener( typename TaoTikhonovProblemListener<ForwardProblem>::Ptr listener) {
    m_listeners.push_back(listener);
  }

  virtual StateVec &stateSolution() {
    return m_forward.solution();
  }

  virtual DesignVec &designSolution() {
    return m_d;
  }

  virtual PetscErrorCode connect(TaoSolver tao);

  virtual PetscErrorCode monitorTao(TaoSolver tao);

  virtual PetscErrorCode convergenceTest(TaoSolver tao); 

  virtual PetscErrorCode formInitialGuess(Vec *v, TerminationReason::Ptr & reason) {
    *v = m_dGlobal.get_vec();
    reason = GenericTerminationReason::success();
    return 0;
  }

protected:

  virtual PetscErrorCode construct();
  
  IceGrid *m_grid;
  
  ForwardProblem &m_forward;

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

  std::vector<typename TaoTikhonovProblemListener<ForwardProblem>::Ptr> m_listeners;

  PetscReal m_tikhonov_atol;
  PetscReal m_tikhonov_rtol;

};

template<class ForwardProblem> TaoTikhonovProblem<ForwardProblem>::TaoTikhonovProblem( ForwardProblem &forward,
                 DesignVec &d0, StateVec &u_obs, PetscReal eta,
                 Functional<DesignVec> &designFunctional, Functional<StateVec> &stateFunctional ):
                  m_forward(forward), m_d0(d0), m_u_obs(u_obs), m_eta(eta),
                  m_designFunctional(designFunctional), m_stateFunctional(stateFunctional)
{
  PetscErrorCode ierr = this->construct();
  CHKERRCONTINUE(ierr);
  assert(ierr == 0);
}


template<class ForwardProblem> PetscErrorCode TaoTikhonovProblem<ForwardProblem>::construct() {
  PetscErrorCode ierr;

  m_grid = m_d0.get_grid();

  m_tikhonov_atol = m_grid->config.get("tikhonov_atol");
  m_tikhonov_rtol = m_grid->config.get("tikhonov_rtol");

  PetscInt design_stencil_width = m_d0.get_stencil_width();
  PetscInt state_stencil_width = m_u_obs.get_stencil_width();
  ierr = m_d.create(*m_grid, "design variable", kHasGhosts, design_stencil_width); CHKERRQ(ierr);
  ierr = m_dGlobal.create(*m_grid, "design variable (global)", kNoGhosts, design_stencil_width); CHKERRQ(ierr);
  ierr = m_dGlobal.copy_from(m_d0); CHKERRQ(ierr);

  ierr = m_u_diff.create( *m_grid, "state residual", kHasGhosts, state_stencil_width); CHKERRQ(ierr);
  ierr = m_d_diff.create( *m_grid, "design residual", kHasGhosts, design_stencil_width); CHKERRQ(ierr);

  ierr = m_grad_state.create( *m_grid, "state gradient", kNoGhosts, design_stencil_width); CHKERRQ(ierr);
  ierr = m_grad_design.create( *m_grid, "design gradient", kNoGhosts, design_stencil_width); CHKERRQ(ierr);
  ierr = m_grad.create( *m_grid, "gradient", kNoGhosts, design_stencil_width); CHKERRQ(ierr);

  ierr = m_adjointRHS.create(*m_grid,"work vector", kNoGhosts, design_stencil_width); CHKERRQ(ierr);

  return 0;
}

template<class ForwardProblem> TaoTikhonovProblem<ForwardProblem>::~TaoTikhonovProblem() {}

template<class ForwardProblem> PetscErrorCode TaoTikhonovProblem<ForwardProblem>::connect(TaoSolver tao) {
  PetscErrorCode ierr;
  typedef TaoObjGradCallback<TaoTikhonovProblem<ForwardProblem>,&TaoTikhonovProblem<ForwardProblem>::evaluateObjectiveAndGradient> ObjGradCallback; 
  ierr = ObjGradCallback::connect(tao,*this); CHKERRQ(ierr);
  ierr = TaoMonitorCallback< TaoTikhonovProblem<ForwardProblem> >::connect(tao,*this); CHKERRQ(ierr);
  ierr = TaoConvergenceCallback< TaoTikhonovProblem<ForwardProblem> >::connect(tao,*this); CHKERRQ(ierr);

  PetscReal fatol = 1e-10, frtol = 1e-20;
  PetscReal gatol = PETSC_DEFAULT, grtol = PETSC_DEFAULT, gttol = PETSC_DEFAULT;
  ierr = TaoSetTolerances(tao, fatol, frtol, gatol, grtol, gttol); CHKERRQ(ierr);

  return 0;
}

template<class ForwardProblem> PetscErrorCode TaoTikhonovProblem<ForwardProblem>::monitorTao(TaoSolver tao) {
  PetscErrorCode ierr;
  
  PetscInt its;
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

template<class ForwardProblem> PetscErrorCode TaoTikhonovProblem<ForwardProblem>::convergenceTest(TaoSolver tao) {
  PetscErrorCode ierr;
  PetscReal designNorm, stateNorm, sumNorm;
  PetscReal dWeight, sWeight;
  dWeight = 1/m_eta;
  sWeight = 1;
  
  ierr = m_grad_design.norm(NORM_2,designNorm); CHKERRQ(ierr);
  ierr = m_grad_state.norm(NORM_2,stateNorm); CHKERRQ(ierr);
  ierr = m_grad.norm(NORM_2,sumNorm); CHKERRQ(ierr);
  designNorm *= dWeight;    
  stateNorm  *= sWeight;
  
  if( sumNorm < m_tikhonov_atol) {
    ierr = TaoSetTerminationReason(tao, TAO_CONVERGED_GATOL); CHKERRQ(ierr);    
  } else if( sumNorm < m_tikhonov_rtol*PetscMax(designNorm,stateNorm) ) {
    ierr = TaoSetTerminationReason(tao,TAO_CONVERGED_USER); CHKERRQ(ierr);
  } else {
    ierr = TaoDefaultConvergenceTest(tao,NULL); CHKERRQ(ierr);
  }
  return 0;
}

template<class ForwardProblem> PetscErrorCode TaoTikhonovProblem<ForwardProblem>::evaluateObjectiveAndGradient(TaoSolver tao, Vec x, PetscReal *value, Vec gradient) {
  PetscErrorCode ierr;

  // Variable 'x' has no ghosts.  We need ghosts for computation with the design variable.
  ierr = m_d.copy_from(x); CHKERRQ(ierr);

  TerminationReason::Ptr reason;
  ierr = m_forward.linearize_at(m_d, reason); CHKERRQ(ierr);
  if(reason->failed()) {
    ierr = verbPrintf(2,m_grid->com,"TaoTikhonovProblem::evaluateObjectiveAndGradient failure in forward solve\n%s\n",reason->description().c_str()); CHKERRQ(ierr);
    ierr = TaoSetTerminationReason(tao,TAO_DIVERGED_USER); CHKERRQ(ierr);
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

  ierr = m_grad.copy_to(gradient); CHKERRQ(ierr);      

  PetscReal valDesign, valState;
  ierr = m_designFunctional.valueAt(m_d_diff,&valDesign); CHKERRQ(ierr);
  ierr = m_stateFunctional.valueAt(m_u_diff,&valState); CHKERRQ(ierr);

  m_val_design = valDesign;
  m_val_state = valState;
  
  *value = valDesign / m_eta + valState;

  return 0;
}



#endif /* end of include guard: TAOTIKHONOVPROBLEM_HH_4NMM724B */
