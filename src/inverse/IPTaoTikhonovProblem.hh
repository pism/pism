// Copyright (C) 2012,2013,2014,2015,2016,2017,2020,2022,2023  David Maxwell and Constantine Khroulev
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

#ifndef PISM_IPTAOTIKHONOVPROBLEM_HH
#define PISM_IPTAOTIKHONOVPROBLEM_HH

#include <memory>

#include "pism/inverse/TaoUtil.hh"
#include "pism/inverse/functional/IPFunctional.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Logger.hh"
#include "pism/util/array/Array.hh" //
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"

namespace pism {

class Grid;

namespace inverse {

template <class ForwardProblem>
class IPTaoTikhonovProblem;

//! Iteration callback class for IPTaoTikhonovProblem
/** A class for objects receiving iteration callbacks from a
 * IPTaoTikhonovProblem. These callbacks can be used to monitor the
 * solution, plot iterations, print diagnostic messages, etc.
 * IPTaoTikhonovProblemListeners are ususally used via a reference
 * counted pointer IPTaoTikhonovProblemListener::Ptr to allow for good
 * memory management when Listeners are created as subclasses of
 * python classes. It would have been better to nest this inside of
 * IPTaoTikhonovProblem, but SWIG has a hard time with nested classes,
 * so it's outer instead.
 */
template <class ForwardProblem>
class IPTaoTikhonovProblemListener {
public:
  typedef std::shared_ptr<IPTaoTikhonovProblemListener> Ptr;

  typedef std::shared_ptr<typename ForwardProblem::DesignVec> DesignVecPtr;
  typedef std::shared_ptr<typename ForwardProblem::StateVec> StateVecPtr;
  typedef std::shared_ptr<typename ForwardProblem::StateVec1> StateVec1Ptr;

  IPTaoTikhonovProblemListener() {
  }
  virtual ~IPTaoTikhonovProblemListener() {
  }

  //! The method called after each minimization iteration.
  virtual void iteration(IPTaoTikhonovProblem<ForwardProblem> &problem, double eta, int iter,
                         double objectiveValue, double designValue, const DesignVecPtr d,
                         const DesignVecPtr diff_d, const DesignVecPtr grad_d, const StateVecPtr u,
                         const StateVecPtr diff_u, const DesignVecPtr grad_u,
                         const DesignVecPtr gradient) = 0;
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
  array::Vector u_obs;     // Set this to the surface velocity observations.
  array::Scalar tauc_0;    // Set this to the initial guess for tauc.
  double eta;           // Set this to the desired penalty parameter.

  typedef InvSSATauc IPTaoTikhonovProblem<SSATaucForwardProblem>;
  InvSSATauc tikhonovProblem(forwardProblem,tauc_0,u_obs,eta,designFunctional,stateFunctional);

  TaoBasicSolver<InvSSATauc> solver(com, "cg", tikhonovProblem);

  TerminationReason::Ptr reason = solver.solve();

  if (reason->succeeded()) {
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
  typedef array::Scalar DesignVec;
  typedef array::Vector StateVec;
  \endcode
  would be appropriate for a map from basal yeild stress to surface velocities.

  <li> A method
  \code
  TerminationReason::Ptr linearize_at(DesignVec &d);
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
  void apply_linearization_transpose(StateVec &du, DesignVec &dzeta);
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
template <class ForwardProblem>
class IPTaoTikhonovProblem {
public:
  typedef typename ForwardProblem::DesignVec DesignVec;
  typedef typename ForwardProblem::StateVec StateVec;
  typedef typename ForwardProblem::StateVec1 StateVec1;

  typedef typename ForwardProblem::DesignVecGhosted DesignVecGhosted;
  typedef std::shared_ptr<typename ForwardProblem::DesignVecGhosted> DesignVecGhostedPtr;

  typedef std::shared_ptr<typename ForwardProblem::DesignVec> DesignVecPtr;
  typedef std::shared_ptr<typename ForwardProblem::StateVec> StateVecPtr;
  typedef std::shared_ptr<typename ForwardProblem::StateVec1> StateVec1Ptr;

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

  IPTaoTikhonovProblem(ForwardProblem &forward, DesignVec &d0, StateVec &u_obs, double eta,
                       IPFunctional<DesignVec> &designFunctional,
                       IPFunctional<StateVec> &stateFunctional);

  virtual ~IPTaoTikhonovProblem();


  //! Sets the initial guess for minimization iterations. If this isn't set explicitly,
  //  the parameter \f$d0\f$ appearing the in the Tikhonov functional will be used.
  virtual void setInitialGuess(DesignVec &d) {
    m_dGlobal.copy_from(d);
  }

  //! Callback provided to TAO for objective evaluation.
  virtual void evaluateObjectiveAndGradient(Tao tao, Vec x, double *value, Vec gradient);

  //! Add an object to the list of objects to be called after each iteration.
  virtual void addListener(typename IPTaoTikhonovProblemListener<ForwardProblem>::Ptr listener) {
    m_listeners.push_back(listener);
  }

  //! Final value of \f$F(d)\f$, where \f$d\f$ is the solution of the minimization.
  virtual StateVecPtr stateSolution() {
    return m_forward.solution();
  }

  //! Value of \f$d\f$, the solution of the minimization problem.
  virtual DesignVecPtr designSolution() {
    return m_d;
  }

  //! Callback from TaoBasicSolver, used to wire the connections between a Tao and
  //  the current class.
  virtual void connect(Tao tao);

  //! Callback from TAO after each iteration.  The call is forwarded to each element of our list of listeners.
  virtual void monitorTao(Tao tao);

  //! Callback from TAO to detect convergence.  Allows us to implement a custom convergence check.
  virtual void convergenceTest(Tao tao);

  //! Callback from TaoBasicSolver to form the starting iterate for the minimization.  See also
  //  setInitialGuess.
  virtual std::shared_ptr<TerminationReason> formInitialGuess(Vec *v) {
    *v = m_dGlobal.vec();
    return GenericTerminationReason::success();
  }

protected:
  std::shared_ptr<const Grid> m_grid;

  ForwardProblem &m_forward;

  /// Current iterate of design parameter
  DesignVecGhostedPtr m_d; // ghosted
  /// Initial iterate of design parameter, stored without ghosts for the benefit of TAO.
  DesignVec m_dGlobal;
  /// A-priori estimate of design parameter
  DesignVec &m_d0;
  /// Storage for (m_d-m_d0)
  DesignVecPtr m_d_diff; // ghosted

  /// State parameter to match via F(d)=u_obs
  StateVec &m_u_obs;
  /// Storage for F(d)-u_obs
  StateVec1Ptr m_u_diff; // ghosted

  /// Temporary storage used in gradient computation.
  StateVec m_adjointRHS;

  /// Gradient of \f$J_D\f$ at the current iterate.
  DesignVecPtr m_grad_design;
  /// Gradient of \f$J_S\f$ at the current iterate.
  DesignVecPtr m_grad_state;
  /** Weighted sum of the design and state gradients corresponding to
   * the gradient of the Tikhonov functional \f$J\f$.
   */
  DesignVecPtr m_grad;

  ///  Penalty parameter/Lagrange multiplier.
  double m_eta;

  ///  Value of \f$J_D\f$ at the current iterate.
  double m_val_design;
  ///  Value of \f$J_S\f$ at the current iterate.
  double m_val_state;

  /// Implementation of \f$J_D\f$.
  IPFunctional<array::Scalar> &m_designFunctional;
  /// Implementation of \f$J_S\f$.
  IPFunctional<array::Vector> &m_stateFunctional;

  /// List of iteration callbacks.
  std::vector<typename IPTaoTikhonovProblemListener<ForwardProblem>::Ptr> m_listeners;

  /// Convergence parameter: convergence stops when \f$||J_D||_2 <\f$ m_tikhonov_rtol.
  double m_tikhonov_atol;

  /** Convergence parameter: convergence stops when \f$||J_D||_2 \f$
   * is less than m_tikhonov_rtol times the maximum of the gradient of
   * \f$J_S\f$ and \f$(1/\eta)J_D\f$. This occurs when the two terms
   * forming the sum of the gradient of \f$J\f$ point in roughly
   * opposite directions with the same magnitude.
  */
  double m_tikhonov_rtol;
};

template <class ForwardProblem>
IPTaoTikhonovProblem<ForwardProblem>::IPTaoTikhonovProblem(
    ForwardProblem &forward, DesignVec &d0, StateVec &u_obs, double eta,
    IPFunctional<DesignVec> &designFunctional, IPFunctional<StateVec> &stateFunctional)
    : m_grid(d0.grid()),
      m_forward(forward),
      m_dGlobal(d0.grid(), "design variable (global)"),
      m_d0(d0),
      m_u_obs(u_obs),
      m_adjointRHS(d0.grid(), "work vector"),
      m_eta(eta),
      m_designFunctional(designFunctional),
      m_stateFunctional(stateFunctional) {

  m_tikhonov_atol = m_grid->ctx()->config()->get_number("inverse.tikhonov.atol");
  m_tikhonov_rtol = m_grid->ctx()->config()->get_number("inverse.tikhonov.rtol");

  m_d = std::make_shared<DesignVecGhosted>(m_grid, "design variable");

  m_dGlobal.copy_from(m_d0);

  m_u_diff = std::make_shared<StateVec1>(m_grid, "state residual");

  m_d_diff = std::make_shared<DesignVecGhosted>(m_grid, "design residual");

  m_grad_state = std::make_shared<DesignVec>(m_grid, "state gradient");

  m_grad_design = std::make_shared<DesignVec>(m_grid, "design gradient");

  m_grad = std::make_shared<DesignVec>(m_grid, "gradient");
}

template<class ForwardProblem>
IPTaoTikhonovProblem<ForwardProblem>::~IPTaoTikhonovProblem() {
  // empty
}

template<class ForwardProblem>
void IPTaoTikhonovProblem<ForwardProblem>::connect(Tao tao) {
  typedef taoutil::TaoObjGradCallback<IPTaoTikhonovProblem<ForwardProblem>,
                             &IPTaoTikhonovProblem<ForwardProblem>::evaluateObjectiveAndGradient> ObjGradCallback; 

  ObjGradCallback::connect(tao,*this);

  taoutil::TaoMonitorCallback< IPTaoTikhonovProblem<ForwardProblem> >::connect(tao,*this);

  taoutil::TaoConvergenceCallback< IPTaoTikhonovProblem<ForwardProblem> >::connect(tao,*this);

  double
    gatol = PETSC_DEFAULT,
    grtol = PETSC_DEFAULT,
    gttol = PETSC_DEFAULT;

#if PETSC_VERSION_LT(3,7,0)
  double fatol = 1e-10, frtol = 1e-20;
  PetscErrorCode ierr = TaoSetTolerances(tao, fatol, frtol, gatol, grtol, gttol);
  PISM_CHK(ierr, "TaoSetTolerances");
#else
  PetscErrorCode ierr = TaoSetTolerances(tao, gatol, grtol, gttol);
  PISM_CHK(ierr, "TaoSetTolerances");
#endif
}

template<class ForwardProblem>
void IPTaoTikhonovProblem<ForwardProblem>::monitorTao(Tao tao) {
  PetscInt its;
  TaoGetSolutionStatus(tao, &its, NULL, NULL, NULL, NULL, NULL);

  int nListeners = m_listeners.size();
  for (int k=0; k<nListeners; k++) {
    m_listeners[k]->iteration(*this, m_eta,
                              its, m_val_design, m_val_state,
                              m_d, m_d_diff, m_grad_design,
                              m_forward.solution(), m_u_diff, m_grad_state,
                              m_grad);
  }
}

template<class ForwardProblem> void IPTaoTikhonovProblem<ForwardProblem>::convergenceTest(Tao tao) {
  double designNorm, stateNorm, sumNorm;
  double dWeight, sWeight;
  dWeight = 1/m_eta;
  sWeight = 1;
  
  designNorm = m_grad_design->norm(NORM_2)[0];
  stateNorm  = m_grad_state->norm(NORM_2)[0];
  sumNorm    = m_grad->norm(NORM_2)[0];

  designNorm *= dWeight;    
  stateNorm  *= sWeight;
  
  if (sumNorm < m_tikhonov_atol) {
    TaoSetConvergedReason(tao, TAO_CONVERGED_GATOL);
  } else if (sumNorm < m_tikhonov_rtol*std::max(designNorm,stateNorm)) {
    TaoSetConvergedReason(tao, TAO_CONVERGED_USER);
  } else {
    TaoDefaultConvergenceTest(tao, NULL);
  }
}

template<class ForwardProblem>
void IPTaoTikhonovProblem<ForwardProblem>::evaluateObjectiveAndGradient(Tao tao, Vec x,
                                                                        double *value, Vec gradient) {
  PetscErrorCode ierr;
  // Variable 'x' has no ghosts.  We need ghosts for computation with the design variable.
  {
    ierr = DMGlobalToLocalBegin(*m_d->dm(), x, INSERT_VALUES, m_d->vec());
    PISM_CHK(ierr, "DMGlobalToLocalBegin");

    ierr = DMGlobalToLocalEnd(*m_d->dm(), x, INSERT_VALUES, m_d->vec());
    PISM_CHK(ierr, "DMGlobalToLocalEnd");
  }

  auto reason = m_forward.linearize_at(*m_d);
  if (reason->failed()) {
    Logger::ConstPtr log = m_grid->ctx()->log();
    log->message(2,
                 "IPTaoTikhonovProblem::evaluateObjectiveAndGradient"
                 " failure in forward solve\n%s\n", reason->description().c_str());
    ierr = TaoSetConvergedReason(tao, TAO_DIVERGED_USER);
    PISM_CHK(ierr, "TaoSetConvergedReason");
    return;
  }

  m_d_diff->copy_from(*m_d);
  m_d_diff->add(-1, m_d0);
  m_designFunctional.gradientAt(*m_d_diff, *m_grad_design);

  m_u_diff->copy_from(*m_forward.solution());
  m_u_diff->add(-1, m_u_obs);

  // The following computes the reduced gradient.
  m_stateFunctional.gradientAt(*m_u_diff, m_adjointRHS);
  m_forward.apply_linearization_transpose(m_adjointRHS, *m_grad_state);

  m_grad->copy_from(*m_grad_design);
  m_grad->scale(1.0 / m_eta);
  m_grad->add(1, *m_grad_state);

  ierr = VecCopy(m_grad->vec(), gradient);
  PISM_CHK(ierr, "VecCopy");

  double valDesign, valState;
  m_designFunctional.valueAt(*m_d_diff, &valDesign);
  m_stateFunctional.valueAt(*m_u_diff, &valState);

  m_val_design = valDesign;
  m_val_state = valState;

  *value = valDesign / m_eta + valState;
}

} // end of namespace inverse
} // end of namespace pism

#endif // PISM_IPTAOTIKHONOVPROBLEM_HH
