// Copyright (C) 2012, 2013, 2014, 2015, 2017  David Maxwell and Constantine Khroulev
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

#ifndef TAOUTIL_HH_W42GJNRO
#define TAOUTIL_HH_W42GJNRO

#include <petsc.h>
#include <stdexcept>
#include <string>

#include "pism/util/pism_utilities.hh"
#include "pism/util/TerminationReason.hh"
#include "pism/util/petscwrappers/Tao.hh"
#include "pism/util/error_handling.hh"

namespace pism {

//! TAO (inverse modeling) utilities.
namespace taoutil {

//! Encapsulate TAO's TaoSolverTerminationReason codes as a PISM TerminationReason.
class TAOTerminationReason: public TerminationReason {
public:
  TAOTerminationReason(TaoConvergedReason r);
  virtual void get_description(std::ostream &desc,int indent_level=0);
};

//! \brief An interface for solving an optimization problem with TAO where the
//! problem itself is defined by a separate Problem class.
/*!  The primary interface to a TAO optimization problem is mediated by a PETSc-style
  `TaoSolver` object. The PISM TaoBasicSolver C++ class wraps a `TaoSolver` and some of 
  its initialization boilierplate, and allows a separate class to define the function to be minimized.

  To use a TaoBasicSolver you create a `Problem` class that defines the objective function and initial
  guess, as well any auxilliary callbacks desired.  The Problem class must define a

  \code
  void Problem::connect(TaoSolver solver);
  \endcode

  method which gives the `Problem` an opportunity to register its methods as callbacks to the solver, 
  perhaps taking advantage of the various `TaoFooCallback` classes provided in TaoUtil.hh to facilitate this.
  For example, a problem class MyProblem that did nothing more than register a combined objective/gradient
  callback could define

  \code
  void MyProblem::connect(TaoSolver tao) {
  typedef TaoObjGradCallback<Problem,&MyProblem::evaluateObjectiveAndGradient> ObjGradCallback; 
  ObjGradCallback::connect(tao,*this);
  }
  \endcode

  In addition to the `connect` method, a `Problem` must define 
  \code
  TerminationReason::Ptr MyProblem::formInitialGuess(Vec *v)
  \endcode
  which allows the problem to set the initial guess for optimization. If the minimization
  is successful, the solution will be found in the same vector that was returned by this method.

  Assuming a `MyProblem` called `problem` has been constructed, solution
  of the minimization is done using, for example, the TAO algorithm
  tao_cg:

  \code
  TaoBasicSolver<MyProblem> solver(com,"tao_cg",problem);

  TerminationReason::Ptr reason = solver.solve();

  if (reason->succeeded()) {
  printf("Success: %s\n",reason->description().c_str());
  } else {
  printf("Failure: %s\n",reason->description().c_str());
  }
  \endcode
*/
template<class Problem>
class TaoBasicSolver {
public:
    
  //! Construct a solver to solve `prob` using TAO algorithm `tao_type`.
  TaoBasicSolver(MPI_Comm comm, const std::string & tao_type, Problem &prob)
    : m_comm(comm), m_problem(prob) {
    PetscErrorCode ierr;

    ierr = TaoCreate(m_comm, m_tao.rawptr());
    PISM_CHK(ierr, "TaoCreate");

    ierr = TaoSetType(m_tao,tao_type.c_str());
    PISM_CHK(ierr, "TaoSetType");

    m_problem.connect(m_tao);

    ierr = TaoSetFromOptions(m_tao);
    PISM_CHK(ierr, "TaoSetFromOptions");
  }
  
  virtual ~TaoBasicSolver() {
    // empty
  }

  //! Solve the minimization problem.
  virtual TerminationReason::Ptr solve() {
    PetscErrorCode ierr;

    /* Solve the application */ 
    Vec x0;
    TerminationReason::Ptr reason = m_problem.formInitialGuess(&x0);
    if (reason->failed()) {
      TerminationReason::Ptr root_cause = reason;
      reason.reset(new GenericTerminationReason(-1, "Unable to form initial guess"));
      reason->set_root_cause(root_cause);
      return reason;
    }
    ierr = TaoSetInitialVector(m_tao, x0);
    PISM_CHK(ierr, "TaoSetInitialVector");

    ierr = TaoSolve(m_tao);
    PISM_CHK(ierr, "TaoSolve");

    TaoConvergedReason tao_reason;
    ierr = TaoGetConvergedReason(m_tao, &tao_reason);
    PISM_CHK(ierr, "TaoGetConvergedReason");

    reason.reset(new TAOTerminationReason(tao_reason));
    return reason;
  }

  virtual void setMaximumIterations(int max_it) {
    PetscErrorCode ierr = TaoSetMaximumIterations(m_tao, max_it);
    PISM_CHK(ierr, "TaoSetMaximumIterations");
  }
  
  virtual Problem &problem() {
    return m_problem;
  }

protected:
  MPI_Comm m_comm;
  petsc::Tao m_tao;
  Problem &m_problem;
};


//! \brief Adaptor to connect a TAO objective function callback to a C++ object method.
/*! The TAO library interfaces with user code via C-style callback functions.
  This class makes it convenient to associate a TAO Objective callback
  with a C++ object method. To assign 
  \code
  void MyObject::evaluateObjective(TaoSolver tao,Vec x, double *value);
  \endcode

  as the objective function to a `TaoSolver` `tao`, 

  \code
  MyObject obj;
  TaoObjectiveCallback<MyObject>::connect(tao,obj);
  \endcode

  The method name `evaluateObjective` for the callback is hard-coded.
  See TaoObjGradCallback for a technique to allow 
  the method name to be specified (at the expense of a little more cumbersome code).
*/
template<class Problem>
class TaoObjectiveCallback {
public:

  static void connect(Tao tao, Problem &p) {
    PetscErrorCode ierr;
    ierr = TaoSetObjectiveRoutine(tao,
                                  TaoObjectiveCallback<Problem>::callback,
                                  &p);
    PISM_CHK(ierr, "TaoSetObjectiveRoutine");
  }

protected:
  static PetscErrorCode callback(Tao tao, Vec x, double *value, void *ctx) {
    try {
      Problem *p = reinterpret_cast<Problem *>(ctx);
      PetscErrorCode ierr = p->evaluateObjective(tao,x,value); CHKERRQ(ierr);
    } catch (...) {
      MPI_Comm com = MPI_COMM_SELF;
      PetscErrorCode ierr = PetscObjectGetComm((PetscObject)tao, &com); CHKERRQ(ierr);
      handle_fatal_errors(com);
      SETERRQ(com, 1, "A PISM callback failed");
    }
    return 0;
  }
};


//! \brief Adaptor to connect a TAO monitoring callback to a C++ object method.
/*! The TAO library interfaces with user code via C-style callback functions.
  This class makes it convenient to associate a TAO Monitor callback
  with a C++ object method. To assign 
  \code
  void MyObject::monitorTao(Tao tao)
  \endcode

  as the objective function to a `Tao` `tao`, 

  \code
  MyObject obj;
  TaoMonitorCallback<MyObject>::connect(tao,obj);
  \endcode

  The method name `monitorTao` for the callback is hard-coded.
  See TaoObjGradCallback for a technique to allow 
  the method name to be specified (at the expense of a little more cumbersome code).
*/
template<class Problem>
class TaoMonitorCallback {
public:
  static void connect(Tao tao, Problem &p) {
    PetscErrorCode ierr = TaoSetMonitor(tao,
                                        TaoMonitorCallback<Problem>::callback,
                                        &p, NULL);
    PISM_CHK(ierr, "TaoSetMonitor");
  }
protected:
  static PetscErrorCode callback(Tao tao, void *ctx) {
    try {
      Problem *p = reinterpret_cast<Problem *>(ctx);
      p->monitorTao(tao);
    } catch (...) {
      MPI_Comm com = MPI_COMM_SELF;
      PetscErrorCode ierr = PetscObjectGetComm((PetscObject)tao, &com); CHKERRQ(ierr);
      handle_fatal_errors(com);
      SETERRQ(com, 1, "A PISM callback failed");
    }
    return 0;
  }
};

//! \brief Adaptor to connect a TAO objective function callback to a C++ object method.
/*! The TAO library interfaces with user code via C-style callback functions.
  This class makes it convenient to associate a TAO VariableBounds callback
  with a C++ object method. To assign 
  \code
  void MyObject::getVariableBounds(Tao tao,Vec lo, Vec hi);
  \endcode

  as the objective function to a `Tao` `tao`, 

  \code
  MyObject obj;
  TaoGetVariableBoundsCallback<MyObject>::connect(tao,obj);
  \endcode

  The method name `getVariableBounds` for the callback is hard-coded.
  See TaoObjGradCallback for a technique to allow 
  the method name to be specified (at the expense of a little more cumbersome code).
*/
template<class Problem>
class TaoGetVariableBoundsCallback {
public:
  static void connect(Tao tao, Problem &p) {
    PetscErrorCode ierr;
    ierr = TaoSetVariableBoundsRoutine(tao,
                                       TaoGetVariableBoundsCallback<Problem>::callback,
                                       &p);
    PISM_CHK(ierr, "TaoSetVariableBoundsRoutine");
  }

protected:
  static PetscErrorCode callback(Tao tao, Vec lo, Vec hi, void *ctx) {
    try {
      Problem *p = reinterpret_cast<Problem *>(ctx);
      p->getVariableBounds(tao,lo,hi);
    } catch (...) {
      MPI_Comm com = MPI_COMM_SELF;
      PetscErrorCode ierr = PetscObjectGetComm((PetscObject)tao, &com); CHKERRQ(ierr);
      handle_fatal_errors(com);
      SETERRQ(com, 1, "A PISM callback failed");
    }
    return 0;
  }
};

//! \brief Adaptor to connect a TAO objective gradient callback to a C++ object method.
/*! The TAO library interfaces with user code via C-style callback functions.
  This class makes it convenient to associate a TAO Objective Gradient callback
  with a C++ object method. To assign 
  \code
  void MyObject::evaluateGradient(Tao tao,Vec x, Vec gradient);
  \endcode

  as the objective function to a `Tao` `tao`, 

  \code
  MyObject obj;
  TaoGradientCallback<MyObject>::connect(tao,obj);
  \endcode

  The method name `evaluateGradient` for the callback is hard-coded.
  See TaoObjGradCallback for a technique to allow 
  the method name to be specified (at the expense of a little more cumbersome code).
*/
template<class Problem>
class TaoGradientCallback {
public:
  static void connect(Tao tao, Problem &p) {
    PetscErrorCode ierr;
    ierr = TaoSetGradientRoutine(tao,
                                 TaoGradientCallback<Problem>::callback,
                                 &p);
    PISM_CHK(ierr, "TaoSetGradientRoutine");
  }

protected:
  static PetscErrorCode callback(Tao tao, Vec x, Vec gradient, void *ctx) {
    try {
      Problem *p = reinterpret_cast<Problem *>(ctx);
      p->evaluateGradient(tao,x,gradient);
    } catch (...) {
      MPI_Comm com = MPI_COMM_SELF;
      PetscErrorCode ierr = PetscObjectGetComm((PetscObject)tao, &com); CHKERRQ(ierr);
      handle_fatal_errors(com);
      SETERRQ(com, 1, "A PISM callback failed");
    }
    return 0;
  }
};

//! \brief Adaptor to connect a TAO objective function callback to a C++ object method.
/*! The TAO library interfaces with user code via C-style callback functions.
  This class makes it convenient to associate a TAO convergence monitoring callback
  with a C++ object method. To assign 
  \code
  void MyObject::convergenceTest(Tao tao);
  \endcode

  as the convergence test function to a `Tao` `tao`, 

  \code
  MyObject obj;
  TaoConvergenceCallback<MyObject>::connect(tao,obj);
  \endcode

  The method name `convergenceTest` for the callback is hard-coded.
  See TaoObjGradCallback for a technique to allow 
  the method name to be specified (at the expense of a little more cumbersome code).
*/
template<class Problem>
class TaoConvergenceCallback {
public:
  static void connect(Tao tao, Problem &p) {
    PetscErrorCode ierr;
    ierr = TaoSetConvergenceTest(tao,
                                 TaoConvergenceCallback<Problem>::callback,
                                 &p);
    PISM_CHK(ierr, "TaoSetConvergenceTest");
  }
protected:
  static PetscErrorCode callback(Tao tao, void *ctx) {
    try {
      Problem *p = reinterpret_cast<Problem *>(ctx);
      p->convergenceTest(tao);
    } catch (...) {
      MPI_Comm com = MPI_COMM_SELF;
      PetscErrorCode ierr = PetscObjectGetComm((PetscObject)tao, &com); CHKERRQ(ierr);
      handle_fatal_errors(com);
      SETERRQ(com, 1, "A PISM callback failed");
    }
    return 0;
  }
};


//! \brief Adaptor to connect a TAO objective and gradient function callback to a C++ object method.
/*! The TAO library interfaces with user code via C-style callback functions.
  This class makes it convenient to associate a TAO combined objective value and gradient 
  callback with a C++ object method. To assign 
  \code
  void MyObject::someObjectiveFunction(Tao tao,Vec x, double *value, Vec gradient);
  \endcode

  as the convergence test function to a `Tao` `tao`, 

  \code
  MyObject obj;
  typedef TaoObjGradCallback<MyObject,&MyObject::someObjectiveFunction> ObjGradCallback;
  ObjGradCallback::connect(tao,obj);
  \endcode

  Note that the method name for the callback must be specified explicitly via a template argument.
*/
template<class Problem, void (Problem::*Callback)(Tao,Vec,double*,Vec) >
class TaoObjGradCallback {
public:
  static void connect(Tao tao, Problem &p) {
    PetscErrorCode ierr;
    ierr = TaoSetObjectiveAndGradientRoutine(tao,
                                             TaoObjGradCallback<Problem,Callback>::callback,
                                             &p);
    PISM_CHK(ierr, "TaoSetObjectiveAndGradientRoutine");
  }
protected:
  static PetscErrorCode callback(Tao tao, Vec x, double *value, Vec gradient, void *ctx) {
    try {
      Problem *p = reinterpret_cast<Problem *>(ctx);
      (p->*Callback)(tao,x,value,gradient);
    } catch (...) {
      MPI_Comm com = MPI_COMM_SELF;
      PetscErrorCode ierr = PetscObjectGetComm((PetscObject)tao, &com); CHKERRQ(ierr);
      handle_fatal_errors(com);
      SETERRQ(com, 1, "A PISM callback failed");
    }
    return 0;
  }
};

//! \brief Adaptor to connect a TAO objective function callback to a C++ object method.
/*! The TAO library interfaces with user code via C-style callback functions.
  This class makes it convenient to associate a TAO Linearly Constrained Augmented Lagrangian (LCL)
  callbacks with C++ object methods. To assign 
  \code
  void MyObject::evaluateConstraints(Tao tao,Vec x,Vec c);
  void MyObject::evaluateConstraintsJacobianState(Tao tao, Vec x, Mat *J, Mat *Jpc, Mat *Jinv, MatStructure *structure);
  void MyObject::evaluateConstraintsJacobianDesign(Tao tao, Vec x, Mat *J);
  \endcode
  as the LCL callbacks to a `Tao` `tao`,

  \code
  MyObject obj;
  TaoLCLCallback<MyObject>::connect(tao,obj);
  \endcode

  The method names for the callback (`evaluateConstraints`, etc.) are hard-coded.
*/
template<class Problem>
class TaoLCLCallbacks {
public:
  static void connect(Tao tao, Problem &p, Vec c, Mat Jc, Mat Jd, Mat Jcpc=NULL, Mat Jcinv=NULL) {
    PetscErrorCode ierr;

    ierr = TaoSetConstraintsRoutine(tao, c,
                                    TaoLCLCallbacks<Problem>::constraints_callback, &p);
    PISM_CHK(ierr, "TaoSetConstraintsRoutine");

    if (Jcpc == NULL) {
      Jcpc = Jc;
    }

    ierr = TaoSetJacobianStateRoutine(tao, Jc, Jcpc, Jcinv,
                                      TaoLCLCallbacks<Problem>::jacobian_state_callback, &p);
    PISM_CHK(ierr, "TaoSetJacobianStateRoutine");

    ierr = TaoSetJacobianDesignRoutine(tao, Jd,
                                       TaoLCLCallbacks<Problem>::jacobian_design_callback, &p);
    PISM_CHK(ierr, "TaoSetJacobianDesignRoutine");
  }
protected:
  static PetscErrorCode constraints_callback(Tao tao, Vec x, Vec c, void*ctx) {
    try {
      Problem *p = reinterpret_cast<Problem *>(ctx);
      p->evaluateConstraints(tao, x, c);
    } catch (...) {
      MPI_Comm com = MPI_COMM_SELF;
      PetscErrorCode ierr = PetscObjectGetComm((PetscObject)tao, &com); CHKERRQ(ierr);
      handle_fatal_errors(com);
      SETERRQ(com, 1, "A PISM callback failed");
    }
    return 0;
  }

  static PetscErrorCode jacobian_state_callback(Tao tao, Vec x, Mat J, Mat Jpc,
                                                Mat Jinv, void*ctx) {
    try {
      Problem *p = reinterpret_cast<Problem *>(ctx);
      // The MatStructure argument is not used in PETSc 3.5, but I want
      // to preserve the signature of
      // evaluateConstraintsJacobianState(...) for now -- (CK)
      MatStructure structure;
      p->evaluateConstraintsJacobianState(tao, x, J, Jpc, Jinv, &structure);
    } catch (...) {
      MPI_Comm com = MPI_COMM_SELF;
      PetscErrorCode ierr = PetscObjectGetComm((PetscObject)tao, &com); CHKERRQ(ierr);
      handle_fatal_errors(com);
      SETERRQ(com, 1, "A PISM callback failed");
    }
    return 0;
  }

  static PetscErrorCode jacobian_design_callback(Tao tao, Vec x, Mat J, void*ctx) {
    try {
      Problem *p = reinterpret_cast<Problem *>(ctx);
      p->evaluateConstraintsJacobianDesign(tao, x, J);
    } catch (...) {
      MPI_Comm com = MPI_COMM_SELF;
      PetscErrorCode ierr = PetscObjectGetComm((PetscObject)tao, &com); CHKERRQ(ierr);
      handle_fatal_errors(com);
      SETERRQ(com, 1, "A PISM callback failed");
    }
    return 0;
  }
};

} // end of namespace taoutil
} // end of namespace pism

#endif /* end of include guard: TAOUTIL_HH_W42GJNRO */
