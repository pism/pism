// Copyright (C) 2012, 2013  David Maxwell
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

#ifndef TAOUTIL_HH_W42GJNRO
#define TAOUTIL_HH_W42GJNRO

#include <tao.h>
#include <string>
#include "pism_const.hh"
#include "TerminationReason.hh"

extern const char *const* TaoConvergedReasons;

//! Class to initialize the TAO library using the Resource Allocation Is Initialization (RAII) paradigm.
/*! Declare a TaoInitializer on the stack to initialize the library in, e.g. \c main. When its destructor is called,
    the TAO library will be finalized.
*/
class TaoInitializer {
public:
  TaoInitializer(int *argc, char ***argv, char *file, char *help);
  TaoInitializer(int *argc, char ***argv, char *help);
  TaoInitializer(int *argc, char ***argv);
  ~TaoInitializer();
private:
};

//! Encapsulate TAO's TaoSolverTerminationReason codes as a PISM TerminationReason.
class TAOTerminationReason: public TerminationReason {
public:
  TAOTerminationReason( TaoSolverTerminationReason r);
  virtual void get_description( std::ostream &desc,int indent_level=0);
};

//! \brief An interface for solving an optimization problem with TAO where the
//! problem itself is defined by a separate Problem class.
/*!  The primary interface to a TAO optimization problem is mediated by a PETSc-style
\c TaoSolver object. The PISM TaoBasicSolver C++ class wraps a \c TaoSolver and some of 
its initialization boilierplate, and allows a separate class to define the function to be minimized.

To use a TaoBasicSolver you create a \c Problem class that defines the objective function and initial
guess, as well any auxilliary callbacks desired.  The Problem class must define a

\code
PetscErrorCode Problem::connect(TaoSolver solver);
\endcode

method which gives the \c Problem an opportunity to register its methods as callbacks to the solver, 
perhaps taking advantage of the various \c TaoFooCallback classes provided in TaoUtil.hh to facilitate this.
For example, a problem class MyProblem that did nothing more than register a combined objective/gradient
callback could define

\code
PetscErrorCode MyProblem::connect(TaoSolver tao) {
  PetscErrorCode ierr;
  typedef TaoObjGradCallback<Problem,&MyProblem::evaluateObjectiveAndGradient> ObjGradCallback; 
  ierr = ObjGradCallback::connect(tao,*this); CHKERRQ(ierr);
  return 0;
}
\endcode

In addition to the \c connect method, a \c Problem must define 
\code
PetscErrorCode MyProblem::formInitialGuess(Vec *v, TerminationReason::Ptr &reason)
\endcode
which allows the problem to set the initial guess for optimization. If the minimization
is successful, the solution will be found in the same vector that was returned by this method.

Assuming a \c MyProblem called \c problem has been constructed, solution
of the minimization is done using, for example, the TAO algorithm
tao_cg:

\code
TaoBasicSolver<MyProblem> solver(com,"tao_cg",problem);

TerminationReason::Ptr reason;
solver.solve(reason);

if(reason->succeeded()) {
  printf("Success: %s\n",reason->description().c_str());
} else {
  printf("Failure: %s\n",reason->description().c_str());
}
\endcode
*/
template<class Problem>
class TaoBasicSolver {
public:
    
  //! Construct a solver to solve \c prob using TAO algorithm \c tao_type.
  TaoBasicSolver(MPI_Comm comm, const char* tao_type, Problem &prob):
  m_comm(comm), m_problem(prob)
   {
    PetscErrorCode ierr;
    ierr = this->construct(tao_type);
    if(ierr) {
      CHKERRCONTINUE(ierr);
      PetscPrintf(m_comm, "FATAL ERROR: TaoBasicProblem allocation failed.\n");
      PISMEnd();
    }    
  }
  
  virtual ~TaoBasicSolver() {
    PetscErrorCode ierr;
    ierr = this->destruct();
    if(ierr) {
      CHKERRCONTINUE(ierr);
      PetscPrintf(m_comm, "FATAL ERROR: TaoBasicProblem deallocation failed.\n");
      PISMEnd();
    }
  };

  //! Solve the minimization problem.
  virtual PetscErrorCode solve(TerminationReason::Ptr &reason) {
    PetscErrorCode ierr;

    /* Solve the application */ 
    Vec x0;
    ierr = m_problem.formInitialGuess(&x0,reason); CHKERRQ(ierr);
    if(reason->failed()) {
      TerminationReason::Ptr root_cause = reason;
      reason.reset(new GenericTerminationReason(-1,"Unable to form initial guess"));
      reason->set_root_cause(root_cause);
      return 0;
    }
    ierr = TaoSetInitialVector(m_tao, x0); CHKERRQ(ierr);
    ierr = TaoSolve(m_tao); CHKERRQ(ierr);  

    TaoSolverTerminationReason tao_reason;
    ierr = TaoGetTerminationReason(m_tao, &tao_reason); CHKERRQ(ierr);
    reason.reset(new TAOTerminationReason(tao_reason));

    return 0;
  }

  virtual Problem &problem() {
    return m_problem;
  }

protected:

  //! Initialize the TaoSolver and allow the Problem to connect its callbacks.
  virtual PetscErrorCode construct(const char* tao_type) {
    PetscErrorCode ierr;
    ierr = TaoCreate(m_comm ,&m_tao); CHKERRQ(ierr); 
    ierr = TaoSetType(m_tao,tao_type); CHKERRQ(ierr);    
    ierr = m_problem.connect(m_tao); CHKERRQ(ierr);
    ierr = TaoSetFromOptions(m_tao); CHKERRQ(ierr);
    return 0;
  }

  //! Finalize the TaoSolver.
  virtual PetscErrorCode destruct() {
    PetscErrorCode ierr;
    if(m_tao) {
      ierr = TaoDestroy(&m_tao); CHKERRQ(ierr);
    }
    return 0; 
  }
  
  MPI_Comm m_comm;
  TaoSolver m_tao;
  Problem  &m_problem;
};


//! \brief Adaptor to connect a TAO objective function callback to a C++ object method.
/*! The TAO library interfaces with user code via C-style callback functions.
This class makes it convenient to associate a TAO Objective callback
with a C++ object method. To assign 
\code
PetscErrorCode MyObject::evaluateObjective(TaoSolver tao,Vec x, PetscReal *value);
\endcode

as the objective function to a \c TaoSolver \c tao, 

\code
MyObject obj;
TaoObjectiveCallback<MyObject>::connect(tao,obj);
\endcode

The method name \c evaluateObjective for the callback is hard-coded.
See TaoObjGradCallback for a technique to allow 
the method name to be specified (at the expense of a little more cumbersome code).
*/
template<class Problem>
class TaoObjectiveCallback {
public:

  static PetscErrorCode connect(TaoSolver tao, Problem &p) {
    PetscErrorCode ierr;
    ierr = TaoSetObjectiveRoutine(tao,
      TaoObjectiveCallback<Problem>::evaluateObjectiveCallback,
      &p ); CHKERRQ(ierr);
    return 0;
  }

protected:

  static PetscErrorCode evaluateObjectiveCallback(TaoSolver tao,
                                 Vec x, PetscReal *value, void *ctx ) {
    PetscErrorCode ierr;
    Problem *p = reinterpret_cast<Problem *>(ctx);
    ierr = p->evaluateObjective(tao,x,value); CHKERRQ(ierr);
    return 0;
  }
};


//! \brief Adaptor to connect a TAO monitoring callback to a C++ object method.
/*! The TAO library interfaces with user code via C-style callback functions.
This class makes it convenient to associate a TAO Monitor callback
with a C++ object method. To assign 
\code
PetscErrorCode MyObject::monitorTao(TaoSolver tao)
\endcode

as the objective function to a \c TaoSolver \c tao, 

\code
MyObject obj;
TaoMonitorCallback<MyObject>::connect(tao,obj);
\endcode

The method name \c monitorTao for the callback is hard-coded.
See TaoObjGradCallback for a technique to allow 
the method name to be specified (at the expense of a little more cumbersome code).
*/
template<class Problem>
class TaoMonitorCallback {
public:

  static PetscErrorCode connect(TaoSolver tao, Problem &p) {
    PetscErrorCode ierr;
    ierr = TaoSetMonitor(tao,
      TaoMonitorCallback<Problem>::monitorTao,
      &p, NULL ); CHKERRQ(ierr);
    return 0;
  }

protected:

  static PetscErrorCode monitorTao(TaoSolver tao, void *ctx ) {
    PetscErrorCode ierr;
    Problem *p = reinterpret_cast<Problem *>(ctx);
    ierr = p->monitorTao(tao); CHKERRQ(ierr);
    return 0;
  }
};

//! \brief Adaptor to connect a TAO objective function callback to a C++ object method.
/*! The TAO library interfaces with user code via C-style callback functions.
This class makes it convenient to associate a TAO VariableBounds callback
with a C++ object method. To assign 
\code
PetscErrorCode MyObject::getVariableBounds(TaoSolver tao,Vec lo, Vec hi);
\endcode

as the objective function to a \c TaoSolver \c tao, 

\code
MyObject obj;
TaoGetVariableBoundsCallback<MyObject>::connect(tao,obj);
\endcode

The method name \c getVariableBounds for the callback is hard-coded.
See TaoObjGradCallback for a technique to allow 
the method name to be specified (at the expense of a little more cumbersome code).
*/
template<class Problem>
class TaoGetVariableBoundsCallback {
public:

  static PetscErrorCode connect(TaoSolver tao, Problem &p) {
    PetscErrorCode ierr;
    ierr = TaoSetVariableBoundsRoutine(tao,
      TaoGetVariableBoundsCallback<Problem>::getVariableBounds,
      &p); CHKERRQ(ierr);
    return 0;
  }

protected:

  static PetscErrorCode getVariableBounds(TaoSolver tao, Vec lo, Vec hi, void *ctx ) {
    PetscErrorCode ierr;
    Problem *p = reinterpret_cast<Problem *>(ctx);
    ierr = p->getVariableBounds(tao,lo,hi); CHKERRQ(ierr);
    return 0;
  }
};

//! \brief Adaptor to connect a TAO objective gradient callback to a C++ object method.
/*! The TAO library interfaces with user code via C-style callback functions.
This class makes it convenient to associate a TAO Objective Gradient callback
with a C++ object method. To assign 
\code
PetscErrorCode MyObject::evaluateGradient(TaoSolver tao,Vec x, Vec gradient);
\endcode

as the objective function to a \c TaoSolver \c tao, 

\code
MyObject obj;
TaoGradientCallback<MyObject>::connect(tao,obj);
\endcode

The method name \c evaluateGradient for the callback is hard-coded.
See TaoObjGradCallback for a technique to allow 
the method name to be specified (at the expense of a little more cumbersome code).
*/
template<class Problem>
class TaoGradientCallback {
public:

  static PetscErrorCode connect(TaoSolver tao, Problem &p) {
    PetscErrorCode ierr;
    ierr = TaoSetGradientRoutine(tao,
      TaoGradientCallback<Problem>::evaluateGradient,
      &p ); CHKERRQ(ierr);
    return 0;
  }

protected:

  static PetscErrorCode evaluateGradient(TaoSolver tao,
                                 Vec x, Vec gradient, void *ctx ) {
    PetscErrorCode ierr;
    Problem *p = reinterpret_cast<Problem *>(ctx);
    ierr = p->evaluateGradient(tao,x,gradient); CHKERRQ(ierr);
    return 0;
  }
};

//! \brief Adaptor to connect a TAO objective function callback to a C++ object method.
/*! The TAO library interfaces with user code via C-style callback functions.
This class makes it convenient to associate a TAO convergence monitoring callback
with a C++ object method. To assign 
\code
PetscErrorCode MyObject::convergenceTest(TaoSolver tao);
\endcode

as the convergence test function to a \c TaoSolver \c tao, 

\code
MyObject obj;
TaoConvergenceCallback<MyObject>::connect(tao,obj);
\endcode

The method name \c convergenceTest for the callback is hard-coded.
See TaoObjGradCallback for a technique to allow 
the method name to be specified (at the expense of a little more cumbersome code).
*/
template<class Problem>
class TaoConvergenceCallback {
public:

  static PetscErrorCode connect(TaoSolver tao, Problem &p) {
    PetscErrorCode ierr;
    ierr = TaoSetConvergenceTest(tao,
      TaoConvergenceCallback<Problem>::convergenceTestCallback,
      &p ); CHKERRQ(ierr);
    return 0;
  }

protected:

  static PetscErrorCode convergenceTestCallback(TaoSolver tao, void *ctx ) {
    PetscErrorCode ierr;
    Problem *p = reinterpret_cast<Problem *>(ctx);
    ierr = p->convergenceTest(tao); CHKERRQ(ierr);
    return 0;
  }
};


//! \brief Adaptor to connect a TAO objective and gradient function callback to a C++ object method.
/*! The TAO library interfaces with user code via C-style callback functions.
This class makes it convenient to associate a TAO combined objective value and gradient 
callback with a C++ object method. To assign 
\code
PetscErrorCode MyObject::someObjectiveFunction(TaoSolver tao,Vec x, PetscReal *value, Vec gradient);
\endcode

as the convergence test function to a \c TaoSolver \c tao, 

\code
MyObject obj;
typedef TaoObjGradCallback<MyObject,&MyObject::someObjectiveFunction> ObjGradCallback;
ObjGradCallback::connect(tao,obj);
\endcode

Note that the method name for the callback must be specified explicitly via a template argument.
*/
template<class Problem, PetscErrorCode (Problem::*Callback)(TaoSolver,Vec,PetscReal*,Vec) >
class TaoObjGradCallback {
public:

  static PetscErrorCode connect(TaoSolver tao, Problem &p) {
    PetscErrorCode ierr;
    ierr = TaoSetObjectiveAndGradientRoutine(tao,
      TaoObjGradCallback<Problem,Callback>::evaluateObjectiveAndGradientCallback,
      &p ); CHKERRQ(ierr);
    return 0;
  }
  
protected:

  static PetscErrorCode evaluateObjectiveAndGradientCallback(TaoSolver tao,
                                 Vec x, PetscReal *value, Vec gradient, void *ctx ) {
    PetscErrorCode ierr;
    Problem *p = reinterpret_cast<Problem *>(ctx);
    ierr = (p->*Callback)(tao,x,value,gradient); CHKERRQ(ierr);
    return 0;
  }
};

//! \brief Adaptor to connect a TAO objective function callback to a C++ object method.
/*! The TAO library interfaces with user code via C-style callback functions.
This class makes it convenient to associate a TAO Linearly Constrained Augmented Lagrangian (LCL)
callbacks with C++ object methods. To assign 
\code
PetscErrorCode MyObject::evaluateConstraints(TaoSolver tao,Vec x,Vec c);
PetscErrorCode MyObject::evaluateConstraintsJacobianState(TaoSolver tao, Vec x, Mat *J, Mat *Jpc, Mat *Jinv, MatStructure *structure);
PetscErrorCode MyObject::evaluateConstraintsJacobianDesign(TaoSolver tao, Vec x, Mat *J);
\endcode
as the LCL callbacks to a \c TaoSolver \c tao, 

\code
MyObject obj;
TaoLCLCallback<MyObject>::connect(tao,obj);
\endcode

The method names for the callback (\c evaluateConstraints, etc.) are hard-coded.
*/
template<class Problem>
class TaoLCLCallbacks {
public:
  static PetscErrorCode connect(TaoSolver tao, Problem &p, Vec c, Mat Jc, Mat Jd, Mat Jcpc=NULL, Mat Jcinv=NULL) {
    PetscErrorCode ierr;
    ierr = TaoSetConstraintsRoutine(tao,c,TaoLCLCallbacks<Problem>::evaluateConstraintsCallback,&p); CHKERRQ(ierr);
    if(Jcpc==NULL) Jcpc = Jc;
    ierr = TaoSetJacobianStateRoutine(tao,Jc,Jcpc,Jcinv,TaoLCLCallbacks<Problem>::evaluateJacobianStateCallback,&p); CHKERRQ(ierr);
    ierr = TaoSetJacobianDesignRoutine(tao,Jd,TaoLCLCallbacks<Problem>::evaluateJacobianDesignCallback,&p); CHKERRQ(ierr);
    return 0;
  }
protected:
  static PetscErrorCode evaluateConstraintsCallback(TaoSolver tao, Vec x,Vec c, void*ctx){
    PetscErrorCode ierr;
    Problem *p = reinterpret_cast<Problem *>(ctx);
    ierr = p->evaluateConstraints(tao,x,c); CHKERRQ(ierr);
    return 0;
  }

  static PetscErrorCode evaluateJacobianStateCallback(TaoSolver tao, Vec x, Mat *J, Mat *Jpc, Mat *Jinv, MatStructure *structure, void*ctx){
    PetscErrorCode ierr;
    Problem *p = reinterpret_cast<Problem *>(ctx);
    ierr = p->evaluateConstraintsJacobianState(tao,x,J,Jpc,Jinv,structure);
      CHKERRQ(ierr);
    return 0;
  }

  static PetscErrorCode evaluateJacobianDesignCallback(TaoSolver tao, Vec x, Mat *J, void*ctx){
    PetscErrorCode ierr;
    Problem *p = reinterpret_cast<Problem *>(ctx);
    ierr = p->evaluateConstraintsJacobianDesign(tao,x,J); CHKERRQ(ierr);
    return 0;
  }
};

#endif /* end of include guard: TAOUTIL_HH_W42GJNRO */
