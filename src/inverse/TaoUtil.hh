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

#ifndef TAOUTIL_HH_W42GJNRO
#define TAOUTIL_HH_W42GJNRO

#include <tao.h>
#include <string>
#include "pism_const.hh"

extern const char *const* TaoConvergedReasons;

class TaoInitializer {
public:
  TaoInitializer(int *argc, char ***argv, char *file, char *help);
  TaoInitializer(int *argc, char ***argv, char *help);
  TaoInitializer(int *argc, char ***argv);
  ~TaoInitializer();
private:
};

//typedef PetscObjectGuard<TaoSolver,TaoDestroy> TaoSolverGuard;
typedef enum {
  kTaoBasicCallbacksCombined,
  kTaoBasicCallbacksSeparated
} TaoBasicCallbackType;


template<class Problem>
class TaoBasicSolver {
public:
    
  TaoBasicSolver(MPI_Comm comm, const char* tao_type, Problem &prob, TaoBasicCallbackType cbType = kTaoBasicCallbacksCombined ):
  m_comm(comm), m_problem(prob)
   {
    PetscErrorCode ierr;
    ierr = this->construct(m_comm, tao_type);
    if(ierr) {
      PetscPrintf(m_comm, "FATAL ERROR: TaoBasicProblem allocation failed.\n");
      PISMEnd();
    }    
  }
  
  virtual ~TaoBasicSolver() {
    PetscErrorCode ierr;
    ierr = this->destruct();
    if(ierr) {
      PetscPrintf(m_comm, "FATAL ERROR: TaoBasicProblem deallocation failed.\n");
      PISMEnd();
    }
  };

  virtual PetscErrorCode solve(bool &success) {
    PetscErrorCode ierr;
    m_success=false; 
    m_reason = TAO_CONTINUE_ITERATING;
    m_reasonDescription.clear();

    /* Solve the application */ 
    ierr = TaoSetInitialVector(m_tao, m_problem.formInitialGuess()); CHKERRQ(ierr);
    ierr = TaoSolve(m_tao); CHKERRQ(ierr);  

    ierr = TaoGetTerminationReason(m_tao, &m_reason); CHKERRQ(ierr);

    if(reason>0){
      m_success = true;
    } 
    m_reasonDescription = TaoConvergedReasons[reason];
    success =  m_success;
    return 0;
  }

  virtual TaoSolverTerminationReason reason() {
    return m_reason;
  }
  virtual std::string const &reasonDescription() {
    return m_reasonDescription;
  }

  virtual Problem &problem() {
    return m_problem;
  }

protected:

  virtual PetscErrorCode construct(const char* tao_type, TaoBasicCallbackType cbType  ) {
    PetscErrorCode ierr;
    ierr = TaoCreate(m_comm ,&m_tao); CHKERRQ(ierr); 
    ierr = TaoSetType(m_tao,tao_type); CHKERRQ(ierr);    
    if(cbType == kTaoBasicCallbacksCombined) {
      // ierr = TaoObjectiveAndGradientCallbacks::connectCombined(m_tao,&m_problem); CHKERRQ(ierr);
    } else {
      // ierr = TaoObjectiveAndGradientCallbacks::connectSeparate(m_tao,&m_problem); CHKERRQ(ierr);    
    }
    ierr = TaoSetFromOptions(m_tao); CHKERRQ(ierr);
    return 0;
  }

  virtual PetscErrorCode destruct() {
    PetscErrorCode ierr;
    if(m_tao) {
      ierr = TaoDestroy(&m_tao); CHKERRQ(ierr);
    }
    return 0; 
  }
  
  // virtual void evaluateGradient(TaoSolver *tao, Vec x, PetscReal *value);
  
  MPI_Comm m_comm;
  TaoSolver m_tao;
  Problem  &m_problem;

  bool m_success;
  TaoSolverTerminationReason m_reason;
  std::string m_reasonDescription;
};

/*
A Problem for a TaoBasicSolver must implement:
virtual Vec formInitialGuess() = 0;
and either
virtual void evaluateObjective(TaoSolver tao, Vec x, PetscReal *value);
virtual void evaluateGradient(TaoSolver tao, Vec x, Vec gradient);
or
virtual void evaluateObjectiveAndGradient(TaoSolver tao, Vec x, PetscReal *value, Vec gradient);
*/


template<class Problem>
class TaoCombinedObjectiveAndGradientCallback {
public:

  static PetscErrorCode attach(TaoSolver tao, Problem &p) {
    PetscErrorCode ierr;
    ierr = TaoSetObjectiveAndGradientRoutine(tao,
      TaoCombinedObjectiveAndGradientCallback<Problem>::evaluateObjectiveAndGradientCallback,
      &p ); CHKERRQ(ierr);
  }
  
protected:

  static PetscErrorCode evaluateObjectiveAndGradientCallback(TaoSolver tao,
                                 Vec x, PetscReal *value, Vec gradient, Problem *ctx ) {
    PetscErrorCode ierr;
    Problem *p = reinterpret_cast<Problem *>(ctx);
    ierr = p->evaluateObjectiveAndGradient(tao,x,value,gradient); CHKERRQ(ierr);
    return 0;
  }
};

template<class Problem>
class TaoSeparateObjectiveAndGradientCallbacks {
public:

  static PetscErrorCode attach(TaoSolver tao, Problem &p) {
    PetscErrorCode ierr;
    ierr = TaoSetObjectiveRoutine(tao,
      TaoSeparateObjectiveAndGradientCallbacks<Problem>::evaluateObjective
      &p ); CHKERRQ(ierr);

    ierr = TaoSetGradientRoutine(tao,
      TaoSeparateObjectiveAndGradientCallbacks<Problem>::evaluateGradient
      &p ); CHKERRQ(ierr);
  }

protected:

  static PetscErrorCode evaluateObjectiveCallback(TaoSolver tao,
                                 Vec x, PetscReal *value, void *ctx ) {
    PetscErrorCode ierr;
    Problem *p = reinterpret_cast<Problem *>(ctx);
    ierr = p->evaluateObjective(tao,x,value); CHKERRQ(ierr);
    return 0;
  }

  static PetscErrorCode evaluateGradientCallback(TaoSolver tao,
                                 Vec x, Vec gradient, void *ctx ) {
    PetscErrorCode ierr;
    Problem *p = reinterpret_cast<Problem *>(ctx);
    ierr = p->evaluateGradient(tao,x,gradient); CHKERRQ(ierr);
    return 0;
  }
};

template<class Problem>
class TaoLCLCallbacks {
public:
  static PetscErrorCode connect(TaoSolver tao, Problem &p, Vec c, Mat Jc, Mat Jd, Mat Jcpc=NULL, Mat Jcinv=NULL, Mat Jdpc=NULL) {
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
