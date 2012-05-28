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

#include "TaoUtil.hh"

class ObjectiveFunction {
public:
  ObjectiveFunction(PetscInt N);
  ~ObjectiveFunction();

  PetscErrorCode create(PetscInt N);

  PetscErrorCode evaluateObjectiveAndGradient(TaoSolver /*tao*/, Vec x, PetscReal *val, Vec gradient);
  
  Vec formInitialGuess() {
    return m_x;
  }

  Vec get_x() {
    return m_x;
  }

private:
  Vec m_x;
  Vec m_x0;
  Vec m_work;
};

ObjectiveFunction::ObjectiveFunction(PetscInt N) {
  PetscErrorCode ierr = this->create(N);
  if(ierr) {
    CHKERRCONTINUE(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"ObjectiveFunction creation failed.");
  }
}

ObjectiveFunction::~ObjectiveFunction() {
  VecDestroy(&m_x);
  VecDestroy(&m_x0);
  VecDestroy(&m_work);
}

PetscErrorCode ObjectiveFunction::create(PetscInt N) {
  PetscErrorCode ierr;
  ierr = VecCreate(PETSC_COMM_WORLD, &m_x0); CHKERRQ(ierr);
  ierr = VecSetFromOptions(m_x0); CHKERRQ(ierr);
  ierr = VecSetSizes(m_x0,PETSC_DECIDE,N); CHKERRQ(ierr);
  ierr = VecSet(m_x0,1); CHKERRQ(ierr);
  
  PetscInt i0, i1;
  ierr = VecGetOwnershipRange(m_x0,&i0,&i1);
  PetscInt n = i1-i0;

  ierr = VecCreate(PETSC_COMM_WORLD, &m_x); CHKERRQ(ierr);
  ierr = VecSetFromOptions(m_x); CHKERRQ(ierr);
  ierr = VecSetSizes(m_x,n,PETSC_DECIDE); CHKERRQ(ierr);
  ierr = VecSet(m_x,0); CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD, &m_work); CHKERRQ(ierr);
  ierr = VecSetFromOptions(m_work); CHKERRQ(ierr);
  ierr = VecSetSizes(m_work,n,PETSC_DECIDE); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode ObjectiveFunction::evaluateObjectiveAndGradient(TaoSolver, Vec x, PetscReal *val, Vec gradient) {

  PetscErrorCode ierr;
  ierr = VecCopy(x,m_work); CHKERRQ(ierr);
  ierr = VecAXPY(m_work,-1,m_x0); CHKERRQ(ierr);
  ierr = VecCopy(m_work,gradient); CHKERRQ(ierr);
  ierr = VecScale(gradient,2.); CHKERRQ(ierr);

  ierr = VecNorm(m_work,NORM_2,val);
  return 0;
}

// typedef TaoBasicSolver<ObjectiveFunction, TaoCombinedObjectiveAndGradientCallback<ObjectiveFunction> > TaoObjGradSolver;

typedef TaoSolverCategory<ObjectiveFunction>::ObjGrad TaoObjGradSolver;

int main(int argc, char** argv) {
  PetscInitialize(&argc,&argv,NULL,NULL);
  {
    TaoInitializer ti(&argc,&argv);

    PetscBool flag;
    PetscInt N = 2;
    PetscErrorCode ierr;
    ierr = PetscOptionsInt("-N","Problem size","",N,&N,&flag); CHKERRQ(ierr);

    ObjectiveFunction f(N);
    TaoObjGradSolver solver(PETSC_COMM_WORLD,"tao_cg", f);
  
    bool success;
    ierr = solver.solve(success); CHKERRQ(ierr);
    if(success) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Success: %s\n",solver.reasonDescription().c_str()); CHKERRQ(ierr);
      ierr = VecView(f.get_x(),PETSC_VIEWER_STDOUT_WORLD);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Failure: %s\n",solver.reasonDescription().c_str()); CHKERRQ(ierr);    
    }
  }
  PetscFinalize();
}