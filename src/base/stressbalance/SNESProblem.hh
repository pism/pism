// Copyright (C) 2011 David Maxwell
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


#ifndef _SNESPROBLEM_H_
#define _SNESPROBLEM_H_

#include "grid.hh"
#include "iceModelVec.hh"

template<int DOF, class U> class SNESProblem{
public:
  SNESProblem(IceGrid &g);
  
  virtual ~SNESProblem();

  virtual PetscErrorCode solve();
  
  virtual const char *name();

  virtual Vec solution()
  {
    return m_X;
  }

protected:

  virtual PetscErrorCode initialize();

  virtual PetscErrorCode setFromOptions();

  virtual PetscErrorCode finalize();
  
  virtual PetscErrorCode compute_local_function(DALocalInfo *info, const U **xg, U **yg) = 0;  

  virtual PetscErrorCode compute_local_jacobian(DALocalInfo *info, const U **x,  Mat B) = 0;

  IceGrid     &m_grid;

  Mat          m_J;
  Vec          m_X;
  Vec          m_F;
  SNES         m_snes;
  DA           m_DA;

private:

  struct SNESProblemCallbackData {
    DA           da;
    SNESProblem<DOF,U> *solver;
  };

  SNESProblemCallbackData m_callbackData;

  static PetscErrorCode LocalFunction(DALocalInfo *info, const U **x, U **f, SNESProblemCallbackData *);
  static PetscErrorCode LocalJacobian(DALocalInfo *info, const U **x, Mat B, SNESProblemCallbackData *);

};

typedef SNESProblem<1,PetscScalar> SNESScalarProblem;
typedef SNESProblem<2,PISMVector2> SNESVectorProblem;



template<int DOF, class U>
PetscErrorCode SNESProblem<DOF,U>::LocalFunction(DALocalInfo *info,
                             const U **x, U **f,
                             SNESProblem<DOF,U>::SNESProblemCallbackData *cb)
{
  return cb->solver->compute_local_function(info,x,f);
}

template<int DOF, class U>
PetscErrorCode SNESProblem<DOF,U>::LocalJacobian(DALocalInfo *info,
                             const U **x, Mat J,
                             SNESProblem<DOF,U>::SNESProblemCallbackData *cb)
{
  return cb->solver->compute_local_jacobian(info,x,J);
}


template<int DOF, class U>
SNESProblem<DOF,U>::SNESProblem(IceGrid &g) :
m_grid(g)
{
  PetscErrorCode ierr;
  ierr = setFromOptions(); CHKERRABORT(m_grid.com,ierr);
  ierr = initialize();     CHKERRABORT(m_grid.com,ierr);
}
template<int DOF, class U>
SNESProblem<DOF,U>::~SNESProblem()
{
  PetscErrorCode ierr = finalize();
  CHKERRABORT(m_grid.com,ierr);  
}

template<int DOF, class U>
PetscErrorCode SNESProblem<DOF,U>::setFromOptions()
{
  return 0;
}

template<int DOF, class U>
PetscErrorCode SNESProblem<DOF,U>::initialize()
{
  PetscErrorCode ierr;

  // mimic IceGrid::createDA() with TRANSPOSE :
  PetscInt stencil_width=1;
  ierr = DACreate2d(m_grid.com, DA_XYPERIODIC, DA_STENCIL_BOX,
                    m_grid.My, m_grid.Mx,
                    m_grid.Ny, m_grid.Nx,
                    DOF, stencil_width,
                    m_grid.procs_y.data(), m_grid.procs_x.data(),
                    &m_DA); CHKERRQ(ierr);

  ierr = DACreateGlobalVector(m_DA, &m_X); CHKERRQ(ierr);
  ierr = DACreateGlobalVector(m_DA, &m_F); CHKERRQ(ierr);
  ierr = DAGetMatrix(m_DA, "baij",  &m_J); CHKERRQ(ierr);

  ierr = SNESCreate(m_grid.com, &m_snes);CHKERRQ(ierr);

  // Set the SNES callbacks to call into our compute_local_function and compute_local_jacobian
  // methods via SSAFEFunction and SSAFEJ
  ierr = DASetLocalFunction(m_DA,(DALocalFunction1)SNESProblem<DOF,U>::LocalFunction);CHKERRQ(ierr);
  ierr = DASetLocalJacobian(m_DA,(DALocalFunction1)SNESProblem<DOF,U>::LocalJacobian);CHKERRQ(ierr);
  m_callbackData.da = m_DA;  m_callbackData.solver = this;
  ierr = SNESSetFunction(m_snes, m_F, SNESDAFormFunction, &m_callbackData);CHKERRQ(ierr);
  ierr = SNESSetJacobian(m_snes, m_J, m_J, SNESDAComputeJacobian, &m_callbackData);CHKERRQ(ierr);  

  ierr = SNESSetFromOptions(m_snes);CHKERRQ(ierr);

  return 0;
}

//! Undo the allocations of SSAFEM::allocate_fem; called by the destructor.
template<int DOF, class U>
PetscErrorCode SNESProblem<DOF,U>::finalize() {
  PetscErrorCode ierr;

  ierr = SNESDestroy(m_snes);CHKERRQ(ierr);
  ierr = VecDestroy(m_X); CHKERRQ(ierr);
  ierr = VecDestroy(m_F); CHKERRQ(ierr);
  ierr = MatDestroy(m_J); CHKERRQ(ierr);

  return 0;
}

template<int DOF, class U>
const char *SNESProblem<DOF,U>::name()
{
  return "UnnamedProblem";
}

template<int DOF, class U>
PetscErrorCode SNESProblem<DOF,U>::solve()
{
  PetscErrorCode ierr;

  // Solve:
  ierr = SNESSolve(m_snes,NULL,m_X);CHKERRQ(ierr);
  
  // See if it worked.
  SNESConvergedReason reason;
  ierr = SNESGetConvergedReason( m_snes, &reason); CHKERRQ(ierr);
  if(reason < 0)
  {
    SETERRQ2(1, 
      "SNESProblem %s solve failed to converge (SNES reason %s)\n\n", name(), SNESConvergedReasons[reason]);
  }

  verbPrintf(1,m_grid.com,"SNESProblem %s converged (SNES reason %s)\n", name(), SNESConvergedReasons[reason]);
  return 0;
}


#endif
