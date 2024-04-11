// Copyright (C) 2011, 2012, 2014, 2015, 2016, 2017, 2022, 2024 David Maxwell
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


#ifndef PISM_SNESPROBLEM_H
#define PISM_SNESPROBLEM_H

#include "pism/util/Grid.hh" // inline implementation in the header uses Grid
#include "pism/util/Vector2d.hh" // to get Vector2
#include "pism/util/petscwrappers/SNES.hh"
#include "pism/util/Logger.hh"

namespace pism {

template<int DOF, class U> class SNESProblem {
public:
  SNESProblem(std::shared_ptr<const Grid> g);

  virtual ~SNESProblem();

  virtual void solve();

  virtual const std::string& name();

  virtual Vec solution()
  {
    return m_X;
  }

protected:

  virtual void compute_local_function(DMDALocalInfo *info, const U **xg, U **yg) = 0;
  virtual void compute_local_jacobian(DMDALocalInfo *info, const U **x,  Mat B) = 0;

  std::shared_ptr<const Grid> m_grid;

  petsc::Vec m_X;
  petsc::SNES m_snes;
  petsc::DM::Ptr m_DA;

private:

  struct CallbackData {
    DM da;
    SNESProblem<DOF,U> *solver;
  };

  CallbackData m_callbackData;

  static PetscErrorCode function_callback(DMDALocalInfo *info, const U **x, U **f,
                                          CallbackData *);
  static PetscErrorCode jacobian_callback(DMDALocalInfo *info, const U **x, Mat B,
                                          CallbackData *);
};

typedef SNESProblem<1,double> SNESScalarProblem;
typedef SNESProblem<2,Vector2d> SNESVectorProblem;

template<int DOF, class U>
PetscErrorCode SNESProblem<DOF,U>::function_callback(DMDALocalInfo *info,
                                                     const U **x, U **f,
                                                     SNESProblem<DOF,U>::CallbackData *cb) {
  try {
    cb->solver->compute_local_function(info,x,f);
  } catch (...) {
    MPI_Comm com = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)cb->da, &com); CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

template<int DOF, class U>
PetscErrorCode SNESProblem<DOF,U>::jacobian_callback(DMDALocalInfo *info,
                                                     const U **x, Mat J,
                                                     SNESProblem<DOF,U>::CallbackData *cb) {
  try {
    cb->solver->compute_local_jacobian(info, x, J);
  } catch (...) {
    MPI_Comm com = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)cb->da, &com); CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

template<int DOF, class U>
SNESProblem<DOF, U>::SNESProblem(std::shared_ptr<const Grid> g)
  : m_grid(g) {

  PetscErrorCode ierr;

  int stencil_width=1;
  m_DA = m_grid->get_dm(DOF, stencil_width);

  ierr = DMCreateGlobalVector(*m_DA, m_X.rawptr());
  PISM_CHK(ierr, "DMCreateGlobalVector");

  ierr = SNESCreate(m_grid->com, m_snes.rawptr());
  PISM_CHK(ierr, "SNESCreate");

  // Set the SNES callbacks to call into our compute_local_function and compute_local_jacobian
  m_callbackData.da = *m_DA;
  m_callbackData.solver = this;

  ierr = DMDASNESSetFunctionLocal(*m_DA, INSERT_VALUES,
#if PETSC_VERSION_LT(3,21,0)
                                  (DMDASNESFunction)SNESProblem<DOF, U>::function_callback,
#else
                                  (DMDASNESFunctionFn*)SNESProblem<DOF, U>::function_callback,
#endif
                                  &m_callbackData);
  PISM_CHK(ierr, "DMDASNESSetFunctionLocal");

  ierr = DMDASNESSetJacobianLocal(*m_DA,
#if PETSC_VERSION_LT(3,21,0)
                                  (DMDASNESJacobian)SNESProblem<DOF, U>::jacobian_callback,
#else
                                  (DMDASNESJacobianFn*)SNESProblem<DOF, U>::jacobian_callback,
#endif
                                  &m_callbackData);
  PISM_CHK(ierr, "DMDASNESSetJacobianLocal");

  ierr = DMSetMatType(*m_DA, "baij");
  PISM_CHK(ierr, "DMSetMatType");

  ierr = DMSetApplicationContext(*m_DA, &m_callbackData);
  PISM_CHK(ierr, "DMSetApplicationContext");

  ierr = SNESSetDM(m_snes, *m_DA);
  PISM_CHK(ierr, "SNESSetDM");

  ierr = SNESSetFromOptions(m_snes);
  PISM_CHK(ierr, "SNESSetFromOptions");
}

template<int DOF, class U>
SNESProblem<DOF,U>::~SNESProblem() {
  // empty
}

template<int DOF, class U>
const std::string& SNESProblem<DOF,U>::name() {
  return "UnnamedProblem";
}

template<int DOF, class U>
void SNESProblem<DOF,U>::solve() {
  PetscErrorCode ierr;

  // Solve:
  ierr = SNESSolve(m_snes,NULL,m_X); PISM_CHK(ierr, "SNESSolve");

  // See if it worked.
  SNESConvergedReason reason;
  ierr = SNESGetConvergedReason(m_snes, &reason); PISM_CHK(ierr, "SNESGetConvergedReason");
  if (reason < 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "SNESProblem %s solve failed to converge (SNES reason %s)",
                                  name().c_str(), SNESConvergedReasons[reason]);
  }

  m_grid->ctx()->log()->message(1, "SNESProblem %s converged (SNES reason %s)\n",
                               name().c_str(), SNESConvergedReasons[reason]);
}

} // end of namespace pism

#endif // PISM_SNESPROBLEM_H
