/* Copyright (C) 2020 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>

#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/util/petscwrappers/SNES.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/petscwrappers/Mat.hh"
#include "pism/util/error_handling.hh"

/*
   User-defined application context - contains data needed by the
   application-provided call-back routines, FormJacobian() and
   FormFunction().
*/
typedef struct {
  PetscReal param; /* test problem parameter */
  DM da;           /* distributed array data structure */
} AppCtx;

/*
   User-defined routines
*/
extern PetscErrorCode FormFunctionLocal(SNES, Vec, Vec, void *);
extern PetscErrorCode FormFunction(SNES, Vec, Vec, void *);
extern PetscErrorCode FormInitialGuess(AppCtx *, Vec);
extern PetscErrorCode FormJacobian(SNES, Vec, Mat, Mat, void *);

namespace pism {
namespace stressbalance {

class Poisson : public ShallowStressBalance {
public:
  Poisson(IceGrid::ConstPtr grid);
  virtual ~Poisson();

  void update(const Inputs &inputs, bool);
protected:

  petsc::DM m_da;
  petsc::Vec m_x, m_r;
  petsc::SNES m_snes;
  petsc::Mat m_J;

  MatFDColoring m_coloring;

  struct CallbackData {
    DM da;
    Poisson *solver;
  };

  CallbackData m_callback_data;

  void compute_local_function(DMDALocalInfo *info, const double **xg, double **yg);
  static PetscErrorCode function_callback(DMDALocalInfo *info, const double **x, double **f,
                                          CallbackData *);
};

PetscErrorCode Poisson::function_callback(DMDALocalInfo *info, const double **x, double **f,
                                          CallbackData *data) {
  try {
    data->solver->compute_local_function(info, x, f);
  } catch (...) {
    MPI_Comm com = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)data->da, &com); CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

void Poisson::compute_local_function(DMDALocalInfo *info, const double **xg, double **yg) {
  (void) info;
  (void) xg;
  (void) yg;
}


Poisson::Poisson(IceGrid::ConstPtr grid)
  : ShallowStressBalance(grid) {
  int ierr = 0;

  bool
    coloring       = false,
    coloring_ds    = false,
    local_coloring = false;

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESCreate(PETSC_COMM_WORLD, m_snes.rawptr());
  PISM_CHK(ierr, "SNESCreate");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create distributed array (DMDA) to manage parallel grid and vectors
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMDACreate3d(PETSC_COMM_WORLD,
                      DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                      DMDA_STENCIL_STAR,
                      4, 4, 4,
                      PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                      1,        // dof
                      1,        // stencil width
                      NULL, NULL, NULL,
                      m_da.rawptr());
  PISM_CHK(ierr, "DMDACreate3d");

  ierr = DMSetFromOptions(m_da);
  PISM_CHK(ierr, "DMSetFromOptions");

  ierr = DMSetUp(m_da);
  PISM_CHK(ierr, "DMSetUp");

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Extract global vectors from DMDA; then duplicate for remaining
      vectors that are the same types
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMCreateGlobalVector(m_da, m_x.rawptr());
  PISM_CHK(ierr, "DMCreateGlobalVector");

  ierr = VecDuplicate(m_x, m_r.rawptr());
  PISM_CHK(ierr, "VecDuplicate");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Set function evaluation routine and vector
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = DMDASNESSetFunctionLocal(m_da, INSERT_VALUES,
                                  (DMDASNESFunction)function_callback,
                                  &m_callback_data);
  PISM_CHK(ierr, "DMDASNESSetFunctionLocal");

  if (true) {
    ierr = DMSetMatType(m_da, MATAIJ);
    PISM_CHK(ierr, "DMSetMatType");

    ierr = DMCreateMatrix(m_da, m_J.rawptr());
    PISM_CHK(ierr, "DMCreateMatrix");

    if (coloring) {
      ISColoring iscoloring;
      if (!local_coloring) {
        ierr = DMCreateColoring(m_da, IS_COLORING_GLOBAL, &iscoloring);
        PISM_CHK(ierr, "DMCreateColoring");

        ierr = MatFDColoringCreate(m_J, iscoloring, &m_coloring);
        PISM_CHK(ierr, "MatFDColoringCreate");

        ierr = MatFDColoringSetFunction(m_coloring, (PetscErrorCode(*)(void))FormFunction,
                                        &m_callback_data);
        PISM_CHK(ierr, "MatFDColoringSetFunction");
      } else {
        ierr = DMCreateColoring(m_da, IS_COLORING_LOCAL, &iscoloring);
        PISM_CHK(ierr, "DMCreateColoring");

        ierr = MatFDColoringCreate(m_J, iscoloring, &m_coloring);
        PISM_CHK(ierr, "MatFDColoringCreate");

        ierr = MatFDColoringUseDM(m_J, m_coloring);
        PISM_CHK(ierr, "MatFDColoringUseDM");

        ierr = MatFDColoringSetFunction(m_coloring,
                                        (PetscErrorCode(*)(void))function_callback,
                                        &m_callback_data);
        PISM_CHK(ierr, "MatFDColoringSetFunction");
      }
      if (coloring_ds) {
        ierr = MatFDColoringSetType(m_coloring, MATMFFD_DS);
        PISM_CHK(ierr, "MatFDColoringSetType");
      }
      ierr = MatFDColoringSetFromOptions(m_coloring);
      PISM_CHK(ierr, "MatFDColoringSetFromOptions");

      ierr = MatFDColoringSetUp(m_J, iscoloring, m_coloring);
      PISM_CHK(ierr, "MatFDColoringSetUp");

      ierr = SNESSetJacobian(m_snes, m_J, m_J, SNESComputeJacobianDefaultColor, m_coloring);
      PISM_CHK(ierr, "SNESSetJacobian");

      ierr = ISColoringDestroy(&iscoloring);
      PISM_CHK(ierr, "ISColoringDestroy");
    } else {
      ierr = SNESSetJacobian(m_snes, m_J, m_J, FormJacobian, &m_callback_data);
      PISM_CHK(ierr, "SNESSetJacobian");
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Customize nonlinear solver; set runtime options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESSetDM(m_snes, m_da);
  PISM_CHK(ierr, "SNESSetDM");

  ierr = SNESSetFromOptions(m_snes);
  PISM_CHK(ierr, "SNESSetFromOptions");

  // set the initial guess
}

Poisson::~Poisson() {
  MatFDColoringDestroy(&m_coloring);
}

void Poisson::update(const Inputs &inputs, bool) {
  (void) inputs;

  int ierr = 0;
  PetscInt its = 0;
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Solve nonlinear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  ierr = SNESSolve(m_snes, NULL, m_x); PISM_CHK(ierr, "SNESSolve");
  ierr = SNESGetIterationNumber(m_snes, &its); PISM_CHK(ierr, "SNESGetIterationNumber");

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Explicitly check norm of the residual of the solution
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  // ierr = FormFunction(snes, x, r, (void *)&user); CHKERRQ(ierr);
  // ierr = VecNorm(r, NORM_2, &fnorm); CHKERRQ(ierr);
  // ierr = PetscPrintf(PETSC_COMM_WORLD,
  //                    "Number of SNES iterations = %D fnorm %g\n",
  //                    its, (double)fnorm); CHKERRQ(ierr);
}

} // end of namespace stressbalance
} // end of namespace pism

PetscErrorCode FormFunctionLocal(SNES snes, Vec localX, Vec F, void *ptr) {
  AppCtx *user = (AppCtx *)ptr;
  PetscErrorCode ierr;
  PetscInt i, j, k, Mx, My, Mz, xs, ys, zs, xm, ym, zm;
  PetscReal two = 2.0, lambda, hx, hy, hz, hxhzdhy, hyhzdhx, hxhydhz, sc;
  PetscScalar u_north, u_south, u_east, u_west, u_up, u_down, u, u_xx, u_yy, u_zz, ***x, ***f;
  DM da;

  PetscFunctionBeginUser;
  ierr = SNESGetDM(snes, &da); CHKERRQ(ierr);
  ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
                     PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE); CHKERRQ(ierr);

  lambda  = user->param;
  hx      = 1.0 / (PetscReal)(Mx - 1);
  hy      = 1.0 / (PetscReal)(My - 1);
  hz      = 1.0 / (PetscReal)(Mz - 1);
  sc      = hx * hy * hz * lambda;
  hxhzdhy = hx * hz / hy;
  hyhzdhx = hy * hz / hx;
  hxhydhz = hx * hy / hz;

  /*
     Get pointers to vector data
  */
  ierr = DMDAVecGetArrayRead(da, localX, &x); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, F, &f); CHKERRQ(ierr);

  /*
     Get local grid boundaries
  */
  ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);

  /*
     Compute function over the locally owned part of the grid
  */
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        if (i == 0 || j == 0 || k == 0 || i == Mx - 1 || j == My - 1 || k == Mz - 1) {
          f[k][j][i] = x[k][j][i];
        } else {
          u          = x[k][j][i];
          u_east     = x[k][j][i + 1];
          u_west     = x[k][j][i - 1];
          u_north    = x[k][j + 1][i];
          u_south    = x[k][j - 1][i];
          u_up       = x[k + 1][j][i];
          u_down     = x[k - 1][j][i];
          u_xx       = (-u_east + two * u - u_west) * hyhzdhx;
          u_yy       = (-u_north + two * u - u_south) * hxhzdhy;
          u_zz       = (-u_up + two * u - u_down) * hxhydhz;
          f[k][j][i] = u_xx + u_yy + u_zz - sc * PetscExpScalar(u);
        }
      }
    }
  }

  /*
     Restore vectors
  */
  ierr = DMDAVecRestoreArrayRead(da, localX, &x); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da, F, &f); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode FormFunction(SNES snes, Vec X, Vec F, void *ptr) {
  PetscErrorCode ierr;
  Vec localX;
  DM da;

  PetscFunctionBeginUser;
  ierr = SNESGetDM(snes, &da); CHKERRQ(ierr);
  ierr = DMGetLocalVector(da, &localX); CHKERRQ(ierr);

  /*
     Scatter ghost points to local vector,using the 2-step process
        DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.
  */
  ierr = DMGlobalToLocalBegin(da, X, INSERT_VALUES, localX); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da, X, INSERT_VALUES, localX); CHKERRQ(ierr);

  ierr = FormFunctionLocal(snes, localX, F, ptr); CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da, &localX); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode FormJacobian(SNES snes, Vec X, Mat J, Mat jac, void *ptr) {
  AppCtx *user = (AppCtx *)ptr; /* user-defined application context */
  Vec localX;
  PetscErrorCode ierr;
  PetscInt i, j, k, Mx, My, Mz;
  MatStencil col[7], row;
  PetscInt xs, ys, zs, xm, ym, zm;
  PetscScalar lambda, v[7], hx, hy, hz, hxhzdhy, hyhzdhx, hxhydhz, sc, ***x;
  DM da;

  PetscFunctionBeginUser;
  ierr = SNESGetDM(snes, &da); CHKERRQ(ierr);
  ierr = DMGetLocalVector(da, &localX); CHKERRQ(ierr);
  ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, &My, &Mz, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
                     PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE); CHKERRQ(ierr);

  lambda  = user->param;
  hx      = 1.0 / (PetscReal)(Mx - 1);
  hy      = 1.0 / (PetscReal)(My - 1);
  hz      = 1.0 / (PetscReal)(Mz - 1);
  sc      = hx * hy * hz * lambda;
  hxhzdhy = hx * hz / hy;
  hyhzdhx = hy * hz / hx;
  hxhydhz = hx * hy / hz;

  /*
     Scatter ghost points to local vector, using the 2-step process
        DMGlobalToLocalBegin(), DMGlobalToLocalEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.
  */
  ierr = DMGlobalToLocalBegin(da, X, INSERT_VALUES, localX); CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da, X, INSERT_VALUES, localX); CHKERRQ(ierr);

  /*
     Get pointer to vector data
  */
  ierr = DMDAVecGetArrayRead(da, localX, &x); CHKERRQ(ierr);

  /*
     Get local grid boundaries
  */
  ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm); CHKERRQ(ierr);

  /*
     Compute entries for the locally owned part of the Jacobian.
      - Currently, all PETSc parallel matrix formats are partitioned by
        contiguous chunks of rows across the processors.
      - Each processor needs to insert only elements that it owns
        locally (but any non-local elements will be sent to the
        appropriate processor during matrix assembly).
      - Here, we set all entries for a particular row at once.
      - We can set matrix entries either using either
        MatSetValuesLocal() or MatSetValues(), as discussed above.
  */
  for (k = zs; k < zs + zm; k++) {
    for (j = ys; j < ys + ym; j++) {
      for (i = xs; i < xs + xm; i++) {
        row.k = k;
        row.j = j;
        row.i = i;
        /* boundary points */
        if (i == 0 || j == 0 || k == 0 || i == Mx - 1 || j == My - 1 || k == Mz - 1) {
          v[0] = 1.0;
          ierr = MatSetValuesStencil(jac, 1, &row, 1, &row, v, INSERT_VALUES); CHKERRQ(ierr);
        } else {
          /* interior grid points */
          v[0]     = -hxhydhz;
          col[0].k = k - 1;
          col[0].j = j;
          col[0].i = i;
          v[1]     = -hxhzdhy;
          col[1].k = k;
          col[1].j = j - 1;
          col[1].i = i;
          v[2]     = -hyhzdhx;
          col[2].k = k;
          col[2].j = j;
          col[2].i = i - 1;
          v[3]     = 2.0 * (hyhzdhx + hxhzdhy + hxhydhz) - sc * PetscExpScalar(x[k][j][i]);
          col[3].k = row.k;
          col[3].j = row.j;
          col[3].i = row.i;
          v[4]     = -hyhzdhx;
          col[4].k = k;
          col[4].j = j;
          col[4].i = i + 1;
          v[5]     = -hxhzdhy;
          col[5].k = k;
          col[5].j = j + 1;
          col[5].i = i;
          v[6]     = -hxhydhz;
          col[6].k = k + 1;
          col[6].j = j;
          col[6].i = i;
          ierr     = MatSetValuesStencil(jac, 1, &row, 7, col, v, INSERT_VALUES); CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = DMDAVecRestoreArrayRead(da, localX, &x); CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(da, &localX); CHKERRQ(ierr);

  /*
     Assemble matrix, using the 2-step process:
       MatAssemblyBegin(), MatAssemblyEnd().
  */
  ierr = MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  /*
     Normally since the matrix has already been assembled above; this
     would do nothing. But in the matrix free mode -snes_mf_operator
     this tells the "matrix-free" matrix that a new linear system solve
     is about to be done.
  */
  ierr = MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  /*
     Tell the matrix we will never add a new nonzero location to the
     matrix. If we do, it will generate an error.
  */
  ierr = MatSetOption(jac, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);

  return 0;
}
