/* Copyright (C) 2019, 2020 PISM Authors
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

#include "Poisson.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Context.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"

namespace pism {

Poisson::Poisson(IceGrid::ConstPtr grid)
  : m_grid(grid),
    m_log(grid->ctx()->log()),
    m_b(grid, "poisson_rhs", WITHOUT_GHOSTS),
    m_x(grid, "poisson_x", WITHOUT_GHOSTS),
    m_mask(grid, "poisson_mask", WITH_GHOSTS){

  m_da = m_x.dm();

  // PETSc objects and settings
  {
    PetscErrorCode ierr;
    ierr = DMSetMatType(*m_da, MATAIJ);
    PISM_CHK(ierr, "DMSetMatType");

    ierr = DMCreateMatrix(*m_da, m_A.rawptr());
    PISM_CHK(ierr, "DMCreateMatrix");

    ierr = KSPCreate(m_grid->com, m_KSP.rawptr());
    PISM_CHK(ierr, "KSPCreate");

    ierr = KSPSetOptionsPrefix(m_KSP, "poisson_");
    PISM_CHK(ierr, "KSPSetOptionsPrefix");

    // Process options:
    ierr = KSPSetFromOptions(m_KSP);
    PISM_CHK(ierr, "KSPSetFromOptions");
  }
}

/*!
 * Solve the Poisson equation on the domain defined by `mask == 1` with Dirichlet BC
 * provided in `bc` (used only where `mask == 0`, possibly redundant away from the domain)
 * with the constant right hand side `rhs`.
 *
 * Set the mask to 2 to use zero Neumann BC.
 */
int Poisson::solve(const IceModelVec2Int& mask, const IceModelVec2S& bc, double rhs,
                   bool reuse_matrix) {

  PetscErrorCode ierr;

  // make a ghosted copy of the mask
  m_mask.copy_from(mask);

  if (reuse_matrix) {
    // Use non-zero initial guess. I assume that re-using the matrix means that the BC and
    // RHS provided are close to the ones used in the previous call and the solution
    // stored in m_x is a good initial guess for the current problem.
    ierr = KSPSetInitialGuessNonzero(m_KSP, PETSC_TRUE);
    PISM_CHK(ierr, "KSPSetInitialGuessNonzero");
  } else {
    ierr = KSPSetInitialGuessNonzero(m_KSP, PETSC_FALSE);
    PISM_CHK(ierr, "KSPSetInitialGuessNonzero");

    assemble_matrix(m_mask, m_A);
  }

  assemble_rhs(rhs, m_mask, bc, m_b);

  // Call PETSc to solve linear system by iterative method.
  ierr = KSPSetOperators(m_KSP, m_A, m_A);
  PISM_CHK(ierr, "KSPSetOperator");

  ierr = KSPSolve(m_KSP, m_b.vec(), m_x.vec());
  PISM_CHK(ierr, "KSPSolve");

  // Check if diverged
  KSPConvergedReason  reason;
  ierr = KSPGetConvergedReason(m_KSP, &reason);
  PISM_CHK(ierr, "KSPGetConvergedReason");

  if (reason < 0) {
    // KSP diverged
    m_log->message(1,
                   "PISM ERROR: KSP iteration failed while solving the Poisson equation\n"
                   "            reason = %d = '%s'\n",
                   reason, KSPConvergedReasons[reason]);

    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "KSP iteration failed: %s",
                                  KSPConvergedReasons[reason]);
  }

  // report on KSP success
  PetscInt ksp_iterations = 0;
  ierr = KSPGetIterationNumber(m_KSP, &ksp_iterations);
  PISM_CHK(ierr, "KSPGetIterationNumber");

  return ksp_iterations;
}

const IceModelVec2S& Poisson::solution() const {
  return m_x;
}

// Maxima code deriving the discretization
//
// /* Shift in x and y directions. */
// shift(expr, dx, dy) := op(expr)[args(expr)[1] + dx, args(expr)[2] + dy]$
//
// d_px(var) := (shift(var, 1, 0) - var)$
// d_mx(var) := (var - shift(var, -1, 0))$
//
// d_py(var) := (shift(var, 0, 1) - var)$
// d_my(var) := (var - shift(var, 0, -1))$
//
// constants : [C_x = 1 / (dx^2), C_y = 1 / (dy^2)]$
//
// /* discretization of -\nabla \dot (D \nabla W) */
// L: (E * d_px(W[i, j]) - W * d_mx(W[i, j])) / dx^2 +
//    (N * d_py(W[i, j]) - S * d_my(W[i, j])) / dy^2$
//
// /* discretization of the Poisson equation */
// eq: - L = f;
//
// /* cleaned up equation */
// s : map(lambda([x,y], solve(x, y)[1]), constants, [dx^2, dy^2]);
// eq2 : at(eq, s);
//
// /* compute matrix coefficients */
// for m: -1 thru 1 do (for n: -1 thru 1 do
//   (c[2 - n, m + 2] : factor(ratcoef(lhs(eq2), W[i + m, j + n]))))$
//
// /* print results */
// A : genmatrix(c, 3, 3);
//
// b : rhs(eq2);
//
// print(''out)$
void Poisson::assemble_matrix(const IceModelVec2Int &mask, Mat A) {
  PetscErrorCode ierr = 0;

  const double
    dx = m_grid->dx(),
    dy = m_grid->dy(),
    C_x = 1.0 / (dx * dx),
    C_y = 1.0 / (dy * dy);

  const int
    nrow = 1,
    ncol = 5,
    Mx   = m_grid->Mx(),
    My   = m_grid->My();

  ierr = MatZeroEntries(A); PISM_CHK(ierr, "MatZeroEntries");

  IceModelVec::AccessList list{&mask};

  /* matrix assembly loop */
  ParallelSection loop(m_grid->com);
  try {
    MatStencil row, col[ncol];
    row.c = 0;

    for (int m = 0; m < ncol; m++) {
      col[m].c = 0;
    }

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      /* i indices */
      const int I[] = {i, i - 1,  i,  i + 1, i};

      /* j indices */
      const int J[] = {j + 1, j,  j,  j, j - 1};

      row.i = i;
      row.j = j;

      for (int m = 0; m < ncol; m++) {
        col[m].i = I[m];
        col[m].j = J[m];
      }

      auto M = mask.star(i, j);

      if (M.ij == 1) {
        // Regular location: use coefficients of the discretization of the Laplacian

        // Use zero Neumann BC if a neighbor is marked as a Neumann boundary
        double
          N = M.n == 2 ? 0.0 : 1.0,
          E = M.e == 2 ? 0.0 : 1.0,
          W = M.w == 2 ? 0.0 : 1.0,
          S = M.s == 2 ? 0.0 : 1.0;

        // Use zero Neumann BC at edges of the computational domain
        {
          N = j == My - 1 ? 0.0 : N;
          E = i == Mx - 1 ? 0.0 : E;
          W = i == 0      ? 0.0 : W;
          S = j == 0      ? 0.0 : S;
        }

        // discretization of the Laplacian
        double L[ncol] = {- N * C_y,
                          - W * C_x, (W + E) * C_x + (N + S) * C_y, - E * C_x,
                          - S * C_y};

        ierr = MatSetValuesStencil(A, nrow, &row, ncol, col, L, INSERT_VALUES);
        PISM_CHK(ierr, "MatSetValuesStencil");
      } else {
        // Boundary or outside the domain: assemble trivial equations (1 on the diagonal,
        // 0 elsewhere)
        double D[ncol] = {0.0,
                          0.0, 1.0, 0.0,
                          0.0};

        ierr = MatSetValuesStencil(A, nrow, &row, ncol, col, D, INSERT_VALUES);
        PISM_CHK(ierr, "MatSetValuesStencil");
      }

    } // i,j-loop
  } catch (...) {
    loop.failed();
  }
  loop.check();

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); PISM_CHK(ierr, "MatAssemblyBegin");
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); PISM_CHK(ierr, "MatAssemblyif");

#if (Pism_DEBUG==1)
  ierr = MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  PISM_CHK(ierr, "MatSetOption");
#endif

}

void Poisson::assemble_rhs(double rhs,
                           const IceModelVec2Int &mask,
                           const IceModelVec2S &bc,
                           IceModelVec2S &b) {
  IceModelVec::AccessList list{&mask, &bc, &b};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask.as_int(i, j) == 1) {
      // inside the domain
      b(i, j) = rhs;
    } else {
      // at the boundary or outside the domain
      b(i, j) = bc(i, j);
    }
  }
}

} // end of namespace pism
