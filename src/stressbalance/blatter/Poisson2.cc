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

#include <cassert>

#include "Poisson2.hh"
#include "pism/util/fem/FEM.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace stressbalance {

PetscErrorCode Poisson2::function_callback(DMDALocalInfo *info,
                                           const double **x, double **f,
                                           CallbackData *data) {
  try {
    data->solver->compute_residual(info, x, f);
  } catch (...) {
    MPI_Comm com = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)data->da, &com); CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

static double xy(double L, double delta, int k) {
  return -L + k * delta;
}

static bool dirichlet_node(const DMDALocalInfo *info, const fem::Element3::GlobalIndex& I) {
  return
    (I.i == 0 or I.i == info->mx - 1) or
    (I.j == 0 or I.j == info->my - 1);
}

static double u_bc(double x, double y) {
  return 2.0 * (1 + y) / ((3 + x) * (3 + x)  +  (1 + y) * (1 + y));
}

static double F(double x, double y) {
  (void) x;
  (void) y;
  return 0.0;
}

void Poisson2::compute_residual(DMDALocalInfo *info,
                                const double **x, double **f) {

  // Compute grid spacing from domain dimensions and the grid size
  double
    Lx = m_grid->Lx(),
    Ly = m_grid->Ly(),
    dx = 2.0 * Lx / (info->mx - 1),
    dy = 2.0 * Ly / (info->my - 1);

  fem::Q1Element2 E(*info, dx, dy, fem::Q1Quadrature4());

  // Compute the residual at Dirichlet BC nodes and reset the residual to zero elsewhere.
  //
  // Setting it to zero is necessary because we call DMDASNESSetFunctionLocal() with
  // INSERT_VALUES.
  //
  // here we loop over all the *owned* nodes
  for (int j = info->ys; j < info->ys + info->ym; j++) {
    for (int i = info->xs; i < info->xs + info->xm; i++) {
      if (dirichlet_node(info, {i, j, 0})) {
        double
          xx = xy(Lx, dx, i),
          yy = xy(Ly, dy, j);

        f[j][i] = u_bc(xx, yy) - x[j][i];
      } else {
        f[j][i] = 0.0;
      }
    }
  }

  // values at element nodes
  int Nk = E.n_chi();
  std::vector<double>
    x_nodal(Nk), y_nodal(Nk),
    R_nodal(Nk), u_nodal(Nk);

  // values at quadrature points
  int Nq = E.n_pts();
  std::vector<double> u(Nq), u_x(Nq), u_y(Nq);
  std::vector<double> xq(Nq), yq(Nq);

  // loop over all the elements that have at least one owned node
  for (int j = info->gys; j < info->gys + info->gym - 1; j++) {
    for (int i = info->gxs; i < info->gxs + info->gxm - 1; i++) {

      // reset element residual to zero in preparation
      for (int n = 0; n < Nk; ++n) {
        R_nodal[n] = 0.0;
      }

      E.reset(i, j);

      // compute coordinates of the nodes of this element
      for (int n = 0; n < Nk; ++n) {
        int ii, jj;
        E.local_to_global(n, ii, jj);

        x_nodal[n] = xy(Lx, dx, ii);
        y_nodal[n] = xy(Ly, dy, jj);
      }

      // get nodal values of u
      E.nodal_values(x, u_nodal.data());

      // take care of Dirichlet BC: don't contribute to Dirichlet nodes and set nodal
      // values of the current iterate to the BC value
      for (int n = 0; n < Nk; ++n) {
        int ii, jj;
        E.local_to_global(n, ii, jj);
        if (dirichlet_node(info, {ii, jj, 0})) {
          E.mark_row_invalid(n);
          u_nodal[n] = u_bc(x_nodal[n], y_nodal[n]);
        }
      }

      // evaluate u and its partial derivatives at quadrature points
      E.evaluate(u_nodal.data(), u.data(), u_x.data(), u_y.data());

      // coordinates of quadrature points
      E.evaluate(x_nodal.data(), xq.data());
      E.evaluate(y_nodal.data(), yq.data());

      // loop over all quadrature points
      for (int q = 0; q < Nq; ++q) {
        auto W = E.weight(q);

        // loop over all test functions
        for (int t = 0; t < Nk; ++t) {
          const auto &psi = E.chi(q, t);

          R_nodal[t] += W * (u_x[q] * psi.dx + u_y[q] * psi.dy
                             - F(xq[q], yq[q]) * psi.val);
        }
      }

      E.add_contribution(R_nodal.data(), f);
    } // end of the loop over i
  } // end of the loop over j
}

Poisson2::Poisson2(IceGrid::ConstPtr grid)
  : ShallowStressBalance(grid),
    m_solution(grid, "solution", WITHOUT_GHOSTS),
    m_exact(grid, "exact", WITHOUT_GHOSTS) {
  int ierr = 0;

  {
    auto pism_da = grid->get_dm(1, 0);

    PetscInt dim, Mx, My, Nx, Ny;
    PetscInt
      dof           = 1,
      stencil_width = 1;

    ierr = DMDAGetInfo(*pism_da,
                       &dim,
                       &Mx,
                       &My,
                       NULL, /* Mz */
                       &Nx,  /* number of processors in y-direction */
                       &Ny,  /* number of processors in x-direction */
                       NULL, /* ditto, z-direction */
                       NULL, /* number of degrees of freedom per node */
                       NULL, /* stencil width */
                       NULL, NULL, NULL, /* types of ghost nodes at the boundary */
                       NULL);            /* stencil width */
    PISM_CHK(ierr, "DMDAGetInfo");
    assert(dim == 2);

    const PetscInt *lx, *ly;

    ierr = DMDAGetOwnershipRanges(*pism_da, &lx, &ly, NULL);

    ierr = DMDACreate2d(PETSC_COMM_WORLD,
                        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                        DMDA_STENCIL_BOX,
                        Mx, My,
                        Nx, Ny,
                        dof,           // dof
                        stencil_width, // stencil width
                        lx, ly,
                        m_da.rawptr());
    PISM_CHK(ierr, "DMDACreate3d");

    ierr = DMSetFromOptions(m_da);
    PISM_CHK(ierr, "DMSetFromOptions");

    ierr = DMSetUp(m_da);
    PISM_CHK(ierr, "DMSetUp");
  }

  // SNES
  {
    ierr = SNESCreate(PETSC_COMM_WORLD, m_snes.rawptr());
    PISM_CHK(ierr, "SNESCreate");

    // ierr = SNESSetOptionsPrefix(m_snes, "poi_");
    // PISM_CHK(ierr, "SNESSetOptionsPrefix");

    ierr = SNESSetDM(m_snes, m_da);
    PISM_CHK(ierr, "SNESSetDM");

    m_callback_data.da = m_da;
    m_callback_data.solver = this;

    ierr = DMDASNESSetFunctionLocal(m_da, INSERT_VALUES,
                                    (DMDASNESFunction)function_callback,
                                    &m_callback_data);
    PISM_CHK(ierr, "DMDASNESSetFunctionLocal");

    ierr = SNESSetFromOptions(m_snes);
    PISM_CHK(ierr, "SNESSetFromOptions");
  }

  // set the initial guess
  m_solution.set(0.0);
}

Poisson2::~Poisson2() {
  // empty
}

void Poisson2::exact_solution(IceModelVec2S &result) {

  // Compute grid spacing from domain dimensions and the grid size
  double
    Lx = m_grid->Lx(),
    Ly = m_grid->Ly(),
    dx = 2.0 * Lx / (m_grid->Mx() - 1),
    dy = 2.0 * Ly / (m_grid->My() - 1);

  IceModelVec::AccessList list{&result};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double
      xx = xy(Lx, dx, i),
      yy = xy(Ly, dy, j);

    result(i, j) = u_bc(xx, yy);
  }
}

double Poisson2::error() const {
  IceModelVec2S difference(m_grid, "difference", WITHOUT_GHOSTS);

  m_exact.add(-1.0, m_solution, difference);

  return difference.norm(NORM_INFINITY);
}

void Poisson2::update(const Inputs &inputs, bool) {
  (void) inputs;

  int ierr = 0;

  ierr = SNESSolve(m_snes, NULL, m_solution.vec()); PISM_CHK(ierr, "SNESSolve");

  exact_solution(m_exact);
}

const IceModelVec2S& Poisson2::solution() const {
  return m_solution;
}

const IceModelVec2S& Poisson2::exact() const {
  return m_exact;
}

} // end of namespace stressbalance
} // end of namespace pism


