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

#include "Poisson3.hh"
#include "pism/util/fem/FEM.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace stressbalance {

PetscErrorCode Poisson3::function_callback(DMDALocalInfo *info,
                                           const double ***x, double ***f,
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

void Poisson3::compute_residual(DMDALocalInfo *info,
                                const double ***x, double ***f) {
  assert(info->sw == 1);

  // We start with ghost nodes in x and y directions to include contributions of elements
  // not owned by the current rank. The whole column (in the z direction) is owned by the
  // same rank, so this is not necessary there.
  const int
    gxs = info->gxs,
    xs  = info->xs,
    gys = info->gys,
    ys  = info->ys,
    gzs = info->gzs,
    zs  = info->zs,
    xm  = info->xm,
    ym  = info->ym,
    zm  = info->zm;

  fem::Q1Element3 E(*info, m_grid->dx(), m_grid->dy(), fem::Q13DQuadrature8());



  double
    Mz = info->mz,
    F  = 1.0,                   // right hand side
    b  = 0.0,                   // bottom elevation
    H  = 1000.0;                // thickness

  // loop over elements
  for (int k = gzs; k < zs + zm; k++) {
    for (int j = gys; j < ys + ym; j++) {
      for (int i = gxs; i < xs + xm; i++) {

        // Zero Dirichlet BC at the base
        if (k == 0) {
          f[k][j][i] = x[k][j][i];
          continue;
        }

        // reset residual to zero in preparation
        double R[8];
        for (int n = 0; n < 8; ++n) {
          R[n] = 0.0;
        }

        // compute z-coordinates for the nodes of this element
        std::vector<double> z(8, 0.0);
        {
          for (int n = 0; n < 4; ++n) {
            z[n] = b + H * k / (Mz - 1.0);
            z[n + 4] = H * (k + 1) / (Mz - 1.0);
          }
        }

        E.reset(i, j, k, z);

        // evaluate u, u_x, u_y, u_z at quadrature points
        double u_nodal[8], u[8], u_x[8], u_y[8], u_z[8];
        E.nodal_values(x, u_nodal);
        E.evaluate(u_nodal, u, u_x, u_y, u_z);

        // loop over all quadrature points
        for (int q = 0; q < E.n_pts(); ++q) {
          auto W = E.weight(q);

          // loop over all test functions
          for (int s = 0; s < E.n_chi(); ++s) {
            const auto &psi = E.chi(q, s);

            R[s] += W * (u_x[q] * psi.dx + u_y[q] * psi.dy + u_z[q] * psi.dz - F * psi.val);
          }
        }
        E.add_contribution(R, f);
      }
    }
  }
}

Poisson3::Poisson3(IceGrid::ConstPtr grid)
  : ShallowStressBalance(grid) {
  int ierr = 0;

  auto pism_da = grid->get_dm(1, 0);

  PetscInt dim, Mx, My, Nx, Ny;
  PetscInt Mz = 4, Nz = 1, dof = 1, stencil_width = 1;

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

  // DM
  {
    ierr = DMDACreate3d(PETSC_COMM_WORLD,
                        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                        DMDA_STENCIL_BOX,
                        Mx, My, Mz,
                        Nx, Ny, Nz,
                        dof,           // dof
                        stencil_width, // stencil width
                        lx, ly, NULL,
                        m_da.rawptr());
    PISM_CHK(ierr, "DMDACreate3d");

    ierr = DMSetFromOptions(m_da);
    PISM_CHK(ierr, "DMSetFromOptions");

    ierr = DMSetUp(m_da);
    PISM_CHK(ierr, "DMSetUp");
  }

  // Vecs, Mat
  {
    ierr = DMCreateGlobalVector(m_da, m_x.rawptr());
    PISM_CHK(ierr, "DMCreateGlobalVector");

    ierr = VecDuplicate(m_x, m_r.rawptr());
    PISM_CHK(ierr, "VecDuplicate");

    // ierr = DMCreateMatrix(m_da, m_J.rawptr());
    // PISM_CHK(ierr, "DMCreateMatrix");
  }

  // SNES
  {
    ierr = SNESCreate(PETSC_COMM_WORLD, m_snes.rawptr());
    PISM_CHK(ierr, "SNESCreate");

    ierr = SNESSetOptionsPrefix(m_snes, "poi_");
    PISM_CHK(ierr, "SNESSetOptionsPrefix");

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
  ierr = VecSet(m_x, 1.0);
  PISM_CHK(ierr, "VecSet");

  {
    std::vector<double> sigma(Mz);
    double dz = 1.0 / (Mz - 1);
    for (int i = 0; i < Mz; ++i) {
      sigma[i] = i * dz;
    }
    sigma.back() = 1.0;

    std::map<std::string,std::string> z_attrs =
      {{"axis", "Z"},
       {"long_name", "scaled Z-coordinate in the ice (z_base=0, z_surface=1)"},
       {"units", "1"},
       {"positive", "up"}};

    m_xx.reset(new IceModelVec3Custom(grid, "xx", "z_sigma", sigma, z_attrs));
    m_xx->set_attrs("diagnostic", "solution", "1", "1", "", 0);
  }
}

Poisson3::~Poisson3() {
  // empty
}

void Poisson3::update(const Inputs &inputs, bool) {
  (void) inputs;

  int ierr = 0;

  ierr = SNESSolve(m_snes, NULL, m_x); PISM_CHK(ierr, "SNESSolve");

  int Mz = m_xx->levels().size();

  {
    double ***x = nullptr;
    ierr = DMDAVecGetArray(m_da, m_x, &x); PISM_CHK(ierr, "DMDAVecGetArray");

    IceModelVec::AccessList list{m_xx.get()};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      auto c = m_xx->get_column(i, j);

      for (int k = 0; k < Mz; ++k) {
        c[k] = x[k][j][i];
      }
    }

    ierr = DMDAVecRestoreArray(m_da, m_x, &x); PISM_CHK(ierr, "DMDAVecRestoreArray");
  }
}

IceModelVec3Custom::Ptr Poisson3::x() const {
  return m_xx;
}


} // end of namespace stressbalance
} // end of namespace pism


