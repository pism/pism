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

#include "Poisson3.hh"

#include "pism/util/error_handling.hh"

namespace pism {
namespace stressbalance {

PetscErrorCode Poisson3::function_callback(DMDALocalInfo *info,
                                           const double ***x, double ***f,
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

void Poisson3::compute_local_function(DMDALocalInfo *info,
                                      const double ***x, double ***f) {
  const int
    xs = info->xs,
    ys = info->ys,
    zs = info->zs,
    xm = info->xm,
    ym = info->ym,
    zm = info->zm;

  for (int k = zs; k < zs + zm; k++) {
    for (int j = ys; j < ys + ym; j++) {
      for (int i = xs; i < xs + xm; i++) {
        f[k][j][i]  =  x[k][j][i];
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


