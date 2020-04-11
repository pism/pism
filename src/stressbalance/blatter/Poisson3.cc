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
#include <cmath>                // std::pow, std::fabs

using std::pow;
using std::fabs;

#include "Poisson3.hh"
#include "pism/util/fem/FEM.hh"
#include "pism/util/error_handling.hh"

namespace pism {
namespace stressbalance {

/*!
 * dot product
 */
static double dot(const std::vector<double> &a, const std::vector<double> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/*!
 * x and y coordinates
 *
 * @param[in] L domain half-width
 * @param[in] delta grid spacing
 * @param[in] k node index
 */
static double xy(double L, double delta, int k) {
  return -L + k * delta;
}

/*!
 * z coordinate
 *
 * @param[in] b surface elevation of the bottom of the domain
 * @param[in] H domain thickness
 * @param[in] Mz number of grid points in each vertical column
 * @param[in] k node index in the z direction
 */
static double z(double b, double H, int Mz, int k) {
  return b + H * k / (Mz - 1.0);
}

/*!
 * Returns true if a node is in the Dirichlet part of the boundary, false otherwise.
 */
static bool dirichlet_node(const DMDALocalInfo *info, const fem::Element3::GlobalIndex& I) {
  return (I.i == info->mx - 1) or (I.j == info->my - 1) or (I.k == info->mz - 1);
}

/*!
 * Returns true if a node is in the Neumann part of the boundary, false otherwise.
 */
static bool neumann_node(const DMDALocalInfo *info, const fem::Element3::GlobalIndex& I) {
  (void) info;
  return I.i == 0 or I.j == 0 or I.k == 0;
}

/*! Dirichlet BC and the exact solution

 b : -1 + x + y;
 n_b : [-diff(b, x), -diff(b, y), 1];
 u : x*y*(z+1)^2+(2.0*(y+1))/((y+1)^2+(x+2)^2)$
 grind('u = u);
 grind(F = ratsimp(-(diff(u, x, 2) + diff(u, y, 2) + diff(u, z, 2))));
 grind(u_x = diff(u, x));
 grind(u_y = diff(u, y));
 grind(u_z = diff(u, z));
*/
static double u_exact(double x, double y, double z) {
  return x * y * pow(z + 1, 2.0) + (2.0 * (y + 1)) / (pow(y + 1, 2.0) + pow(x + 2, 2.0));
}

/*!
 * Right hand side
 *
 * F = - (diff(u, x, 2) + diff(u, y, 2) + diff(u, z, 2))
 */
static double F(double x, double y, double z) {
  (void) z;
  return -2.0 * x * y;
}

/*!
 * Neumann BC
 */
static double G(double x, double y, double z, double b) {

  double u_x = (y * pow(z + 1, 2.0) -
                (4.0 * (x + 2) * (y + 1)) / pow(pow(y + 1, 2.0) + pow(x + 2, 2.0), 2.0));
  double u_y = x * pow(z + 1, 2.0) + 2.0 / (pow(y + 1, 2.0) + pow(x + 2, 2.0)) -
    (4.0 * pow(y + 1, 2.0)) / pow(pow(y + 1, 2.0) + pow(x + 2, 2.0), 2.0);
  double u_z = 2 * x * y * (z + 1);

  double eps = 1e-12;
  if (fabs(x - (-1.0)) < eps) {
    return u_x;
  } else if (fabs(y - (-1.0)) < eps) {
    return u_y;
  } else if (fabs(z - b) < eps) {
    // normal to the bottom surface {-b_x, -b_y, 1}
    std::vector<double> n = {-1, -1, 1}; // magnitude: sqrt(3)

    return dot({u_x, u_y, u_z}, n) / sqrt(3); // normalize
  } else {
    // We are not on a Neumann boundary. This value will not be used.
    return 0.0;
  }
}

void Poisson3::compute_residual(DMDALocalInfo *info,
                                const double ***x, double ***R) {
  // Stencil width of 1 is not very important, but if info->sw > 1 will lead to more
  // redundant computation (we would be looping over elements that don't contribute to any
  // owned nodes).
  assert(info->sw == 1);

  // Compute grid spacing from domain dimensions and the grid size
  double
    Lx = m_grid->Lx(),
    Ly = m_grid->Ly(),
    dx = 2.0 * Lx / (info->mx - 1),
    dy = 2.0 * Ly / (info->my - 1);

  fem::Q1Element3 E(*info, dx, dy, fem::Q13DQuadrature8());
  fem::Q1Element3Face E_face(dx, dy, fem::Q1Quadrature4());

  Parameters **P;
  Vec X;

  begin_2d_access(info->da, true, &X, &P);

  // Compute the residual at Dirichlet BC nodes and reset the residual to zero elsewhere.
  //
  // Setting it to zero is necessary because we call DMDASNESSetFunctionLocal() with
  // INSERT_VALUES.
  //
  // here we loop over all the *owned* nodes
  for (int k = info->zs; k < info->zs + info->zm; k++) {
    for (int j = info->ys; j < info->ys + info->ym; j++) {
      for (int i = info->xs; i < info->xs + info->xm; i++) {
        if (dirichlet_node(info, {i, j, k})) {
          double
            xx = xy(Lx, dx, i),
            yy = xy(Ly, dy, j),
            b  = P[j][i].bed,
            H  = P[j][i].thickness,
            zz = z(b, H, info->mz, k);

          // FIXME: scaling goes here
          R[k][j][i] = u_exact(xx, yy, zz) - x[k][j][i];
        } else {
          R[k][j][i] = 0.0;
        }
      }
    }
  }

  // values at element nodes
  int Nk = E.n_chi();

  // hard-wired maximum number of nodes per element
  const int Nk_max = 8;
  assert(Nk <= Nk_max);
  double
    x_nodal[Nk_max], y_nodal[Nk_max],
    R_nodal[Nk_max], u_nodal[Nk_max], b_nodal[Nk_max];
  std::vector<double> z_nodal(Nk);

  // values at quadrature points
  int Nq = E.n_pts();
  // hard-wired maximum number of quadrature points per element
  const int Nq_max = 16;
  assert(Nq <= Nq_max);
  double u[Nq_max], u_x[Nq_max], u_y[Nq_max], u_z[Nq_max];
  double xq[Nq_max], yq[Nq_max], zq[Nq_max], bq[Nq_max];

  // make sure that xq, yq, zq and big enough for quadrature points on element faces
  assert(E_face.n_pts() <= Nq_max);

  // loop over all the elements that have at least one owned node
  for (int k = info->gzs; k < info->gzs + info->gzm - 1; k++) {
    for (int j = info->gys; j < info->gys + info->gym - 1; j++) {
      for (int i = info->gxs; i < info->gxs + info->gxm - 1; i++) {

        // Reset element residual to zero in preparation.
        for (int n = 0; n < Nk; ++n) {
          R_nodal[n] = 0.0;
        }

        // Compute coordinates of the nodes of this element.
        for (int n = 0; n < Nk; ++n) {
          auto I = E.local_to_global(i, j, k, n);

          double H = P[I.j][I.i].thickness;

          b_nodal[n] = P[I.j][I.i].bed;

          x_nodal[n] = xy(Lx, dx, I.i);
          y_nodal[n] = xy(Ly, dy, I.j);
          z_nodal[n] = z(b_nodal[n], H, info->mz, I.k);
        }

        E.reset(i, j, k, z_nodal);

        // Get nodal values of u.
        E.nodal_values(x, u_nodal);

        // Take care of Dirichlet BC: don't contribute to Dirichlet nodes and set nodal
        // values of the current iterate to Dirichler BC values.
        for (int n = 0; n < Nk; ++n) {
          auto I = E.local_to_global(n);
          if (dirichlet_node(info, I)) {
            E.mark_row_invalid(n);
            u_nodal[n] = u_exact(x_nodal[n], y_nodal[n], z_nodal[n]);
          }
        }

        // evaluate u and its partial derivatives at quadrature points
        E.evaluate(u_nodal, u, u_x, u_y, u_z);

        // coordinates of quadrature points
        E.evaluate(x_nodal, xq);
        E.evaluate(y_nodal, yq);
        E.evaluate(z_nodal.data(), zq);

        // loop over all quadrature points
        for (int q = 0; q < Nq; ++q) {
          auto W = E.weight(q);

          // loop over all test functions
          for (int t = 0; t < Nk; ++t) {
            const auto &psi = E.chi(q, t);

            // FIXME: scaling goes here
            R_nodal[t] += W * (u_x[q] * psi.dx + u_y[q] * psi.dy + u_z[q] * psi.dz
                               - F(xq[q], yq[q], zq[q]) * psi.val);
          }
        }

        // loop over all faces
        for (int face = 0; face < fem::q13d::n_faces; ++face) {
          auto nodes = fem::q13d::incident_nodes[face];
          // loop over all nodes corresponding to this face. A face is a part of the
          // Neumann boundary if all four nodes are Neumann nodes. If a node is *both* a
          // Neumann and a Dirichlet node (this may happen), then we treat it as a Neumann
          // node here: add_contribution() will do the right thing later.
          bool neumann = true;
          for (int n = 0; n < 4; ++n) {
            auto I = E.local_to_global(nodes[n]);
            if (not neumann_node(info, I)) {
              neumann = false;
              break;
            }
          }

          if (not neumann) {
            continue;
          }

          E_face.reset(face, z_nodal);

          // compute physical coordinates of quadrature points on this face
          E_face.evaluate(x_nodal, xq);
          E_face.evaluate(y_nodal, yq);
          E_face.evaluate(z_nodal.data(), zq);

          // Compute the bed elevation at (below) quadrature points. This is needed to
          // compute G() below
          E_face.evaluate(b_nodal, bq);

          // loop over all quadrature points
          for (int q = 0; q < E_face.n_pts(); ++q) {
            auto W = E_face.weight(q);

            // loop over all test functions
            for (int t = 0; t < Nk; ++t) {
              auto psi = E_face.chi(q, t);

              // FIXME: scaling goes here
              R_nodal[t] += W * psi * G(xq[q], yq[q], zq[q], bq[q]);
            }
          }
        } // end of the loop over element faces

        E.add_contribution(R_nodal, R);
      } // end of the loop over i
    } // end of the loop over j
  } // end of the loop over k

  end_2d_access(info->da, true, &X, &P);
}

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

Poisson3::Poisson3(IceGrid::ConstPtr grid, int Mz)
  : ShallowStressBalance(grid), m_Mz(Mz) {
  int ierr = 0;

  // DM
  {
    auto pism_da = grid->get_dm(1, 0);

    PetscInt dim, Mx, My, Nx, Ny;
    PetscInt
      Nz            = 1,
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

    ierr = DMDACreate3d(PETSC_COMM_WORLD,
                        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                        DMDA_STENCIL_BOX,
                        Mx, My, m_Mz,
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

    setup_level(m_da);
  }

  // Vecs, Mat
  {
    ierr = DMCreateGlobalVector(m_da, m_x.rawptr());
    PISM_CHK(ierr, "DMCreateGlobalVector");

    // ierr = DMCreateMatrix(m_da, m_J.rawptr());
    // PISM_CHK(ierr, "DMCreateMatrix");
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
  ierr = VecSet(m_x, 0.0);
  PISM_CHK(ierr, "VecSet");

  {
    std::vector<double> sigma(m_Mz);
    double dz = 1.0 / (m_Mz - 1);
    for (int i = 0; i < m_Mz; ++i) {
      sigma[i] = i * dz;
    }
    sigma.back() = 1.0;

    std::map<std::string,std::string> z_attrs =
      {{"axis", "Z"},
       {"long_name", "scaled Z-coordinate in the ice (z_base=0, z_surface=1)"},
       {"units", "1"},
       {"positive", "up"}};

    m_solution.reset(new IceModelVec3Custom(grid, "solution", "z_sigma", sigma, z_attrs));
    m_solution->set_attrs("diagnostic", "solution", "1", "1", "", 0);

    m_exact.reset(new IceModelVec3Custom(grid, "exact", "z_sigma", sigma, z_attrs));
    m_exact->set_attrs("diagnostic", "exact", "1", "1", "", 0);
  }
}

void Poisson3::setup_level(DM dm)
{
  int ierr;

  MPI_Comm comm;
  ierr = PetscObjectGetComm((PetscObject)dm, &comm);
  PISM_CHK(ierr, "PetscObjectGetComm");

  // Get grid information
  PetscInt Mx, My, Mz, mx, my, mz, stencil_width;
  DMDAStencilType stencil_type;
  const PetscInt *lx, *ly, *lz;
  {
    ierr = DMDAGetInfo(dm,
                       NULL,       // dimensions
                       &Mx, &My, &Mz, // grid size
                       &mx, &my, &mz, // number of processors in each direction
                       NULL,           // number of degrees of freedom
                       &stencil_width,
                       NULL, NULL, NULL, // types of ghost nodes at the boundary
                       &stencil_type);
    PISM_CHK(ierr, "DMDAGetInfo");

    ierr = DMDAGetOwnershipRanges(dm, &lx, &ly, &lz);
    PISM_CHK(ierr, "DMDAGetOwnershipRanges");
  }

  // Create a 2D DMDA and a global Vec, then stash them in the dm passed to this method.
  {
    // compute the number of parameters per map-plane location
    int dof = sizeof(Parameters)/sizeof(double);

    DM  da;
    Vec parameters;

    ierr = DMDACreate2d(comm,
                        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                        stencil_type,
                        Mx, My,
                        mx, my,
                        dof,
                        stencil_width,
                        lx, ly,
                        &da);
    PISM_CHK(ierr, "DMDACreate2d");

    ierr = DMSetUp(da);
    PISM_CHK(ierr, "DMSetUp");

    ierr = DMCreateGlobalVector(da, &parameters);
    PISM_CHK(ierr, "DMCreateGlobalVector");

    ierr = PetscObjectCompose((PetscObject)dm, "2D_DM", (PetscObject)da);
    PISM_CHK(ierr, "PetscObjectCompose");
    ierr = PetscObjectCompose((PetscObject)dm, "2D_Vec", (PetscObject)parameters);
    PISM_CHK(ierr, "PetscObjectCompose");

    ierr = DMDestroy(&da);
    PISM_CHK(ierr, "DMDestroy");

    ierr = VecDestroy(&parameters);
    PISM_CHK(ierr, "VecDestroy");
  }

  // Create a 3D DMDA and a global Vec, then stash them in the dm passed to this method.
  {
    DM  da;
    Vec parameters;
    int dof = 1;

    ierr = DMDACreate3d(comm,
                        DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                        stencil_type,
                        Mx, My, Mz,
                        mx, my, mz,
                        dof,
                        stencil_width,
                        lx, ly, lz,
                        &da);
    PISM_CHK(ierr, "DMDACreate3d");

    ierr = DMSetUp(da);
    PISM_CHK(ierr, "DMSetUp");

    ierr = DMCreateGlobalVector(da, &parameters);
    PISM_CHK(ierr, "DMCreateGlobalVector");

    ierr = PetscObjectCompose((PetscObject)dm, "3D_DM", (PetscObject)da);
    PISM_CHK(ierr, "PetscObjectCompose");
    ierr = PetscObjectCompose((PetscObject)dm, "3D_Vec", (PetscObject)parameters);
    PISM_CHK(ierr, "PetscObjectCompose");

    ierr = DMDestroy(&da);
    PISM_CHK(ierr, "DMDestroy");

    ierr = VecDestroy(&parameters);
    PISM_CHK(ierr, "VecDestroy");
  }

  // get refinement level
  PetscInt level = 0;
  {
    PetscInt refinelevel, coarsenlevel;
    ierr = DMGetRefineLevel(dm, &refinelevel);
    PISM_CHK(ierr, "DMGetRefineLevel");
    ierr = DMGetCoarsenLevel(dm, &coarsenlevel);
    PISM_CHK(ierr, "DMGetCoarsenLevel");
    level = refinelevel - coarsenlevel;
  }

  // report
  {
    double
      Lx = 2.0 * m_grid->Lx(),
      Ly = 2.0 * m_grid->Ly();
    PetscPrintf(comm,
                "Level %D domain size (m) %8.2g x %8.2g, num elements %3d x %3d x %3d (%8d), size (m) %g x %g\n",
                level, Lx, Ly, Mx, My, Mz, Mx*My*Mz, Lx / (Mx - 1), Ly / (My - 1));
  }
}

void Poisson3::begin_2d_access(DM da, bool local, Vec *X_out, Parameters ***prm) {
  int ierr;

  DM  da_2d;
  Vec X;

  ierr = PetscObjectQuery((PetscObject)da, "2D_DM", (PetscObject*)&da_2d);
  PISM_CHK(ierr, "PetscObjectQuery");

  if (!da_2d) {
    throw RuntimeError(PISM_ERROR_LOCATION, "No 2D_DM composed with given DMDA");
  }

  ierr = PetscObjectQuery((PetscObject)da, "2D_Vec", (PetscObject*)&X);
  PISM_CHK(ierr, "PetscObjectQuery");

  if (!X) {
    throw RuntimeError(PISM_ERROR_LOCATION, "No DMDA_2D_Vec composed with given DMDA");
  }

  if (local) {
    ierr = DMGetLocalVector(da_2d, X_out);
    PISM_CHK(ierr, "DMGetLocalVector");

    ierr = DMGlobalToLocalBegin(da_2d, X, INSERT_VALUES, *X_out);
    PISM_CHK(ierr, "DMGlobalToLocalBegin");

    ierr = DMGlobalToLocalEnd(da_2d, X, INSERT_VALUES, *X_out);
    PISM_CHK(ierr, "DMGlobalToLocalEnd");
  } else {
    *X_out = X;
  }

  ierr = DMDAVecGetArray(da_2d, *X_out, prm);
  PISM_CHK(ierr, "DMDAVecGetArray");
}

void Poisson3::end_2d_access(DM da, bool local, Vec *X_out, Parameters ***prm) {
  int ierr;

  DM da_2d;

  ierr = PetscObjectQuery((PetscObject)da, "2D_DM", (PetscObject*)&da_2d);
  PISM_CHK(ierr, "PetscObjectQuery");

  if (!da_2d) {
    throw RuntimeError(PISM_ERROR_LOCATION, "No 2D_DM composed with given DMDA");
  }

  ierr = DMDAVecRestoreArray(da_2d, *X_out, prm);
  PISM_CHK(ierr, "DMDAVecRestoreArray");

  if (local) {
    ierr = DMRestoreLocalVector(da_2d, X_out);
    PISM_CHK(ierr, "DMRestoreLocalVector");
  }
}

void Poisson3::begin_3d_access(DM da, bool local, Vec *X_out, double ****prm) {
  int ierr;

  DM  da_3d;
  Vec X;

  ierr = PetscObjectQuery((PetscObject)da, "3D_DM", (PetscObject*)&da_3d);
  PISM_CHK(ierr, "PetscObjectQuery");

  if (!da_3d) {
    throw RuntimeError(PISM_ERROR_LOCATION, "No 3D_DM composed with given DMDA");
  }

  ierr = PetscObjectQuery((PetscObject)da, "3D_Vec", (PetscObject*)&X);
  PISM_CHK(ierr, "PetscObjectQuery");

  if (!X) {
    throw RuntimeError(PISM_ERROR_LOCATION, "No 3D_Vec composed with given DMDA");
  }

  if (local) {
    ierr = DMGetLocalVector(da_3d, X_out);
    PISM_CHK(ierr, "DMGetLocalVector");

    ierr = DMGlobalToLocalBegin(da_3d, X, INSERT_VALUES, *X_out);
    PISM_CHK(ierr, "DMGlobalToLocalBegin");

    ierr = DMGlobalToLocalEnd(da_3d, X, INSERT_VALUES, *X_out);
    PISM_CHK(ierr, "DMGlobalToLocalEnd");
  } else {
    *X_out = X;
  }

  ierr = DMDAVecGetArray(da_3d, *X_out, prm);
  PISM_CHK(ierr, "DMDAVecGetArray");
}

void Poisson3::end_3d_access(DM da, bool local, Vec *X_out, double ****prm) {
  int ierr;

  DM da_3d;

  ierr = PetscObjectQuery((PetscObject)da, "3D_DM", (PetscObject*)&da_3d);
  PISM_CHK(ierr, "PetscObjectQuery");

  if (!da_3d) {
    throw RuntimeError(PISM_ERROR_LOCATION, "No 3D_DM composed with given DMDA");
  }

  ierr = DMDAVecRestoreArray(da_3d, *X_out, prm);
  PISM_CHK(ierr, "DMDAVecRestoreArray");

  if (local) {
    ierr = DMRestoreLocalVector(da_3d, X_out);
    PISM_CHK(ierr, "DMRestoreLocalVector");
  }
}

/*!
 * Bottom surface elevation
 */
static double b(double x, double y) {
  (void) x;
  (void) y;
  return -1.0 + x + y;
}

/*!
 * Domain thickness
 */
static double H(double x, double y) {
  return 1.0 + x*x + y*y;
}

/*!
 * Set 2D parameters on the finest grid.
 */
void Poisson3::init_2d_parameters() {

  DMDALocalInfo info;
  int ierr = DMDAGetLocalInfo(m_da, &info);
  PISM_CHK(ierr, "DMDAGetLocalInfo");

  // Compute grid spacing from domain dimensions and the grid size
  double
    Lx = m_grid->Lx(),
    Ly = m_grid->Ly(),
    dx = 2.0 * Lx / (info.mx - 1),
    dy = 2.0 * Ly / (info.my - 1);

  Parameters **P;
  Vec X;

  begin_2d_access(m_da, false, &X, &P);

  for (int j = info.ys; j < info.ys + info.ym; j++) {
    for (int i = info.xs; i < info.xs + info.xm; i++) {
      double x = xy(Lx, dx, i);
      double y = xy(Ly, dy, j);

      P[j][i].bed = b(x, y);
      P[j][i].thickness = H(x, y);
    }
  }

  end_2d_access(m_da, false, &X, &P);
}

/*!
 * Set 2D parameters on the finest grid.
 */
void Poisson3::init_3d_parameters() {

  DMDALocalInfo info;
  int ierr = DMDAGetLocalInfo(m_da, &info);
  PISM_CHK(ierr, "DMDAGetLocalInfo");

  // Compute grid spacing from domain dimensions and the grid size
  double
    Lx = m_grid->Lx(),
    Ly = m_grid->Ly(),
    dx = 2.0 * Lx / (info.mx - 1),
    dy = 2.0 * Ly / (info.my - 1);

  Parameters **P2;
  double ***P3;
  Vec X2, X3;

  begin_2d_access(m_da, false, &X2, &P2);
  begin_3d_access(m_da, false, &X3, &P3);

  for (int k = info.zs; k < info.zs + info.zm; k++) {
    for (int j = info.ys; j < info.ys + info.ym; j++) {
      for (int i = info.xs; i < info.xs + info.xm; i++) {
        double
          xx = xy(Lx, dx, i),
          yy = xy(Ly, dy, j),
          b  = P2[j][i].bed,
          H  = P2[j][i].thickness,
          zz = z(b, H, info.mz, k);

        P3[k][j][i] = F(xx, yy, zz);
      }
    }
  }

  end_3d_access(m_da, false, &X3, &P3);
  end_2d_access(m_da, false, &X2, &P2);
}

Poisson3::~Poisson3() {
  // empty
}

void Poisson3::exact_solution(IceModelVec3Custom &result) {
  IceModelVec::AccessList list{&result};

  // Compute grid spacing from domain dimensions and the grid size
  double
    Lx = m_grid->Lx(),
    Ly = m_grid->Ly(),
    dx = 2.0 * Lx / (m_grid->Mx() - 1),
    dy = 2.0 * Ly / (m_grid->My() - 1);

  Parameters **P;
  Vec X;

  begin_2d_access(m_da, false, &X, &P);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double
      xx = xy(Lx, dx, i),
      yy = xy(Ly, dy, j),
      b = P[j][i].bed,
      H = P[j][i].thickness;

    auto c = result.get_column(i, j);

    for (int k = 0; k < m_Mz; ++k) {
      double zz = z(b, H, m_Mz, k);

      c[k] = u_exact(xx, yy, zz);
    }
  }

  end_2d_access(m_da, false, &X, &P);
}

double Poisson3::error() const {
  IceModelVec3Custom difference(m_grid, "difference", "z_sigma",
                                m_exact->levels(), {});

  difference.copy_from(*m_exact);
  difference.add(-1.0, *m_solution);

  return difference.norm(NORM_INFINITY);
}

void Poisson3::update(const Inputs &inputs, bool) {
  (void) inputs;

  init_2d_parameters();
  init_3d_parameters();

  int ierr = 0;

  ierr = SNESSolve(m_snes, NULL, m_x); PISM_CHK(ierr, "SNESSolve");

  exact_solution(*m_exact);

  {
    double ***x = nullptr;
    ierr = DMDAVecGetArray(m_da, m_x, &x); PISM_CHK(ierr, "DMDAVecGetArray");

    IceModelVec::AccessList list{m_solution.get()};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      auto c = m_solution->get_column(i, j);

      for (int k = 0; k < m_Mz; ++k) {
        c[k] = x[k][j][i];
      }
    }

    ierr = DMDAVecRestoreArray(m_da, m_x, &x); PISM_CHK(ierr, "DMDAVecRestoreArray");
  }
}

IceModelVec3Custom::Ptr Poisson3::solution() const {
  return m_solution;
}

IceModelVec3Custom::Ptr Poisson3::exact() const {
  return m_exact;
}

} // end of namespace stressbalance
} // end of namespace pism
