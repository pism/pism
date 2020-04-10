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

static double xy(double L, double delta, int k) {
  return -L + k * delta;
}

static double z(double b, double H, int Mz, int k) {
  return b + H * k / (Mz - 1.0);
}

static bool dirichlet_node(const DMDALocalInfo *info, const fem::Element3::GlobalIndex& I) {
  return
    (I.i == 0 or I.i == info->mx - 1) or
    (I.j == 0 or I.j == info->my - 1) or
    (I.k == info->mz - 1);
}

static bool neumann_node(const DMDALocalInfo *info, const fem::Element3::GlobalIndex& I) {
  (void) info;
  return I.k == 0;
}

// Dirichlet BC and the exact solution
static double u_bc(double x, double y, double z) {
  return 2.0 * (1 + y) / ((3 + x) * (3 + x)  +  (1 + y) * (1 + y)) + (z + 1) * (z + 1);
}

// right hand side
static double F(double x, double y, double z) {
  (void) x;
  (void) y;
  (void) z;
  return -2.0;
}

// Neumann BC
static double G(double x, double y, double z) {
  (void) x;
  (void) y;
  (void) z;
  return 2.0;
}

// Bottom surface elevatipn
static double b(double x, double y) {
  (void) x;
  (void) y;
  return 0.0;
}

// Thickness
static double H(double x, double y) {
  return 1.0 + x*x + y*y;
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
            zz = z(b(xx, yy), H(xx, yy), info->mz, k);

          R[k][j][i] = u_bc(xx, yy, zz) - x[k][j][i];
        } else {
          R[k][j][i] = 0.0;
        }
      }
    }
  }

  // values at element nodes
  int Nk = E.n_chi();
  std::vector<double>
    x_nodal(Nk), y_nodal(Nk), z_nodal(Nk),
    R_nodal(Nk), u_nodal(Nk);

  // values at quadrature points
  int Nq = E.n_pts();
  std::vector<double> u(Nq), u_x(Nq), u_y(Nq), u_z(Nq);
  std::vector<double> xq(Nq), yq(Nq), zq(Nq);

  // loop over all the elements that have at least one owned node
  for (int k = info->gzs; k < info->gzs + info->gzm - 1; k++) {
    for (int j = info->gys; j < info->gys + info->gym - 1; j++) {
      for (int i = info->gxs; i < info->gxs + info->gxm - 1; i++) {

        // reset element residual to zero in preparation
        for (int n = 0; n < Nk; ++n) {
          R_nodal[n] = 0.0;
        }

        // compute coordinates of the nodes of this element
        for (int n = 0; n < Nk; ++n) {
          auto I = E.local_to_global(i, j, k, n);

          x_nodal[n] = xy(Lx, dx, I.i);
          y_nodal[n] = xy(Ly, dy, I.j);
          z_nodal[n] = z(b(x_nodal[n], y_nodal[n]),
                         H(x_nodal[n], y_nodal[n]),
                         info->mz, I.k);
        }

        E.reset(i, j, k, z_nodal);

        // get nodal values of u
        E.nodal_values(x, u_nodal.data());

        // take care of Dirichlet BC: don't contribute to Dirichlet nodes and set nodal
        // values of the current iterate to the BC value
        for (int n = 0; n < Nk; ++n) {
          auto I = E.local_to_global(n);
          if (dirichlet_node(info, I)) {
            E.mark_row_invalid(n);
            u_nodal[n] = u_bc(x_nodal[n], y_nodal[n], z_nodal[n]);
          }
        }

        // evaluate u and its partial derivatives at quadrature points
        E.evaluate(u_nodal.data(), u.data(), u_x.data(), u_y.data(), u_z.data());

        // coordinates of quadrature points
        E.evaluate(x_nodal.data(), xq.data());
        E.evaluate(y_nodal.data(), yq.data());
        E.evaluate(z_nodal.data(), zq.data());

        // loop over all quadrature points
        for (int q = 0; q < Nq; ++q) {
          auto W = E.weight(q);

          // loop over all test functions
          for (int t = 0; t < Nk; ++t) {
            const auto &psi = E.chi(q, t);

            R_nodal[t] += W * (u_x[q] * psi.dx + u_y[q] * psi.dy + u_z[q] * psi.dz
                               - F(xq[q], yq[q], zq[q]) * psi.val);
          }
        }

        // compute the Neumann BC function at the nodes of this element
        //
        // I could compute it directly at quadrature points but I don't want to compute
        // x,y,z coordinates of quadrature points on each face
        std::vector<double> g_nodal(Nk);
        for (int n = 0; n < Nk; ++n) {
          g_nodal[n] = G(x_nodal[n], y_nodal[n], z_nodal[n]);
        }

        // loop over all faces
        for (int face = 0; face < fem::q13d::n_faces; ++face) {
          auto nodes = fem::q13d::incident_nodes[face];
          // loop over all nodes corresponding to this face. A face is a part of the
          // Neumann boundary if all four nodes are Neumann nodes. If a node is *both* a
          // Neumann and a Dirichlet node (this may happen), then we treat it as a Neumann
          // node here: add_contribution() will do the right thing.
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

          std::vector<double> g(E_face.n_pts());
          E_face.evaluate(g_nodal.data(), g.data());

          for (int q = 0; q < E_face.n_pts(); ++q) {
            auto W = E_face.weight(q);

            for (int t = 0; t < Nk; ++t) {
              auto psi = E_face.chi(q, t);

              R_nodal[t] += W * psi * g[q];
            }

          }
        } // end of the loop over element faces

        E.add_contribution(R_nodal.data(), R);
      } // end of the loop over i
    } // end of the loop over j
  } // end of the loop over k
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

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double
      xx = xy(Lx, dx, i),
      yy = xy(Ly, dy, j);

    auto c = result.get_column(i, j);

    for (int k = 0; k < m_Mz; ++k) {
      double zz = z(b(xx, yy), H(xx, yy), m_Mz, k);

      c[k] = u_bc(xx, yy, zz);
    }
  }
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


