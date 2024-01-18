// Copyright (C) 2004--2024 Constantine Khroulev, Ed Bueler and Jed Brown
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

#include <cassert>
#include <stdexcept>

#include "pism/basalstrength/basal_resistance.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/stressbalance/ssa/SSAFD.hh"
#include "pism/stressbalance/ssa/SSAFD_diagnostics.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace stressbalance {

SSAFD::KSPFailure::KSPFailure(const char* reason)
  : RuntimeError(ErrorLocation(), std::string("SSAFD KSP (linear solver) failed: ") + reason){
  // empty
}

SSAFD::PicardFailure::PicardFailure(const std::string &message)
  : RuntimeError(ErrorLocation(), "SSAFD Picard iterations failed: " + message) {
  // empty
}

SSA* SSAFDFactory(std::shared_ptr<const Grid> g) {
  return new SSAFD(g);
}

/*!
Because the FD implementation of the SSA uses Picard iteration, a PETSc KSP
and Mat are used directly.  In particular we set up \f$A\f$
(Mat m_A) and a \f$b\f$ (= Vec m_b) and iteratively solve
linear systems
  \f[ A x = b \f]
where \f$x\f$ (= Vec SSAX).  A PETSc SNES object is never created.
 */
SSAFD::SSAFD(std::shared_ptr<const Grid> grid)
  : SSA(grid),
    m_hardness(grid, "hardness"),
    m_nuH(grid, "nuH"),
    m_nuH_old(grid, "nuH_old"),
    m_work(grid, "work_vector", array::WITH_GHOSTS,
           1 /* stencil width */),
    m_mask(m_grid, "ssafd_cell_type"),
    m_rhs(grid, "right_hand_side"),
    m_taud(m_grid, "taud"),
    m_velocity_old(grid, "velocity_old"),
    m_scaling(1e9)  // comparable to typical beta for an ice stream;
{

  m_velocity_old.metadata(0)
      .long_name("old SSA velocity field; used for re-trying with a different epsilon")
      .units("m s-1");

  m_hardness.metadata(0)
      .long_name("vertically-averaged ice hardness")
      .set_units_without_validation(pism::printf("Pa s^(1/%f)", m_flow_law->exponent()));

  m_nuH.metadata(0)
      .long_name("ice thickness times effective viscosity")
      .units("Pa s m");

  m_nuH_old.metadata(0)
      .long_name("ice thickness times effective viscosity (before an update)")
      .units("Pa s m");

  m_taud.metadata(0)
      .long_name("X-component of the driving shear stress at the base of ice")
      .units("Pa");
  m_taud.metadata(1)
      .long_name("Y-component of the driving shear stress at the base of ice")
      .units("Pa");

  m_work.metadata(0).long_name("temporary storage used to compute nuH");

  // grounded_dragging_floating integer mask
  m_mask.metadata(0)
      .long_name("ice-type (ice-free/grounded/floating/ocean) integer mask");
  m_mask.metadata()["flag_values"]   = { MASK_ICE_FREE_BEDROCK, MASK_GROUNDED, MASK_FLOATING,
                                         MASK_ICE_FREE_OCEAN };
  m_mask.metadata()["flag_meanings"] = "ice_free_bedrock grounded_ice floating_ice ice_free_ocean";

  // The nuH viewer:
  m_view_nuh        = false;
  m_nuh_viewer_size = 300;

  // PETSc objects and settings
  {
    auto dm = m_velocity_global.dm();

    PetscErrorCode ierr;
    ierr = DMSetMatType(*dm, MATAIJ);
    PISM_CHK(ierr, "DMSetMatType");

    ierr = DMCreateMatrix(*dm, m_A.rawptr());
    PISM_CHK(ierr, "DMCreateMatrix");

    ierr = KSPCreate(m_grid->com, m_KSP.rawptr());
    PISM_CHK(ierr, "KSPCreate");

    ierr = KSPSetOptionsPrefix(m_KSP, "ssafd_");
    PISM_CHK(ierr, "KSPSetOptionsPrefix");

    // Use non-zero initial guess (i.e. SSA velocities from the last
    // solve() call).
    ierr = KSPSetInitialGuessNonzero(m_KSP, PETSC_TRUE);
    PISM_CHK(ierr, "KSPSetInitialGuessNonzero");

    // Use the initial residual norm.
    ierr = KSPConvergedDefaultSetUIRNorm(m_KSP);
    PISM_CHK(ierr, "KSPConvergedDefaultSetUIRNorm");
  }
}

const array::Vector& SSAFD::driving_stress() const {
  return m_taud;
}

//! @note Uses `PetscErrorCode` *intentionally*.
void SSAFD::pc_setup_bjacobi() {
  PetscErrorCode ierr;
  PC pc;

  ierr = KSPSetType(m_KSP, KSPGMRES);
  PISM_CHK(ierr, "KSPSetType");

  ierr = KSPSetOperators(m_KSP, m_A, m_A);
  PISM_CHK(ierr, "KSPSetOperators");

  // Get the PC from the KSP solver:
  ierr = KSPGetPC(m_KSP, &pc);
  PISM_CHK(ierr, "KSPGetPC");

  // Set the PC type:
  ierr = PCSetType(pc, PCBJACOBI);
  PISM_CHK(ierr, "PCSetType");

  // Process options:
  ierr = KSPSetFromOptions(m_KSP);
  PISM_CHK(ierr, "KSPSetFromOptions");
}

//! @note Uses `PetscErrorCode` *intentionally*.
void SSAFD::pc_setup_asm() {
  PetscErrorCode ierr;
  PC pc, sub_pc;

  // Set parameters equivalent to
  // -ksp_type gmres -ksp_norm_type unpreconditioned -ksp_pc_side right -pc_type asm -sub_pc_type lu

  ierr = KSPSetType(m_KSP, KSPGMRES);
  PISM_CHK(ierr, "KSPSetType");

  ierr = KSPSetOperators(m_KSP, m_A, m_A);
  PISM_CHK(ierr, "KSPSetOperators");

  // Switch to using the "unpreconditioned" norm.
  ierr = KSPSetNormType(m_KSP, KSP_NORM_UNPRECONDITIONED);
  PISM_CHK(ierr, "KSPSetNormType");

  // Switch to "right" preconditioning.
  ierr = KSPSetPCSide(m_KSP, PC_RIGHT);
  PISM_CHK(ierr, "KSPSetPCSide");

  // Get the PC from the KSP solver:
  ierr = KSPGetPC(m_KSP, &pc);
  PISM_CHK(ierr, "KSPGetPC");

  // Set the PC type:
  ierr = PCSetType(pc, PCASM);
  PISM_CHK(ierr, "PCSetType");

  // Set the sub-KSP object to "preonly"
  KSP *sub_ksp;
  ierr = PCSetUp(pc);
  PISM_CHK(ierr, "PCSetUp");

  ierr = PCASMGetSubKSP(pc, NULL, NULL, &sub_ksp);
  PISM_CHK(ierr, "PCASMGetSubKSP");

  ierr = KSPSetType(*sub_ksp, KSPPREONLY);
  PISM_CHK(ierr, "KSPSetType");

  // Set the PC of the sub-KSP to "LU".
  ierr = KSPGetPC(*sub_ksp, &sub_pc);
  PISM_CHK(ierr, "KSPGetPC");

  ierr = PCSetType(sub_pc, PCLU);
  PISM_CHK(ierr, "PCSetType");

  // Let the user override all this:
  // Process options:
  ierr = KSPSetFromOptions(m_KSP);
  PISM_CHK(ierr, "KSPSetFromOptions");
}

void SSAFD::init_impl() {
  SSA::init_impl();

  m_log->message(2, "  [using the KSP-based finite difference implementation]\n");

  // options
  options::Integer viewer_size("-ssa_nuh_viewer_size", "nuH viewer size", m_nuh_viewer_size);
  m_nuh_viewer_size = viewer_size;
  m_view_nuh        = options::Bool("-ssa_view_nuh", "Enable the SSAFD nuH runtime viewer");

  if (m_config->get_flag("stress_balance.calving_front_stress_bc")) {
    m_log->message(2, "  using PISM-PIK calving-front stress boundary condition ...\n");
  }

  m_default_pc_failure_count     = 0;
  m_default_pc_failure_max_count = 5;
}

//! \brief Computes the right-hand side ("rhs") of the linear problem for the
//! Picard iteration and finite-difference implementation of the SSA equations.
/*!
The right side of the SSA equations is just the driving stress term
   \f[ - \rho g H \nabla h. \f]
The basal stress is put on the left side of the system.  This method builds the
discrete approximation of the right side.  For more about the discretization
of the SSA equations, see comments for assemble_matrix().

The values of the driving stress on the i,j grid come from a call to
compute_driving_stress().

In the case of Dirichlet boundary conditions, the entries on the right-hand side
come from known velocity values.  The fields m_bc_values and m_bc_mask are used for
this.
 */
void SSAFD::assemble_rhs(const Inputs &inputs, const array::CellType1 &cell_type,
                         array::Vector &result) {
  using mask::ice_free;
  using mask::ice_free_land;
  using mask::ice_free_ocean;

  const array::Scalar1 &bed                  = inputs.geometry->bed_elevation;
  const array::Scalar &thickness             = inputs.geometry->ice_thickness,
                      &surface               = inputs.geometry->ice_surface_elevation,
                      &sea_level             = inputs.geometry->sea_level_elevation,
                      *water_column_pressure = inputs.water_column_pressure;

  const double dx = m_grid->dx(), dy = m_grid->dy(),
               standard_gravity = m_config->get_number("constants.standard_gravity"),
               rho_ocean        = m_config->get_number("constants.sea_water.density"),
               rho_ice          = m_config->get_number("constants.ice.density");

  // This constant is for debugging: simulations should not depend on the choice of
  // velocity used in ice-free areas.
  const Vector2d ice_free_velocity(0.0, 0.0);

  const bool use_cfbc       = m_config->get_flag("stress_balance.calving_front_stress_bc"),
             flow_line_mode = m_config->get_flag("stress_balance.ssa.fd.flow_line_mode");

  // FIXME: bedrock_boundary is a misleading name
  bool bedrock_boundary = m_config->get_flag("stress_balance.ssa.dirichlet_bc");

  array::AccessScope list{ &m_taud, &result };

  if (inputs.bc_values != nullptr and inputs.bc_mask != nullptr) {
    list.add({ inputs.bc_values, inputs.bc_mask });
  }

  if (use_cfbc) {
    list.add({ &thickness, &bed, &surface, &cell_type, &sea_level });
  }

  if (use_cfbc and (water_column_pressure != nullptr)) {
    list.add(*water_column_pressure);
  }

  result.set(0.0);

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    Vector2d taud = m_taud(i, j);

    if (flow_line_mode) {
      // no cross-flow driving stress in the flow line mode
      taud.v = 0.0;
    }

    if ((inputs.bc_values != nullptr) and inputs.bc_mask->as_int(i, j) == 1) {
      result(i, j).u = m_scaling * (*inputs.bc_values)(i, j).u;
      result(i, j).v = m_scaling * (*inputs.bc_values)(i, j).v;
      continue;
    }

    if (use_cfbc) {
      double H_ij = thickness(i, j);

      auto M = cell_type.star_int(i, j);

      // Note: this sets velocities at both ice-free ocean and ice-free
      // bedrock to zero. This means that we need to set boundary conditions
      // at both ice/ice-free-ocean and ice/ice-free-bedrock interfaces below
      // to be consistent.
      if (ice_free(M.c)) {
        result(i, j) = m_scaling * ice_free_velocity;
        continue;
      }

      if (is_marginal(i, j, bedrock_boundary)) {
        // weights at the west, east, south, and north cell faces
        int W = 0, E = 0, S = 0, N = 0;
        // direct neighbors
        // NOLINTBEGIN(readability-braces-around-statements)
        if (bedrock_boundary) {
          if (ice_free_ocean(M.e))
            E = 1;
          if (ice_free_ocean(M.w))
            W = 1;
          if (ice_free_ocean(M.n))
            N = 1;
          if (ice_free_ocean(M.s))
            S = 1;
        } else {
          if (ice_free(M.e))
            E = 1;
          if (ice_free(M.w))
            W = 1;
          if (ice_free(M.n))
            N = 1;
          if (ice_free(M.s))
            S = 1;
        }
        // NOLINTEND(readability-braces-around-statements)

        double P_ice = 0.5 * rho_ice * standard_gravity * H_ij, P_water = 0.0;

        if (water_column_pressure != nullptr) {
          P_water = (*water_column_pressure)(i, j);
        } else {
          P_water = pism::average_water_column_pressure(H_ij, bed(i, j), sea_level(i, j), rho_ice,
                                                        rho_ocean, standard_gravity);
        }

        double delta_p = H_ij * (P_ice - P_water);

        if (grid::domain_edge(*m_grid, i, j) and not(flow_line_mode or mask::grounded(M.c))) {
          // In regional setups grounded ice may extend to the edge of the domain. This
          // condition ensures that at a domain edge the ice behaves as if it extends past
          // the edge without a change in geometry.
          delta_p = 0.0;
        }

        {
          // fjord walls, nunataks, etc
          //
          // Override weights if we are at the margin and the grid cell on the other side
          // of the interface is ice-free and above the level of the ice surface.
          //
          // This effectively sets the pressure difference at the corresponding interface
          // to zero, which is exactly what we need.
          auto b   = bed.star(i, j);
          double h = surface(i, j);

          if (ice_free(M.n) and b.n > h) {
            N = 0;
          }
          if (ice_free(M.e) and b.e > h) {
            E = 0;
          }
          if (ice_free(M.s) and b.s > h) {
            S = 0;
          }
          if (ice_free(M.w) and b.w > h) {
            W = 0;
          }
        }

        // Note that if the current cell is "marginal" but not a CFBC
        // location, the following two lines are equaivalent to the "usual
        // case" below.
        //
        // Note: signs below (+E, -W, etc) are explained by directions of outward
        // normal vectors at corresponding cell faces.
        result(i, j).u = taud.u + (E - W) * delta_p / dx;
        result(i, j).v = taud.v + (N - S) * delta_p / dy;

        continue;
      } // end of "if (is_marginal(i, j))"

      // If we reached this point, then CFBC are enabled, but we are in the
      // interior of a sheet or shelf. See "usual case" below.

    } // end of "if (use_cfbc)"

    // usual case: use already computed driving stress
    result(i, j) = taud;
  }
}

static void set_diagonal_matrix_entry(Mat A, int i, int j, int component, double value) {
  MatStencil row, col;

  row.i = i;
  row.j = j;
  row.c = component;

  col.i = i;
  col.j = j;
  col.c = component;

  PetscErrorCode ierr = MatSetValuesStencil(A, 1, &row, 1, &col, &value, INSERT_VALUES);
  PISM_CHK(ierr, "MatSetValuesStencil");
}


//! \brief Assemble the left-hand side matrix for the KSP-based, Picard iteration,
//! and finite difference implementation of the SSA equations.
/*!
Recall the SSA equations are
\f{align*}
 - 2 \left[\nu H \left(2 u_x + v_y\right)\right]_x
        - \left[\nu H \left(u_y + v_x\right)\right]_y
        - \tau_{(b)1}  &= - \rho g H h_x, \\
   - \left[\nu H \left(u_y + v_x\right)\right]_x
        - 2 \left[\nu H \left(u_x + 2 v_y\right)\right]_y
        - \tau_{(b)2}  &= - \rho g H h_y,
\f}
where \f$u\f$ is the \f$x\f$-component of the velocity and \f$v\f$ is the
\f$y\f$-component of the velocity.

The coefficient \f$\nu\f$ is the vertically-averaged effective viscosity.
(The product \f$\nu H\f$ is computed by compute_nuH_staggered().)
The Picard iteration idea is that, to solve the nonlinear equations in which
the effective viscosity depends on the velocity, we freeze the effective
viscosity using its value at the current estimate of the velocity and we solve
the linear equations which come from this viscosity.  In abstract symbols, the
Picard iteration replaces the above nonlinear SSA equations by a sequence of
linear problems

\f[ A(U^{(k)}) U^{(k+1)} = b \f]

where \f$A(U)\f$ is the matrix from discretizing the SSA equations supposing
the viscosity is a function of known velocities \f$U\f$, and where \f$U^{(k)}\f$
denotes the \f$k\f$th iterate in the process.  The current method assembles \f$A(U)\f$.

For ice shelves \f$\tau_{(b)i} = 0\f$ [\ref MacAyealetal].
For ice streams with a basal till modelled as a plastic material,
\f$\tau_{(b)i} = - \tau_c u_i/|\mathbf{u}|\f$ where
\f$\mathbf{u} = (u,v)\f$, \f$|\mathbf{u}| = \left(u^2 + v^2\right)^{1/2}\f$,
where \f$\tau_c(t,x,y)\f$ is the yield stress of the till [\ref SchoofStream].
More generally, ice streams can be modeled with a pseudo-plastic basal till;
see IceModel::initBasalTillModel() and IceModel::updateYieldStressUsingBasalWater()
and reference [\ref BKAJS].  The pseudo-plastic till model includes all power law
sliding relations and the linearly-viscous model for sliding,
\f$\tau_{(b)i} = - \beta u_i\f$ where \f$\beta\f$ is the basal drag
(friction) parameter [\ref MacAyeal].  In any case, PISM assumes that the basal shear
stress can be factored this way, *even if the coefficient depends on the
velocity*, \f$\beta(u,v)\f$.  Such factoring is possible even in the case of
(regularized) plastic till.  This scalar coefficient \f$\beta\f$ is what is
returned by IceBasalResistancePlasticLaw::drag().

Note that the basal shear stress appears on the \em left side of the
linear system we actually solve. We believe this is crucial, because
of its effect on the spectrum of the linear approximations of each
stage. The effect on spectrum is clearest in the linearly-viscous till
case but there seems to be an analogous effect in the plastic till case.

This method assembles the matrix for the left side of the above SSA equations.
The numerical method is finite difference.  Suppose we use difference notation
\f$\delta_{+x}f^{i,j} = f^{i+1,j}-f^{i,j}\f$,
\f$\delta_{-x}f^{i,j} = f^{i,j}-f^{i-1,j}\f$, and
\f$\Delta_{x}f^{i,j} = f^{i+1,j}-f^{i-1,j}\f$, and corresponding notation for
\f$y\f$ differences, and that we write \f$N = \nu H\f$ then the first of the
two "concrete" SSA equations above has this discretization:
\f{align*}
- &2 \frac{N^{i+\frac{1}{2},j}}{\Delta x} \left[2\frac{\delta_{+x}u^{i,j}}{\Delta x} + \frac{\Delta_{y} v^{i+1,j} + \Delta_{y} v^{i,j}}{4 \Delta y}\right] + 2 \frac{N^{i-\frac{1}{2},j}}{\Delta x} \left[2\frac{\delta_{-x}u^{i,j}}{\Delta x} + \frac{\Delta_y v^{i,j} + \Delta_y v^{i-1,j}}{4 \Delta y}\right] \\
&\qquad- \frac{N^{i,j+\frac{1}{2}}}{\Delta y} \left[\frac{\delta_{+y} u^{i,j}}{\Delta y} + \frac{\Delta_x v^{i,j+1} + \Delta_x v^{i,j}}{4 \Delta x}\right] + \frac{N^{i,j-\frac{1}{2}}}{\Delta y} \left[\frac{\delta_{-y}u^{i,j}}{\Delta y} + \frac{\Delta_x v^{i,j} + \Delta_x v^{i,j-1}}{4 \Delta x}\right] - \tau_{(b)1}^{i,j} = - \rho g H^{i,j} \frac{\Delta_x h^{i,j}}{2\Delta x}.
\f}
As a picture, see Figure \ref ssastencil.

\image html ssastencil.png "\b ssastencil:  Stencil for our finite difference discretization of the first of the two scalar SSA equations.  Triangles show staggered grid points where N = nu * H is evaluated.  Circles and squares show where u and v are approximated, respectively."
\anchor ssastencil

It follows immediately that the matrix we assemble in the current method has
13 nonzeros entries per row because, for this first SSA equation, there are 5
grid values of \f$u\f$ and 8 grid values of \f$v\f$ used in this scheme.  For
the second equation we also have 13 nonzeros per row.

FIXME:  document use of DAGetMatrix and MatStencil and MatSetValuesStencil

*/
void SSAFD::assemble_matrix(const Inputs &inputs, const array::Vector &velocity,
                            bool include_basal_shear, Mat A) {
  using mask::grounded_ice;
  using mask::ice_free;
  using mask::ice_free_land;
  using mask::ice_free_ocean;
  using mask::icy;

  const int diag_u = 4;
  const int diag_v = 13;

  PetscErrorCode ierr = 0;

  const array::Scalar1 &thickness        = inputs.geometry->ice_thickness,
                       &bed              = inputs.geometry->bed_elevation;
  const array::Scalar &surface           = inputs.geometry->ice_surface_elevation,
                      &grounded_fraction = inputs.geometry->cell_grounded_fraction,
                      &tauc              = *inputs.basal_yield_stress;

  const double dx = m_grid->dx(), dy = m_grid->dy(),
               beta_lateral_margin = m_config->get_number("basal_resistance.beta_lateral_margin"),
               beta_ice_free_bedrock =
                   m_config->get_number("basal_resistance.beta_ice_free_bedrock");

  const bool
      // FIXME: bedrock_boundary is a misleading name
      bedrock_boundary = m_config->get_flag("stress_balance.ssa.dirichlet_bc"),
      flow_line_mode   = m_config->get_flag("stress_balance.ssa.fd.flow_line_mode"),
      use_cfbc         = m_config->get_flag("stress_balance.calving_front_stress_bc"),
      replace_zero_diagonal_entries =
          m_config->get_flag("stress_balance.ssa.fd.replace_zero_diagonal_entries");

  ierr = MatZeroEntries(A);
  PISM_CHK(ierr, "MatZeroEntries");

  array::AccessScope list{ &m_nuH, &tauc, &velocity, &m_mask, &bed, &surface };

  if (inputs.bc_values != nullptr && inputs.bc_mask != nullptr) {
    list.add(*inputs.bc_mask);
  }

  const bool sub_gl = m_config->get_flag("geometry.grounded_cell_fraction");
  if (sub_gl) {
    list.add(grounded_fraction);
  }

  // handles friction of the ice cell along ice-free bedrock margins when bedrock higher than ice
  // surface (in simplified setups)
  bool lateral_drag_enabled = m_config->get_flag("stress_balance.ssa.fd.lateral_drag.enabled");
  if (lateral_drag_enabled) {
    list.add({ &thickness, &bed, &surface });
  }
  double lateral_drag_viscosity =
      m_config->get_number("stress_balance.ssa.fd.lateral_drag.viscosity");
  double HminFrozen = 0.0;

  /* matrix assembly loop */
  ParallelSection loop(m_grid->com);
  try {
    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      // Handle the easy case: provided Dirichlet boundary conditions
      if (inputs.bc_values != nullptr && inputs.bc_mask != nullptr &&
          inputs.bc_mask->as_int(i, j) == 1) {
        // set diagonal entry to one (scaled); RHS entry will be known velocity;
        set_diagonal_matrix_entry(A, i, j, 0, m_scaling);
        set_diagonal_matrix_entry(A, i, j, 1, m_scaling);
        continue;
      }

      /* Provide shorthand for the following staggered coefficients  nu H:
       *      c_n
       *  c_w     c_e
       *      c_s
       */
      // const
      double c_w = m_nuH(i - 1, j, 0);
      double c_e = m_nuH(i, j, 0);
      double c_s = m_nuH(i, j - 1, 1);
      double c_n = m_nuH(i, j, 1);

      if (lateral_drag_enabled) {
        // if option is set, the viscosity at ice-bedrock boundary layer will
        // be prescribed and is a temperature-independent free (user determined) parameter

        // direct neighbors
        auto M   = m_mask.star_int(i, j);
        auto H   = thickness.star(i, j);
        auto b   = bed.star(i, j);
        double h = surface(i, j);

        if (H.c > HminFrozen) {
          if (b.w > h and ice_free_land(M.w)) {
            c_w = lateral_drag_viscosity * 0.5 * (H.c + H.w);
          }
          if (b.e > h and ice_free_land(M.e)) {
            c_e = lateral_drag_viscosity * 0.5 * (H.c + H.e);
          }
          if (b.n > h and ice_free_land(M.n)) {
            c_n = lateral_drag_viscosity * 0.5 * (H.c + H.n);
          }
          if (b.s > h and ice_free_land(M.s)) {
            c_s = lateral_drag_viscosity * 0.5 * (H.c + H.s);
          }
        }
      }

      // We use DAGetMatrix to obtain the SSA matrix, which means that all 18
      // non-zeros get allocated, even though we use only 13 (or 14). The
      // remaining 5 (or 4) coefficients are zeros, but we set them anyway,
      // because this makes the code easier to understand.
      const int n_nonzeros = 18;
      MatStencil row, col[n_nonzeros];

      // |-----+-----+---+-----+-----|
      // | NW  | NNW | N | NNE | NE  |
      // | WNW |     | | |     | ENE |
      // | W   |-----|-o-|-----| E   |
      // | WSW |     | | |     | ESE |
      // | SW  | SSW | S | SSE | SE  |
      // |-----+-----+---+-----+-----|
      //
      // We use compass rose notation for weights corresponding to interfaces between
      // cells around the current one (i, j). Here N corresponds to the interface between
      // the cell (i, j) and the one to the north of it.
      int N = 1, E = 1, S = 1, W = 1;

      // Similarly, we use compass rose notation for weights used to switch between
      // centered and one-sided finite differences. Here NNE is the interface between
      // cells N and NE, ENE - between E and NE, etc.
      int NNW = 1, NNE = 1, SSW = 1, SSE = 1;
      int WNW = 1, ENE = 1, WSW = 1, ESE = 1;

      int M_ij = m_mask.as_int(i, j);

      if (use_cfbc) {
        auto M = m_mask.box_int(i, j);

        // Note: this sets velocities at both ice-free ocean and ice-free
        // bedrock to zero. This means that we need to set boundary conditions
        // at both ice/ice-free-ocean and ice/ice-free-bedrock interfaces below
        // to be consistent.
        if (ice_free(M.c)) {
          set_diagonal_matrix_entry(A, i, j, 0, m_scaling);
          set_diagonal_matrix_entry(A, i, j, 1, m_scaling);
          continue;
        }

        if (is_marginal(i, j, bedrock_boundary)) {
          // If at least one of the following four conditions is "true", we're
          // at a CFBC location.
          // NOLINTBEGIN(readability-braces-around-statements)
          if (bedrock_boundary) {

            if (ice_free_ocean(M.e))
              E = 0;
            if (ice_free_ocean(M.w))
              W = 0;
            if (ice_free_ocean(M.n))
              N = 0;
            if (ice_free_ocean(M.s))
              S = 0;

            // decide whether to use centered or one-sided differences
            if (ice_free_ocean(M.n) || ice_free_ocean(M.ne))
              NNE = 0;
            if (ice_free_ocean(M.e) || ice_free_ocean(M.ne))
              ENE = 0;
            if (ice_free_ocean(M.e) || ice_free_ocean(M.se))
              ESE = 0;
            if (ice_free_ocean(M.s) || ice_free_ocean(M.se))
              SSE = 0;
            if (ice_free_ocean(M.s) || ice_free_ocean(M.sw))
              SSW = 0;
            if (ice_free_ocean(M.w) || ice_free_ocean(M.sw))
              WSW = 0;
            if (ice_free_ocean(M.w) || ice_free_ocean(M.nw))
              WNW = 0;
            if (ice_free_ocean(M.n) || ice_free_ocean(M.nw))
              NNW = 0;

          } else { // if (not bedrock_boundary)

            if (ice_free(M.e))
              E = 0;
            if (ice_free(M.w))
              W = 0;
            if (ice_free(M.n))
              N = 0;
            if (ice_free(M.s))
              S = 0;

            // decide whether to use centered or one-sided differences
            if (ice_free(M.n) || ice_free(M.ne))
              NNE = 0;
            if (ice_free(M.e) || ice_free(M.ne))
              ENE = 0;
            if (ice_free(M.e) || ice_free(M.se))
              ESE = 0;
            if (ice_free(M.s) || ice_free(M.se))
              SSE = 0;
            if (ice_free(M.s) || ice_free(M.sw))
              SSW = 0;
            if (ice_free(M.w) || ice_free(M.sw))
              WSW = 0;
            if (ice_free(M.w) || ice_free(M.nw))
              WNW = 0;
            if (ice_free(M.n) || ice_free(M.nw))
              NNW = 0;

          } // end of the else clause following "if (bedrock_boundary)"
          // NOLINTEND(readability-braces-around-statements)
        } // end of "if (is_marginal(i, j, bedrock_boundary))"
      }   // end of "if (use_cfbc)"

      /* begin Maxima-generated code */
      const double dx2 = dx * dx, dy2 = dy * dy, d4 = 4 * dx * dy, d2 = 2 * dx * dy;

      /* Coefficients of the discretization of the first equation; u first, then v. */
      double eq1[] = {
        0,
        -c_n * N / dy2,
        0,
        -4 * c_w * W / dx2,
        (c_n * N + c_s * S) / dy2 + (4 * c_e * E + 4 * c_w * W) / dx2,
        -4 * c_e * E / dx2,
        0,
        -c_s * S / dy2,
        0,
        c_w * W * WNW / d2 + c_n * NNW * N / d4,
        (c_n * NNE * N - c_n * NNW * N) / d4 + (c_w * W * N - c_e * E * N) / d2,
        -c_e * E * ENE / d2 - c_n * NNE * N / d4,
        (c_w * W * WSW - c_w * W * WNW) / d2 + (c_n * W * N - c_s * W * S) / d4,
        (c_n * E * N - c_n * W * N - c_s * E * S + c_s * W * S) / d4 +
            (c_e * E * N - c_w * W * N - c_e * E * S + c_w * W * S) / d2,
        (c_e * E * ENE - c_e * E * ESE) / d2 + (c_s * E * S - c_n * E * N) / d4,
        -c_w * W * WSW / d2 - c_s * SSW * S / d4,
        (c_s * SSW * S - c_s * SSE * S) / d4 + (c_e * E * S - c_w * W * S) / d2,
        c_e * E * ESE / d2 + c_s * SSE * S / d4,
      };

      /* Coefficients of the discretization of the second equation; u first, then v. */
      double eq2[] = {
        c_w * W * WNW / d4 + c_n * NNW * N / d2,
        (c_n * NNE * N - c_n * NNW * N) / d2 + (c_w * W * N - c_e * E * N) / d4,
        -c_e * E * ENE / d4 - c_n * NNE * N / d2,
        (c_w * W * WSW - c_w * W * WNW) / d4 + (c_n * W * N - c_s * W * S) / d2,
        (c_n * E * N - c_n * W * N - c_s * E * S + c_s * W * S) / d2 +
            (c_e * E * N - c_w * W * N - c_e * E * S + c_w * W * S) / d4,
        (c_e * E * ENE - c_e * E * ESE) / d4 + (c_s * E * S - c_n * E * N) / d2,
        -c_w * W * WSW / d4 - c_s * SSW * S / d2,
        (c_s * SSW * S - c_s * SSE * S) / d2 + (c_e * E * S - c_w * W * S) / d4,
        c_e * E * ESE / d4 + c_s * SSE * S / d2,
        0,
        -4 * c_n * N / dy2,
        0,
        -c_w * W / dx2,
        (4 * c_n * N + 4 * c_s * S) / dy2 + (c_e * E + c_w * W) / dx2,
        -c_e * E / dx2,
        0,
        -4 * c_s * S / dy2,
        0,
      };

      /* i indices */
      const int I[] = {
        i - 1, i, i + 1, i - 1, i, i + 1, i - 1, i, i + 1,
        i - 1, i, i + 1, i - 1, i, i + 1, i - 1, i, i + 1,
      };

      /* j indices */
      const int J[] = {
        j + 1, j + 1, j + 1, j, j, j, j - 1, j - 1, j - 1,
        j + 1, j + 1, j + 1, j, j, j, j - 1, j - 1, j - 1,
      };

      /* component indices */
      const int C[] = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
      };
      /* end Maxima-generated code */

      /* Dragging ice experiences friction at the bed determined by the
       *    IceBasalResistancePlasticLaw::drag() methods.  These may be a plastic,
       *    pseudo-plastic, or linear friction law.  Dragging is done implicitly
       *    (i.e. on left side of SSA eqns).  */
      double beta_u = 0.0, beta_v = 0.0;
      if (include_basal_shear) {
        double beta = 0.0;
        switch (M_ij) {
        case MASK_ICE_FREE_BEDROCK: {
          // apply drag even in this case, to help with margins; note ice free areas may
          // already have a strength extension
          beta = beta_ice_free_bedrock;
          break;
        }
        case MASK_FLOATING: {
          double scaling = sub_gl ? grounded_fraction(i, j) : 0.0;
          beta = scaling * m_basal_sliding_law->drag(tauc(i, j), velocity(i, j).u, velocity(i, j).v);
          break;
        }
        case MASK_GROUNDED: {
          double scaling = sub_gl ? grounded_fraction(i, j) : 1.0;
          beta = scaling * m_basal_sliding_law->drag(tauc(i, j), velocity(i, j).u, velocity(i, j).v);
          break;
        }
        case MASK_ICE_FREE_OCEAN:
        default:
          beta = 0.0;
        }

        beta_u = beta;
        beta_v = beta;
      }

      {
        // Set very high basal drag *in the direction along the boundary* at locations
        // bordering "fjord walls".

        auto M   = m_mask.star_int(i, j);
        auto b   = bed.star(i, j);
        double h = surface(i, j);

        if ((ice_free(M.n) and b.n > h) or (ice_free(M.s) and b.s > h)) {
          beta_u += beta_lateral_margin;
        }

        if ((ice_free(M.e) and b.e > h) or (ice_free(M.w) and b.w > h)) {
          beta_v += beta_lateral_margin;
        }
      }

      // add beta to diagonal entries
      eq1[diag_u] += beta_u;
      eq2[diag_v] += beta_v;

      if (flow_line_mode) {
        // set values corresponding to a trivial equation v = 0
        for (int k = 0; k < n_nonzeros; ++k) {
          eq2[k] = 0.0;
        }
        eq2[diag_v] = m_scaling;
      }

      // check diagonal entries:
      const double eps = 1e-16;
      if (fabs(eq1[diag_u]) < eps) {
        if (replace_zero_diagonal_entries) {
          eq1[diag_u] = beta_ice_free_bedrock;
        } else {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                        "first  (X) equation in the SSAFD system:"
                                        " zero diagonal entry at a regular (not Dirichlet B.C.)"
                                        " location: i = %d, j = %d\n",
                                        i, j);
        }
      }
      if (fabs(eq2[diag_v]) < eps) {
        if (replace_zero_diagonal_entries) {
          eq2[diag_v] = beta_ice_free_bedrock;
        } else {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                        "second (Y) equation in the SSAFD system:"
                                        " zero diagonal entry at a regular (not Dirichlet B.C.)"
                                        " location: i = %d, j = %d\n",
                                        i, j);
        }
      }

      row.i = i;
      row.j = j;
      for (int m = 0; m < n_nonzeros; m++) {
        col[m].i = I[m];
        col[m].j = J[m];
        col[m].c = C[m];
      }

      // set coefficients of the first equation:
      row.c = 0;
      ierr  = MatSetValuesStencil(A, 1, &row, n_nonzeros, col, eq1, INSERT_VALUES);
      PISM_CHK(ierr, "MatSetValuesStencil");

      // set coefficients of the second equation:
      row.c = 1;
      ierr  = MatSetValuesStencil(A, 1, &row, n_nonzeros, col, eq2, INSERT_VALUES);
      PISM_CHK(ierr, "MatSetValuesStencil");
    } // i,j-loop
  } catch (...) {
    loop.failed();
  }
  loop.check();

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  PISM_CHK(ierr, "MatAssemblyBegin");

  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  PISM_CHK(ierr, "MatAssemblyEnd");
#if (Pism_DEBUG == 1)
  ierr = MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  PISM_CHK(ierr, "MatSetOption");
#endif
}

//! \brief Compute the vertically-averaged horizontal velocity from the shallow
//! shelf approximation.
/*!
This is the main procedure in the SSAFD.  It manages the nonlinear solve process
and the Picard iteration.

The outer loop (over index `k`) is the nonlinear iteration.  In this loop the effective
viscosity is computed by compute_nuH_staggered() and then the linear system is
set up and solved.

Specifically, we call the PETSc procedure KSPSolve() to solve the linear system.
Solving the linear system is also a loop, an iteration, but it occurs
inside KSPSolve().  The user has full control of the KSP solve through the PETSc
interface.  The default choicess for KSP type `-ksp_type` and preconditioner type
`-pc_type` are GMRES(30) for the former and block Jacobi with ILU on the
blocks for the latter.  The defaults usually work.  These choices are important
but poorly understood.  The eigenvalues of the linearized
SSA vary with ice sheet geometry and temperature in ways that are not well-studied.
Nonetheless these eigenvalues determine the convergence of
this (inner) linear iteration.  A well-chosen preconditioner can put the eigenvalues
in the right place so that the KSP can converge quickly.

The preconditioner will behave differently on different numbers of
processors.  If the user wants the results of SSA calculations to be
independent of the number of processors, then `-pc_type none` could
be used, but performance will be poor.

If you want to test different KSP methods, it may be helpful to see how many
iterations were necessary.  Use `-ksp_monitor`.
Initial testing implies that CGS takes roughly half the iterations of
GMRES(30), but is not significantly faster because the iterations are each
roughly twice as slow.  The outputs of PETSc options `-ksp_monitor_singular_value`,
`-ksp_compute_eigenvalues` and `-ksp_plot_eigenvalues -draw_pause N`
(the last holds plots for N seconds) may be useful to diagnose.

The outer loop terminates when the effective viscosity times thickness
is no longer changing much, according to the tolerance set by the
option `-ssafd_picard_rtol`. The outer loop also terminates when a maximum
number of iterations is exceeded. We save the velocity from the last
time step in order to have a better estimate of the effective
viscosity than the u=v=0 result.

In truth there is an "outer outer" loop (over index `l`).  This attempts to
over-regularize the effective viscosity if the nonlinear iteration (the "outer"
loop over `k`) is not converging with the default regularization.  The same
over-regularization is attempted if the KSP object reports that it has not
converged.

(An alternative for recovery in the KSP diverged case, suggested by Jed, is to
revert to a direct linear solve, either for the whole domain (not scalable) or
on the subdomains.  This recovery alternative requires a more nontrivial choice
but it may be worthwhile.  Note the user can already do `-pc_type asm
-sub_pc_type lu` at the command line, forcing subdomain direct solves.)

FIXME: update this doxygen comment
*/
void SSAFD::solve(const Inputs &inputs) {

  // Store away old SSA velocity (it might be needed in case a solver
  // fails).
  m_velocity_old.copy_from(m_velocity);

  // These computations do not depend on the solution, so they need to
  // be done only once.
  {
    // update the cell type mask using the ice-free thickness threshold for stress balance
    // computations
    {
      const double H_threshold = m_config->get_number("stress_balance.ice_free_thickness_standard");
      GeometryCalculator gc(*m_config);
      gc.set_icefree_thickness(H_threshold);

      gc.compute_mask(inputs.geometry->sea_level_elevation, inputs.geometry->bed_elevation,
                      inputs.geometry->ice_thickness, //
                      m_mask);                        // output
    }
    compute_driving_stress(inputs.geometry->ice_thickness, inputs.geometry->ice_surface_elevation,
                           m_mask, inputs.no_model_mask, //
                           m_taud);                      // output
    assemble_rhs(inputs, m_mask, //
                 m_rhs);         // output
    compute_hardav_staggered(inputs, m_hardness);
  }

  for (unsigned int k = 0; k < 3; ++k) {
    try {
      if (k == 0) {
        // default strategy
        picard_iteration(inputs, m_config->get_number("stress_balance.ssa.epsilon"), 1.0);

        break;
      }
      if (k == 1) {
        // try underrelaxing the iteration
        const double underrelax =
            m_config->get_number("stress_balance.ssa.fd.nuH_iter_failure_underrelaxation");
        m_log->message(
            1, "  re-trying with effective viscosity under-relaxation (parameter = %.2f) ...\n",
            underrelax);
        picard_iteration(inputs, m_config->get_number("stress_balance.ssa.epsilon"), underrelax);

        break;
      }
      if (k == 2) {
        // try over-regularization
        picard_strategy_regularization(inputs);

        break;
      }

      // if we reached this, then all strategies above failed
      write_system_petsc("all_strategies_failed");
      throw RuntimeError(PISM_ERROR_LOCATION, "all SSAFD strategies failed");
    } catch (PicardFailure &f) {
      // proceed to the next strategy
    }
  }

  if (m_config->get_flag("stress_balance.ssa.fd.extrapolate_at_margins")) {
    extrapolate_velocity(inputs.geometry->cell_type, m_velocity);
  }

  // Post-process velocities if the user asked for it:
  if (m_config->get_flag("stress_balance.ssa.fd.brutal_sliding")) {
    const double brutal_sliding_scaleFactor =
        m_config->get_number("stress_balance.ssa.fd.brutal_sliding_scale");
    m_velocity.scale(brutal_sliding_scaleFactor);

    m_velocity.update_ghosts();
  }
}

void SSAFD::picard_iteration(const Inputs &inputs, double nuH_regularization,
                             double nuH_iter_failure_underrelax) {

  if (m_default_pc_failure_count < m_default_pc_failure_max_count) {
    // Give BJACOBI another shot if we haven't tried it enough yet

    try {
      pc_setup_bjacobi();
      picard_manager(inputs, nuH_regularization, nuH_iter_failure_underrelax);

    } catch (KSPFailure &f) {

      m_default_pc_failure_count += 1;

      m_log->message(1, "  re-trying using the Additive Schwarz preconditioner...\n");

      pc_setup_asm();

      m_velocity.copy_from(m_velocity_old);

      picard_manager(inputs, nuH_regularization, nuH_iter_failure_underrelax);
    }

  } else {
    // otherwise use ASM
    pc_setup_asm();

    picard_manager(inputs, nuH_regularization, nuH_iter_failure_underrelax);
  }
}

//! \brief Manages the Picard iteration loop.
void SSAFD::picard_manager(const Inputs &inputs, double nuH_regularization,
                           double nuH_iter_failure_underrelax) {
  PetscErrorCode ierr;
  double nuH_norm, nuH_norm_change;
  // ksp_iterations should be a PetscInt because it is used in the
  // KSPGetIterationNumber() call below
  PetscInt ksp_iterations, ksp_iterations_total = 0, outer_iterations;
  KSPConvergedReason reason;

  int max_iterations =
      static_cast<int>(m_config->get_number("stress_balance.ssa.fd.max_iterations"));
  double ssa_relative_tolerance =
      m_config->get_number("stress_balance.ssa.fd.relative_convergence");
  bool verbose = m_log->get_threshold() >= 2, very_verbose = m_log->get_threshold() > 2;

  // set the initial guess:
  m_velocity_global.copy_from(m_velocity);

  m_stdout_ssa.clear();

  bool use_cfbc = m_config->get_flag("stress_balance.calving_front_stress_bc");

  if (use_cfbc) {
    compute_nuH_staggered_cfbc(inputs.geometry->ice_thickness, m_mask, m_velocity, m_hardness,
                               nuH_regularization, m_nuH);
  } else {
    compute_nuH_staggered(inputs.geometry->ice_thickness, m_velocity, m_hardness,
                          nuH_regularization, m_nuH);
  }
  update_nuH_viewers();

  // outer loop
  for (int k = 0; k < max_iterations; ++k) {

    if (very_verbose) {
      m_stdout_ssa += pism::printf("  %2d:", k);
    }

    // in preparation of measuring change of effective viscosity:
    m_nuH_old.copy_from(m_nuH);

    // assemble (or re-assemble) matrix, which depends on updated viscosity
    assemble_matrix(inputs, m_velocity, true, m_A);

    {
      array::Vector residual(m_grid, "ssa_residual");
      residual.copy_from(m_rhs);
      residual.scale(-1.0);

      ierr = MatMultAdd(m_A, m_velocity_global.vec(), residual.vec(), residual.vec());
      PISM_CHK(ierr, "MatMultAdd");

      auto filename = pism::printf("ssa_residual_%d.nc", k);
      residual.dump(filename.c_str());

      m_rhs.dump("ssa_rhs.nc");
    }

    if (very_verbose) {
      m_stdout_ssa += "A:";
    }

    // Call PETSc to solve linear system by iterative method; "inner iteration":
    ierr = KSPSetOperators(m_KSP, m_A, m_A);
    PISM_CHK(ierr, "KSPSetOperator");

    ierr = KSPSolve(m_KSP, m_rhs.vec(), m_velocity_global.vec());
    PISM_CHK(ierr, "KSPSolve");

    // Check if diverged; report to standard out about iteration
    ierr = KSPGetConvergedReason(m_KSP, &reason);
    PISM_CHK(ierr, "KSPGetConvergedReason");

    if (reason < 0) {
      // KSP diverged
      m_log->message(1, "PISM WARNING:  KSPSolve() reports 'diverged'; reason = %d = '%s'\n",
                     reason, KSPConvergedReasons[reason]);

      write_system_petsc("kspdivergederror");

      // Tell the caller that we failed. (The caller might try again,
      // though.)
      throw KSPFailure(KSPConvergedReasons[reason]);
    }

    // report on KSP success; the "inner" iteration is done
    ierr = KSPGetIterationNumber(m_KSP, &ksp_iterations);
    PISM_CHK(ierr, "KSPGetIterationNumber");

    ksp_iterations_total += ksp_iterations;

    if (very_verbose) {
      m_stdout_ssa += pism::printf("S:%d,%d: ", (int)ksp_iterations, reason);
    }

    // limit ice speed
    {
      auto max_speed = m_config->get_number("stress_balance.ssa.fd.max_speed", "m second-1");
      int high_speed_counter = 0;

      array::AccessScope list{ &m_velocity_global };

      for (auto p = m_grid->points(); p; p.next()) {
        const int i = p.i(), j = p.j();

        auto speed = m_velocity_global(i, j).magnitude();

        if (speed > max_speed) {
          m_velocity_global(i, j) *= max_speed / speed;
          high_speed_counter += 1;
        }
      }

      high_speed_counter = GlobalSum(m_grid->com, high_speed_counter);

      if (high_speed_counter > 0) {
        m_log->message(2, "  SSA speed was capped at %d locations\n", high_speed_counter);
      }
    }

    // Communicate so that we have stencil width for evaluation of effective
    // viscosity on next "outer" iteration (and geometry etc. if done):
    // Note that copy_from() updates ghosts of m_velocity.
    m_velocity.copy_from(m_velocity_global);

    // update viscosity and check for viscosity convergence
    if (use_cfbc) {
      compute_nuH_staggered_cfbc(inputs.geometry->ice_thickness, m_mask, m_velocity, m_hardness,
                                 nuH_regularization, m_nuH);
    } else {
      compute_nuH_staggered(inputs.geometry->ice_thickness, m_velocity, m_hardness,
                            nuH_regularization, m_nuH);
    }

    if (nuH_iter_failure_underrelax != 1.0) {
      m_nuH.scale(nuH_iter_failure_underrelax);
      m_nuH.add(1.0 - nuH_iter_failure_underrelax, m_nuH_old);
    }
    compute_nuH_norm(nuH_norm, nuH_norm_change);

    update_nuH_viewers();

    if (very_verbose) {
      m_stdout_ssa += pism::printf("|nu|_2, |Delta nu|_2/|nu|_2 = %10.3e %10.3e\n", nuH_norm,
                                   nuH_norm_change / nuH_norm);

      // assume that high verbosity shows interest in immediate
      // feedback about SSA iterations
      m_log->message(2, m_stdout_ssa);

      m_stdout_ssa.clear();
    }

    outer_iterations = k + 1;

    if (nuH_norm == 0 || nuH_norm_change / nuH_norm < ssa_relative_tolerance) {
      goto done;
    }

  } // outer loop (k)

  // If we're here, it means that we exceeded max_iterations and still
  // failed.

  throw PicardFailure(pism::printf("effective viscosity not converged after %d iterations\n"
                                   "with nuH_regularization=%8.2e.",
                                   max_iterations, nuH_regularization));

done:

  if (very_verbose) {
    auto tempstr =
        pism::printf("... =%5d outer iterations, ~%3.1f KSP iterations each\n",
                     (int)outer_iterations, ((double)ksp_iterations_total) / outer_iterations);
    m_stdout_ssa += tempstr;
  } else if (verbose) {
    // at default verbosity, just record last nuH_norm_change and iterations
    auto tempstr =
        pism::printf("%5d outer iterations, ~%3.1f KSP iterations each\n", (int)outer_iterations,
                     ((double)ksp_iterations_total) / outer_iterations);

    m_stdout_ssa += tempstr;
  }

  if (verbose) {
    m_stdout_ssa = "  SSA: " + m_stdout_ssa;
  }
}

//! Old SSAFD recovery strategy: increase the SSA regularization parameter.
void SSAFD::picard_strategy_regularization(const Inputs &inputs) {
  // this has no units; epsilon goes up by this ratio when previous value failed
  const double DEFAULT_EPSILON_MULTIPLIER_SSA = 4.0;
  double nuH_regularization                   = m_config->get_number("stress_balance.ssa.epsilon");
  unsigned int k = 0, max_tries = 5;

  if (nuH_regularization <= 0.0) {
    throw PicardFailure("will not attempt over-regularization of nuH\n"
                        "because nuH_regularization == 0.0.");
  }

  while (k < max_tries) {
    m_velocity.copy_from(m_velocity_old);
    m_log->message(1, "  re-trying with nuH_regularization multiplied by %8.2f...\n",
                   DEFAULT_EPSILON_MULTIPLIER_SSA);

    nuH_regularization *= DEFAULT_EPSILON_MULTIPLIER_SSA;

    try {
      // 1.0 is the under-relaxation parameter
      picard_iteration(inputs, nuH_regularization, 1.0);
      // if this call succeeded, stop over-regularizing
      break;
    } catch (PicardFailure &f) {
      k += 1;

      if (k == max_tries) {
        throw PicardFailure("exceeded the max. number of nuH over-regularization attempts");
      }
    }
  }
}

//! \brief Compute the norm of nu H and the change in nu H.
/*!
Verification and PST experiments
suggest that an \f$L^1\f$ criterion for convergence is best.  For verification
there seems to be little difference, presumably because the solutions are smooth
and the norms are roughly equivalent on a subspace of smooth functions.  For PST,
the \f$L^1\f$ criterion gives faster runs with essentially the same results.
Presumably that is because rapid (temporal and spatial) variation in
\f$\nu H\f$ occurs at margins, occupying very few horizontal grid cells.
For the significant (e.g.~in terms of flux) parts of the flow, it is o.k. to ignore
a bit of bad behavior at these few places, and \f$L^1\f$ ignores it more than
\f$L^2\f$ (much less \f$L^\infty\f$, which might not work at all).
 */
void SSAFD::compute_nuH_norm(double &norm, double &norm_change) {

  const double area      = m_grid->cell_area();
  const NormType MY_NORM = NORM_1;

  // Test for change in nu
  m_nuH_old.add(-1, m_nuH);

  std::vector<double> nuNorm = m_nuH.norm(MY_NORM), nuChange = m_nuH_old.norm(MY_NORM);

  nuChange[0] *= area;
  nuChange[1] *= area;
  nuNorm[0] *= area;
  nuNorm[1] *= area;

  norm_change = sqrt(PetscSqr(nuChange[0]) + PetscSqr(nuChange[1]));
  norm        = sqrt(PetscSqr(nuNorm[0]) + PetscSqr(nuNorm[1]));
}

//! \brief Computes vertically-averaged ice hardness on the staggered grid.
void SSAFD::compute_hardav_staggered(const Inputs &inputs,
                                     array::Staggered &result) {
  const array::Scalar &thickness = inputs.geometry->ice_thickness;

  const array::Array3D &enthalpy = *inputs.enthalpy;

  const double *E_ij = NULL, *E_offset = NULL;

  auto Mz = m_grid->Mz();
  std::vector<double> E(Mz);

  array::AccessScope list{ &thickness, &enthalpy, &result, &m_mask };

  ParallelSection loop(m_grid->com);
  try {
    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      E_ij = enthalpy.get_column(i, j);
      for (int o = 0; o < 2; o++) {
        const int oi = 1 - o, oj = o;
        double H;

        if (m_mask.icy(i, j) && m_mask.icy(i + oi, j + oj)) {
          H = 0.5 * (thickness(i, j) + thickness(i + oi, j + oj));
        } else if (m_mask.icy(i, j)) {
          H = thickness(i, j);
        } else {
          H = thickness(i + oi, j + oj);
        }

        if (H == 0) {
          result(i, j, o) = -1e6; // an obviously impossible value
          continue;
        }

        E_offset = enthalpy.get_column(i + oi, j + oj);
        // build a column of enthalpy values a the current location:
        for (unsigned int k = 0; k < Mz; ++k) {
          E[k] = 0.5 * (E_ij[k] + E_offset[k]);
        }

        result(i, j, o) = rheology::averaged_hardness(*m_flow_law, H, m_grid->kBelowHeight(H),
                                                      m_grid->z().data(), E.data());
      } // o
    }   // loop over points
  } catch (...) {
    loop.failed();
  }
  loop.check();

  fracture_induced_softening(inputs.fracture_density);
}


/*! @brief Correct vertically-averaged hardness using a
    parameterization of the fracture-induced softening.

  See T. Albrecht, A. Levermann; Fracture-induced softening for
  large-scale ice dynamics; (2013), The Cryosphere Discussions 7;
  4501-4544; DOI:10.5194/tcd-7-4501-2013

  Note that this paper proposes an adjustment of the enhancement factor:

  \f[E_{\text{effective}} = E \cdot (1 - (1-\epsilon) \phi)^{-n}.\f]

  Let \f$E_{\text{effective}} = E\cdot C\f$, where \f$C\f$ is the
  factor defined by the formula above.

  Recall that the effective viscosity is defined by

  \f[\nu(D) = \frac12 B D^{(1-n)/(2n)}\f]

  and the viscosity form of the flow law is

  \f[\sigma'_{ij} = E_{\text{effective}}^{-\frac1n}2\nu(D) D_{ij}.\f]

  Then

  \f[\sigma'_{ij} = E_{\text{effective}}^{-\frac1n}BD^{(1-n)/(2n)}D_{ij}.\f]

  Using the fact that \f$E_{\text{effective}} = E\cdot C\f$, this can be rewritten as

  \f[\sigma'_{ij} = E^{-\frac1n} \left(C^{-\frac1n}B\right) D^{(1-n)/(2n)}D_{ij}.\f]

  So scaling the enhancement factor by \f$C\f$ is equivalent to scaling
  ice hardness \f$B\f$ by \f$C^{-\frac1n}\f$.
*/
void SSAFD::fracture_induced_softening(const array::Scalar *fracture_density) {
  if (fracture_density == nullptr) {
    return;
  }

  const double epsilon = m_config->get_number("fracture_density.softening_lower_limit"),
               n_glen  = m_flow_law->exponent();

  array::AccessScope list{ &m_hardness, fracture_density };

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    for (int o = 0; o < 2; o++) {
      const int oi = 1 - o, oj = o;

      const double
          // fracture density on the staggered grid:
          phi = 0.5 * ((*fracture_density)(i, j) + (*fracture_density)(i + oi, j + oj)),
          // the line below implements equation (6) in the paper
          softening = pow((1.0 - (1.0 - epsilon) * phi), -n_glen);

      m_hardness(i, j, o) *= pow(softening, -1.0 / n_glen);
    }
  }
}

//! \brief Compute the product of ice thickness and effective viscosity (on the
//! staggered grid).
/*!
In PISM the product \f$\nu H\f$ can be
  - constant, or
  - can be computed with a constant ice hardness \f$\bar B\f$ (temperature-independent)
    but with dependence of the viscosity on the strain rates, or
  - it can depend on the strain rates \e and have a vertically-averaged ice
    hardness depending on temperature or enthalpy.

The flow law in ice stream and ice shelf regions must, for now, be a
(temperature-dependent) Glen law.  This is the only flow law we know how to
invert to ``viscosity form''.  More general forms like Goldsby-Kohlstedt are
not yet inverted.

The viscosity form of a Glen law is
   \f[ \nu(T^*,D) = \frac{1}{2} B(T^*) D^{(1/n)-1}\, D_{ij} \f]
where
   \f[  D_{ij} = \frac{1}{2} \left(\frac{\partial U_i}{\partial x_j} +
                                   \frac{\partial U_j}{\partial x_i}\right) \f]
is the strain rate tensor and \f$B\f$ is an ice hardness related to
the ice softness \f$A(T^*)\f$ by
   \f[ B(T^*)=A(T^*)^{-1/n}  \f]
in the case of a temperature dependent Glen-type law.  (Here \f$T^*\f$ is the
pressure-adjusted temperature.)

The effective viscosity is then
   \f[ \nu = \frac{\bar B}{2} \left[\left(\frac{\partial u}{\partial x}\right)^2 +
                               \left(\frac{\partial v}{\partial y}\right)^2 +
                               \frac{\partial u}{\partial x} \frac{\partial v}{\partial y} +
                               \frac{1}{4} \left(\frac{\partial u}{\partial y}
                                                 + \frac{\partial v}{\partial x}\right)^2
                               \right]^{(1-n)/(2n)}  \f]
where in the temperature-dependent case
   \f[ \bar B = \frac{1}{H}\,\int_b^h B(T^*)\,dz\f]
This integral is approximately computed by the trapezoid rule.

The result of this procedure is \f$\nu H\f$, not just \f$\nu\f$, this it is
a vertical integral, not average, of viscosity.

The resulting effective viscosity times thickness is regularized by ensuring that
its minimum is at least \f$\epsilon\f$.  This regularization constant is an argument.

In this implementation we set \f$\nu H\f$ to a constant anywhere the ice is
thinner than a certain minimum. See SSAStrengthExtension and compare how this
issue is handled when -cfbc is set.
*/
void SSAFD::compute_nuH_staggered(const array::Scalar1 &ice_thickness,
                                  const array::Vector1 &velocity, const array::Staggered &hardness,
                                  double nuH_regularization, array::Staggered &result) {

  const array::Vector &uv = velocity; // shortcut

  array::AccessScope list{ &result, &uv, &hardness, &ice_thickness };

  double n_glen                 = m_flow_law->exponent(),
         nu_enhancement_scaling = 1.0 / pow(m_e_factor, 1.0 / n_glen);

  const double dx = m_grid->dx(), dy = m_grid->dy();

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    auto V = velocity.box(i, j);

    for (int o = 0; o < 2; ++o) {
      const int oi = 1 - o, oj = o;

      const double H = 0.5 * (ice_thickness(i, j) + ice_thickness(i + oi, j + oj));

      if (H < strength_extension->get_min_thickness()) {
        result(i, j, o) = strength_extension->get_notional_strength();
        continue;
      }

      double u_x, u_y, v_x, v_y;
      // Check the offset to determine how to differentiate velocity
      if (o == 0) {
        u_x = (V.e.u - V.c.u) / dx;
        v_x = (V.e.v - V.c.v) / dx;
        u_y = (V.n.u + V.ne.u - V.s.u - V.se.u) / (4 * dy);
        v_y = (V.n.v + V.ne.v - V.s.v - V.se.v) / (4 * dy);
      } else {
        u_x = (V.e.u + V.ne.u - V.w.u - V.nw.u) / (4 * dx);
        v_x = (V.e.v + V.ne.v - V.w.v - V.nw.v) / (4 * dx);
        u_y = (V.n.u - V.c.u) / dy;
        v_y = (V.n.v - V.c.v) / dy;
      }

      double nu = 0.0;
      m_flow_law->effective_viscosity(hardness(i, j, o),
                                      secondInvariant_2D({ u_x, v_x }, { u_y, v_y }), &nu, nullptr);

      result(i, j, o) = nu * H;

      // include the SSA enhancement factor; in most cases m_e_factor is 1
      result(i, j, o) *= nu_enhancement_scaling;

      // We ensure that nuH is bounded below by a positive constant.
      result(i, j, o) += nuH_regularization;

    } // o-loop
  }   // i,j-loop

  // Some communication
  result.update_ghosts();
}

/**
 * @brief Compute the product of ice viscosity and thickness on the
 * staggered grid. Used when CFBC is enabled.
 *
 * @param[out] result nu*H product
 * @param[in] nuH_regularization regularization parameter (added to nu*H to keep it away from zero)
 *
 * @return 0 on success
 */
void SSAFD::compute_nuH_staggered_cfbc(const array::Scalar1 &ice_thickness,
                                       const array::CellType2 &cell_type,
                                       const array::Vector2 &velocity,
                                       const array::Staggered &hardness, double nuH_regularization,
                                       array::Staggered &result) {

  double n_glen                 = m_flow_law->exponent(),
         nu_enhancement_scaling = 1.0 / pow(m_e_factor, 1.0 / n_glen);

  const double dx = m_grid->dx(), dy = m_grid->dy();

  array::AccessScope list{ &cell_type, &m_work, &velocity };

  assert(m_work.stencil_width() >= 1);

  for (auto p = m_grid->points(1); p; p.next()) {
    const int i = p.i(), j = p.j();

    // x-derivative, i-offset
    {
      if (cell_type.icy(i, j) && cell_type.icy(i + 1, j)) {
        m_work(i, j).u_x = (velocity(i + 1, j).u - velocity(i, j).u) / dx; // u_x
        m_work(i, j).v_x = (velocity(i + 1, j).v - velocity(i, j).v) / dx; // v_x
        m_work(i, j).w_i = 1.0;
      } else {
        m_work(i, j).u_x = 0.0;
        m_work(i, j).v_x = 0.0;
        m_work(i, j).w_i = 0.0;
      }
    }

    // y-derivative, j-offset
    {
      if (cell_type.icy(i, j) && cell_type.icy(i, j + 1)) {
        m_work(i, j).u_y = (velocity(i, j + 1).u - velocity(i, j).u) / dy; // u_y
        m_work(i, j).v_y = (velocity(i, j + 1).v - velocity(i, j).v) / dy; // v_y
        m_work(i, j).w_j = 1.0;
      } else {
        m_work(i, j).u_y = 0.0;
        m_work(i, j).v_y = 0.0;
        m_work(i, j).w_j = 0.0;
      }
    }
  }

  list.add({ &result, &hardness, &ice_thickness });

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double u_x, u_y, v_x, v_y, H, nu, W;
    // i-offset
    {
      if (cell_type.icy(i, j) && cell_type.icy(i + 1, j)) {
        H = 0.5 * (ice_thickness(i, j) + ice_thickness(i + 1, j));
      } else if (cell_type.icy(i, j)) {
        H = ice_thickness(i, j);
      } else {
        H = ice_thickness(i + 1, j);
      }

      if (H >= strength_extension->get_min_thickness()) {
        u_x = m_work(i, j).u_x;
        v_x = m_work(i, j).v_x;

        W = m_work(i, j).w_j + m_work(i, j - 1).w_j + m_work(i + 1, j - 1).w_j +
            m_work(i + 1, j).w_j;
        if (W > 0) {
          u_y = 1.0 / W *
                (m_work(i, j).u_y + m_work(i, j - 1).u_y + m_work(i + 1, j - 1).u_y +
                 m_work(i + 1, j).u_y);
          v_y = 1.0 / W *
                (m_work(i, j).v_y + m_work(i, j - 1).v_y + m_work(i + 1, j - 1).v_y +
                 m_work(i + 1, j).v_y);
        } else {
          u_y = 0.0;
          v_y = 0.0;
        }

        m_flow_law->effective_viscosity(
            hardness(i, j, 0), secondInvariant_2D({ u_x, v_x }, { u_y, v_y }), &nu, nullptr);
        result(i, j, 0) = nu * H;
      } else {
        result(i, j, 0) = strength_extension->get_notional_strength();
      }
    }

    // j-offset
    {
      if (cell_type.icy(i, j) && cell_type.icy(i, j + 1)) {
        H = 0.5 * (ice_thickness(i, j) + ice_thickness(i, j + 1));
      } else if (cell_type.icy(i, j)) {
        H = ice_thickness(i, j);
      } else {
        H = ice_thickness(i, j + 1);
      }

      if (H >= strength_extension->get_min_thickness()) {
        u_y = m_work(i, j).u_y;
        v_y = m_work(i, j).v_y;

        W = m_work(i, j).w_i + m_work(i - 1, j).w_i + m_work(i - 1, j + 1).w_i +
            m_work(i, j + 1).w_i;
        if (W > 0.0) {
          u_x = 1.0 / W *
                (m_work(i, j).u_x + m_work(i - 1, j).u_x + m_work(i - 1, j + 1).u_x +
                 m_work(i, j + 1).u_x);
          v_x = 1.0 / W *
                (m_work(i, j).v_x + m_work(i - 1, j).v_x + m_work(i - 1, j + 1).v_x +
                 m_work(i, j + 1).v_x);
        } else {
          u_x = 0.0;
          v_x = 0.0;
        }

        m_flow_law->effective_viscosity(hardness(i, j, 1),
                                        secondInvariant_2D({ u_x, v_x }, { u_y, v_y }), &nu, NULL);
        result(i, j, 1) = nu * H;
      } else {
        result(i, j, 1) = strength_extension->get_notional_strength();
      }
    }

    // adjustments:
    for (int o = 0; o < 2; ++o) {
      // include the SSA enhancement factor; in most cases ssa_enhancement_factor is 1
      result(i, j, o) *= nu_enhancement_scaling;

      // We ensure that nuH is bounded below by a positive constant.
      result(i, j, o) += nuH_regularization;
    }
  }

  // Some communication
  result.update_ghosts();
}

//! Update the nuH viewer, which shows log10(nu H).
void SSAFD::update_nuH_viewers() {

  if (not m_view_nuh) {
    return;
  }

  array::Scalar tmp(m_grid, "nuH");
  tmp.metadata(0)
      .long_name("log10 of (viscosity * thickness)")
      .units("Pa s m");

  array::AccessScope list{&m_nuH, &tmp};

  for (auto p = m_grid->points(); p; p.next()) {
    const int i = p.i(), j = p.j();

    double avg_nuH = 0.5 * (m_nuH(i,j,0) + m_nuH(i,j,1));
    if (avg_nuH > 1.0e14) {
      tmp(i,j) = log10(avg_nuH);
    } else {
      tmp(i,j) = 14.0;
    }
  }

  if (not m_nuh_viewer) {
    m_nuh_viewer.reset(new petsc::Viewer(m_grid->com, "nuH", m_nuh_viewer_size,
                                         m_grid->Lx(), m_grid->Ly()));
  }

  tmp.view({m_nuh_viewer});
}

//! \brief Checks if a cell is near or at the ice front.
/*!
 * You need to create array::AccessScope object and add mask to it.
 *
 * Note that a cell is a CFBC location of one of four direct neighbors is ice-free.
 *
 * If one of the diagonal neighbors is ice-free we don't use the CFBC, but we
 * do need to compute weights used in the SSA discretization (see
 * assemble_matrix()) to avoid differencing across interfaces between icy
 * and ice-free cells.
 *
 * This method ensures that checks in assemble_rhs() and assemble_matrix() are
 * consistent.
 */
bool SSAFD::is_marginal(int i, int j, bool ssa_dirichlet_bc) {

  auto M = m_mask.box_int(i, j);

  using mask::ice_free;
  using mask::ice_free_ocean;
  using mask::icy;

  if (ssa_dirichlet_bc) {
    return icy(M.c) &&
      (ice_free(M.e) || ice_free(M.w) || ice_free(M.n) || ice_free(M.s) ||
       ice_free(M.ne) || ice_free(M.se) || ice_free(M.nw) || ice_free(M.sw));
  }

  return icy(M.c) &&
    (ice_free_ocean(M.e) || ice_free_ocean(M.w) ||
     ice_free_ocean(M.n) || ice_free_ocean(M.s) ||
     ice_free_ocean(M.ne) || ice_free_ocean(M.se) ||
     ice_free_ocean(M.nw) || ice_free_ocean(M.sw));
}

void SSAFD::write_system_petsc(const std::string &namepart) {
  PetscErrorCode ierr;

  // write a file with a fixed filename; avoid zillions of files
  std::string filename = "SSAFD_" + namepart + ".petsc";
  m_log->message(1,
             "  writing linear system to PETSc binary file %s ...\n", filename.c_str());

  petsc::Viewer viewer;       // will be destroyed automatically
  ierr = PetscViewerBinaryOpen(m_grid->com, filename.c_str(), FILE_MODE_WRITE,
                               viewer.rawptr());
  PISM_CHK(ierr, "PetscViewerBinaryOpen");

  ierr = MatView(m_A, viewer);
  PISM_CHK(ierr, "MatView");

  ierr = VecView(m_rhs.vec(), viewer);
  PISM_CHK(ierr, "VecView");
}

SSAFD_nuH::SSAFD_nuH(const SSAFD *m) : Diag<SSAFD>(m) {
  m_vars = { { m_sys, "nuH[0]" }, { m_sys, "nuH[1]" } };
  m_vars[0]
      .long_name("ice thickness times effective viscosity, i-offset")
      .units("Pa s m")
      .output_units("kPa s m");
  m_vars[1]
      .long_name("ice thickness times effective viscosity, j-offset")
      .units("Pa s m")
      .output_units("kPa s m");
}

std::shared_ptr<array::Array> SSAFD_nuH::compute_impl() const {
  auto result = allocate<array::Staggered>("nuH");

  result->copy_from(model->integrated_viscosity());

  return result;
}


//! @brief Computes the driving shear stress at the base of ice
//! (diagnostically).
/*! This is *not* a duplicate of SSB_taud: SSAFD_taud::compute() uses
  SSAFD::compute_driving_stress(), which tries to be smarter near ice margins.
*/
class SSAFD_taud : public Diag<SSAFD> {
public:
  SSAFD_taud(const SSAFD *m) : Diag<SSAFD>(m) {

    // set metadata:
    m_vars = { { m_sys, "taud_x" }, { m_sys, "taud_y" } };

    m_vars[0].long_name("X-component of the driving shear stress at the base of ice");
    m_vars[1].long_name("Y-component of the driving shear stress at the base of ice");

    for (auto &v : m_vars) {
      v.units("Pa");
      v["comment"] = "this is the driving stress used by the SSAFD solver";
    }
  }

protected:
  std::shared_ptr<array::Array> compute_impl() const {
    auto result = allocate<array::Vector>("taud");

    result->copy_from(model->driving_stress());

    return result;
  }
};


//! @brief Computes the magnitude of the driving shear stress at the base of
//! ice (diagnostically).
class SSAFD_taud_mag : public Diag<SSAFD> {
public:
  SSAFD_taud_mag(const SSAFD *m) : Diag<SSAFD>(m) {

    // set metadata:
    m_vars = { { m_sys, "taud_mag" } };

    m_vars[0].long_name("magnitude of the driving shear stress at the base of ice").units("Pa");
    m_vars[0]["comment"] = "this is the magnitude of the driving stress used by the SSAFD solver";
  }

protected:
  virtual std::shared_ptr<array::Array> compute_impl() const {
    auto result = allocate<array::Scalar>("taud_mag");

    compute_magnitude(model->driving_stress(), *result);

    return result;
  }
};

DiagnosticList SSAFD::diagnostics_impl() const {
  DiagnosticList result = ShallowStressBalance::diagnostics_impl();

  // replace these diagnostics
  result["taud"] = Diagnostic::Ptr(new SSAFD_taud(this));
  result["taud_mag"] = Diagnostic::Ptr(new SSAFD_taud_mag(this));
  result["nuH"] = Diagnostic::Ptr(new SSAFD_nuH(this));

  return result;
}

const array::Staggered & SSAFD::integrated_viscosity() const {
  return m_nuH;
}


} // end of namespace stressbalance
} // end of namespace pism
