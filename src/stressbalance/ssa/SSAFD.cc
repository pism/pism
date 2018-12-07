// Copyright (C) 2004--2018 Constantine Khroulev, Ed Bueler and Jed Brown
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

#include "SSAFD.hh"
#include "SSAFD_diagnostics.hh"
#include "pism/util/Mask.hh"
#include "pism/basalstrength/basal_resistance.hh"
#include "pism/util/pism_options.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/util/Vars.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Time.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace stressbalance {

using namespace pism::mask;

SSAFD::KSPFailure::KSPFailure(const char* reason)
  : RuntimeError(ErrorLocation(), std::string("SSAFD KSP (linear solver) failed: ") + reason){
  // empty
}

SSAFD::PicardFailure::PicardFailure(const std::string &message)
  : RuntimeError(ErrorLocation(), "SSAFD Picard iterations failed: " + message) {
  // empty
}

SSA* SSAFDFactory(IceGrid::ConstPtr g) {
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
SSAFD::SSAFD(IceGrid::ConstPtr g)
  : SSA(g) {
  m_b.create(m_grid, "right_hand_side", WITHOUT_GHOSTS);

  m_velocity_old.create(m_grid, "velocity_old", WITH_GHOSTS);
  m_velocity_old.set_attrs("internal",
                           "old SSA velocity field; used for re-trying with a different epsilon",
                           "m s-1", "");

  m_hardness.create(m_grid, "hardness", WITHOUT_GHOSTS);
  m_hardness.set_attrs("diagnostic",
                       "vertically-averaged ice hardness",
                       pism::printf("Pa s%f", 1.0 / m_flow_law->exponent()),
                       "");

  m_nuH.create(m_grid, "nuH", WITH_GHOSTS);
  m_nuH.set_attrs("internal",
                  "ice thickness times effective viscosity",
                  "Pa s m", "");

  m_nuH_old.create(m_grid, "nuH_old", WITH_GHOSTS);
  m_nuH_old.set_attrs("internal",
                      "ice thickness times effective viscosity (before an update)",
                      "Pa s m", "");

  m_work.create(m_grid, "m_work", WITH_GHOSTS,
                2, /* stencil width */
                6  /* dof */);
  m_work.set_attrs("internal",
                   "temporary storage used to compute nuH",
                   "", "");

  m_scaling = 1.0e9;  // comparable to typical beta for an ice stream;

  // The nuH viewer:
  m_view_nuh = false;
  m_nuh_viewer_size = 300;

  // PETSc objects and settings
  {
    PetscErrorCode ierr;
    ierr = DMSetMatType(*m_da, MATAIJ);
    PISM_CHK(ierr, "DMSetMatType");

    ierr = DMCreateMatrix(*m_da, m_A.rawptr());
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

SSAFD::~SSAFD() {
  // empty
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

  m_log->message(2,
             "  [using the KSP-based finite difference implementation]\n");

  // options
  options::Integer viewer_size("-ssa_nuh_viewer_size", "nuH viewer size",
                               m_nuh_viewer_size);
  m_nuh_viewer_size = viewer_size;
  m_view_nuh = options::Bool("-ssa_view_nuh", "Enable the SSAFD nuH runtime viewer");

  if (m_config->get_boolean("stress_balance.calving_front_stress_bc")) {
    m_log->message(2,
               "  using PISM-PIK calving-front stress boundary condition ...\n");
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
void SSAFD::assemble_rhs(const Inputs &inputs) {
  const IceModelVec2S
    &thickness             = inputs.geometry->ice_thickness,
    &bed                   = inputs.geometry->bed_elevation,
    &sea_level             = inputs.geometry->sea_level_elevation,
    *melange_back_pressure = inputs.melange_back_pressure;

  const double
    dx                        = m_grid->dx(),
    dy                        = m_grid->dy(),
    ice_free_default_velocity = 0.0,
    standard_gravity          = m_config->get_double("constants.standard_gravity"),
    rho_ocean                 = m_config->get_double("constants.sea_water.density"),
    rho_ice                   = m_config->get_double("constants.ice.density");

  const bool
    use_cfbc          = m_config->get_boolean("stress_balance.calving_front_stress_bc"),
    is_dry_simulation = m_config->get_boolean("ocean.always_grounded");

  // FIXME: bedrock_boundary is a misleading name
  bool bedrock_boundary = m_config->get_boolean("stress_balance.ssa.dirichlet_bc");

  compute_driving_stress(*inputs.geometry, m_taud);

  IceModelVec::AccessList list{&m_taud, &m_b};

  if (inputs.bc_values && inputs.bc_mask) {
    list.add({inputs.bc_values, inputs.bc_mask});
  }

  if (use_cfbc) {
    list.add({&thickness, &bed, &m_mask, &sea_level});
  }

  if (use_cfbc && melange_back_pressure != NULL) {
    list.add(*melange_back_pressure);
  }

  m_b.set(0.0);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (inputs.bc_values != NULL &&
        inputs.bc_mask->as_int(i, j) == 1) {
      m_b(i, j).u = m_scaling * (*inputs.bc_values)(i, j).u;
      m_b(i, j).v = m_scaling * (*inputs.bc_values)(i, j).v;
      continue;
    }

    if (use_cfbc) {
      double H_ij = thickness(i,j);

      int
        M_ij = m_mask.as_int(i,j),
        M_e  = m_mask.as_int(i + 1,j),
        M_w  = m_mask.as_int(i - 1,j),
        M_n  = m_mask.as_int(i,j + 1),
        M_s  = m_mask.as_int(i,j - 1);

      // Note: this sets velocities at both ice-free ocean and ice-free
      // bedrock to zero. This means that we need to set boundary conditions
      // at both ice/ice-free-ocean and ice/ice-free-bedrock interfaces below
      // to be consistent.
      if (ice_free(M_ij)) {
        m_b(i, j).u = m_scaling * ice_free_default_velocity;
        m_b(i, j).v = m_scaling * ice_free_default_velocity;
        continue;
      }

      if (is_marginal(i, j, bedrock_boundary)) {
        int aMM = 1, aPP = 1, bMM = 1, bPP = 1;
        // direct neighbors
        if (bedrock_boundary) {
          if (ice_free_ocean(M_e))
            aPP = 0;
          if (ice_free_ocean(M_w))
            aMM = 0;
          if (ice_free_ocean(M_n))
            bPP = 0;
          if (ice_free_ocean(M_s))
            bMM = 0;
        } else {
          if (ice_free(M_e))
            aPP = 0;
          if (ice_free(M_w))
            aMM = 0;
          if (ice_free(M_n))
            bPP = 0;
          if (ice_free(M_s))
            bMM = 0;
        }

        double ocean_pressure = ocean_pressure_difference(ocean(M_ij), is_dry_simulation,
                                                          H_ij, bed(i,j), sea_level(i, j),
                                                          rho_ice, rho_ocean, standard_gravity);

        if (melange_back_pressure != NULL) {
          double lambda = (*melange_back_pressure)(i, j);

          // adjust the "pressure difference term" using the provided
          // "melange back pressure fraction".
          ocean_pressure *= (1.0 - lambda);
        }

        // Note that if the current cell is "marginal" but not a CFBC
        // location, the following two lines are equaivalent to the "usual
        // case" below.
        m_b(i, j).u = m_taud(i,j).u + (aMM - aPP) * ocean_pressure / dx;
        m_b(i, j).v = m_taud(i,j).v + (bMM - bPP) * ocean_pressure / dy;

        continue;
      } // end of "if (is_marginal(i, j))"

        // If we reached this point, then CFBC are enabled, but we are in the
        // interior of a sheet or shelf. See "usual case" below.

    }   // end of "if (use_cfbc)"

    // usual case: use already computed driving stress
    m_b(i, j).u = m_taud(i, j).u;
    m_b(i, j).v = m_taud(i, j).v;
  }
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
void SSAFD::assemble_matrix(const Inputs &inputs,
                            bool include_basal_shear, Mat A) {
  PetscErrorCode ierr = 0;

  // shortcut:
  const IceModelVec2V &vel = m_velocity;

  const IceModelVec2S
    &thickness         = inputs.geometry->ice_thickness,
    &bed               = inputs.geometry->bed_elevation,
    &surface           = inputs.geometry->ice_surface_elevation,
    &grounded_fraction = inputs.geometry->cell_grounded_fraction,
    &tauc              = *inputs.basal_yield_stress;

  const double
    dx                    = m_grid->dx(),
    dy                    = m_grid->dy(),
    beta_ice_free_bedrock = m_config->get_double("basal_resistance.beta_ice_free_bedrock");

  const bool use_cfbc     = m_config->get_boolean("stress_balance.calving_front_stress_bc");
  const bool replace_zero_diagonal_entries =
    m_config->get_boolean("stress_balance.ssa.fd.replace_zero_diagonal_entries");

  // FIXME: bedrock_boundary is a misleading name
  const bool bedrock_boundary = m_config->get_boolean("stress_balance.ssa.dirichlet_bc");

  ierr = MatZeroEntries(A);
  PISM_CHK(ierr, "MatZeroEntries");

  IceModelVec::AccessList list{&m_nuH, &tauc, &vel, &m_mask};

  if (inputs.bc_values && inputs.bc_mask) {
    list.add(*inputs.bc_mask);
  }

  const bool sub_gl = m_config->get_boolean("geometry.grounded_cell_fraction");
  if (sub_gl) {
    list.add(grounded_fraction);
  }

  // handles friction of the ice cell along ice-free bedrock margins when bedrock higher than ice
  // surface (in simplified setups)
  bool lateral_drag_enabled=m_config->get_boolean("stress_balance.ssa.fd.lateral_drag.enabled");
  if (lateral_drag_enabled) {
    list.add({&thickness, &bed, &surface});
  }
  double lateral_drag_viscosity=m_config->get_double("stress_balance.ssa.fd.lateral_drag.viscosity");
  double HminFrozen=0.0;

  /* matrix assembly loop */
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // Handle the easy case: provided Dirichlet boundary conditions
      if (inputs.bc_values && inputs.bc_mask && inputs.bc_mask->as_int(i,j) == 1) {
        // set diagonal entry to one (scaled); RHS entry will be known velocity;
        set_diagonal_matrix_entry(A, i, j, m_scaling);
        continue;
      }

      /* Provide shorthand for the following staggered coefficients  nu H:
       *      c_n
       *  c_w     c_e
       *      c_s
       */
      // const
      double c_w = m_nuH(i-1,j,0);
      double c_e = m_nuH(i,j,0);
      double c_s = m_nuH(i,j-1,1);
      double c_n = m_nuH(i,j,1);

      if (lateral_drag_enabled) {
        // if option is set, the viscosity at ice-bedrock boundary layer will
        // be prescribed and is a temperature-independent free (user determined) parameter

        // direct neighbors
        int
          M_e = m_mask.as_int(i + 1,j),
          M_w = m_mask.as_int(i - 1,j),
          M_n = m_mask.as_int(i,j + 1),
          M_s = m_mask.as_int(i,j - 1);

        if (thickness(i,j) > HminFrozen) {
          if (bed(i-1,j) > surface(i,j) && ice_free_land(M_w)) {
            c_w = lateral_drag_viscosity * 0.5 * (thickness(i,j)+thickness(i-1,j));
          }
          if (bed(i+1,j) > surface(i,j) && ice_free_land(M_e)) {
            c_e = lateral_drag_viscosity * 0.5 * (thickness(i,j)+thickness(i+1,j));
          }
          if (bed(i,j+1) > surface(i,j) && ice_free_land(M_n)) {
            c_n = lateral_drag_viscosity * 0.5 * (thickness(i,j)+thickness(i,j+1));
          }
          if (bed(i,j-1) > surface(i,j) && ice_free_land(M_s)) {
            c_s = lateral_drag_viscosity * 0.5 * (thickness(i,j)+thickness(i+1,j));
          }
        }
      }

      // We use DAGetMatrix to obtain the SSA matrix, which means that all 18
      // non-zeros get allocated, even though we use only 13 (or 14). The
      // remaining 5 (or 4) coefficients are zeros, but we set them anyway,
      // because this makes the code easier to understand.
      const int sten = 18;
      MatStencil row, col[sten];

      int aMn = 1, aPn = 1, aMM = 1, aPP = 1, aMs = 1, aPs = 1;
      int bPw = 1, bPP = 1, bPe = 1, bMw = 1, bMM = 1, bMe = 1;

      int M_ij = m_mask.as_int(i,j);

      if (use_cfbc) {
        int
          // direct neighbors
          M_e = m_mask.as_int(i + 1,j),
          M_w = m_mask.as_int(i - 1,j),
          M_n = m_mask.as_int(i,j + 1),
          M_s = m_mask.as_int(i,j - 1),
          // "diagonal" neighbors
          M_ne = m_mask.as_int(i + 1,j + 1),
          M_se = m_mask.as_int(i + 1,j - 1),
          M_nw = m_mask.as_int(i - 1,j + 1),
          M_sw = m_mask.as_int(i - 1,j - 1);

        // Note: this sets velocities at both ice-free ocean and ice-free
        // bedrock to zero. This means that we need to set boundary conditions
        // at both ice/ice-free-ocean and ice/ice-free-bedrock interfaces below
        // to be consistent.
        if (ice_free(M_ij)) {
          set_diagonal_matrix_entry(A, i, j, m_scaling);
          continue;
        }

        if (is_marginal(i, j, bedrock_boundary)) {
          // If at least one of the following four conditions is "true", we're
          // at a CFBC location.
          if (bedrock_boundary) {

            if (ice_free_ocean(M_e))
              aPP = 0;
            if (ice_free_ocean(M_w))
              aMM = 0;
            if (ice_free_ocean(M_n))
              bPP = 0;
            if (ice_free_ocean(M_s))
              bMM = 0;

            // decide whether to use centered or one-sided differences
            if (ice_free_ocean(M_n) || ice_free_ocean(M_ne))
              aPn = 0;
            if (ice_free_ocean(M_e) || ice_free_ocean(M_ne))
              bPe = 0;
            if (ice_free_ocean(M_e) || ice_free_ocean(M_se))
              bMe = 0;
            if (ice_free_ocean(M_s) || ice_free_ocean(M_se))
              aPs = 0;
            if (ice_free_ocean(M_s) || ice_free_ocean(M_sw))
              aMs = 0;
            if (ice_free_ocean(M_w) || ice_free_ocean(M_sw))
              bMw = 0;
            if (ice_free_ocean(M_w) || ice_free_ocean(M_nw))
              bPw = 0;
            if (ice_free_ocean(M_n) || ice_free_ocean(M_nw))
              aMn = 0;

          } else {                // if (not bedrock_boundary)

            if (ice_free(M_e))
              aPP = 0;
            if (ice_free(M_w))
              aMM = 0;
            if (ice_free(M_n))
              bPP = 0;
            if (ice_free(M_s))
              bMM = 0;

            // decide whether to use centered or one-sided differences
            if (ice_free(M_n) || ice_free(M_ne))
              aPn = 0;
            if (ice_free(M_e) || ice_free(M_ne))
              bPe = 0;
            if (ice_free(M_e) || ice_free(M_se))
              bMe = 0;
            if (ice_free(M_s) || ice_free(M_se))
              aPs = 0;
            if (ice_free(M_s) || ice_free(M_sw))
              aMs = 0;
            if (ice_free(M_w) || ice_free(M_sw))
              bMw = 0;
            if (ice_free(M_w) || ice_free(M_nw))
              bPw = 0;
            if (ice_free(M_n) || ice_free(M_nw))
              aMn = 0;

          } // end of the else clause following "if (bedrock_boundary)"
        }   // end of "if (is_marginal(i, j, bedrock_boundary))"
      }     // end of "if (use_cfbc)"

      /* begin Maxima-generated code */
      const double dx2 = dx*dx, dy2 = dy*dy, d4 = 4*dx*dy, d2 = 2*dx*dy;

      /* Coefficients of the discretization of the first equation; u first, then v. */
      double eq1[] = {
        0,  -c_n*bPP/dy2,  0,
        -4*c_w*aMM/dx2,  (c_n*bPP+c_s*bMM)/dy2+(4*c_e*aPP+4*c_w*aMM)/dx2,  -4*c_e*aPP/dx2,
        0,  -c_s*bMM/dy2,  0,
        c_w*aMM*bPw/d2+c_n*aMn*bPP/d4,  (c_n*aPn*bPP-c_n*aMn*bPP)/d4+(c_w*aMM*bPP-c_e*aPP*bPP)/d2,  -c_e*aPP*bPe/d2-c_n*aPn*bPP/d4,
        (c_w*aMM*bMw-c_w*aMM*bPw)/d2+(c_n*aMM*bPP-c_s*aMM*bMM)/d4,  (c_n*aPP*bPP-c_n*aMM*bPP-c_s*aPP*bMM+c_s*aMM*bMM)/d4+(c_e*aPP*bPP-c_w*aMM*bPP-c_e*aPP*bMM+c_w*aMM*bMM)/d2,  (c_e*aPP*bPe-c_e*aPP*bMe)/d2+(c_s*aPP*bMM-c_n*aPP*bPP)/d4,
        -c_w*aMM*bMw/d2-c_s*aMs*bMM/d4,  (c_s*aMs*bMM-c_s*aPs*bMM)/d4+(c_e*aPP*bMM-c_w*aMM*bMM)/d2,  c_e*aPP*bMe/d2+c_s*aPs*bMM/d4,
      };

      /* Coefficients of the discretization of the second equation; u first, then v. */
      double eq2[] = {
        c_w*aMM*bPw/d4+c_n*aMn*bPP/d2,  (c_n*aPn*bPP-c_n*aMn*bPP)/d2+(c_w*aMM*bPP-c_e*aPP*bPP)/d4,  -c_e*aPP*bPe/d4-c_n*aPn*bPP/d2,
        (c_w*aMM*bMw-c_w*aMM*bPw)/d4+(c_n*aMM*bPP-c_s*aMM*bMM)/d2,  (c_n*aPP*bPP-c_n*aMM*bPP-c_s*aPP*bMM+c_s*aMM*bMM)/d2+(c_e*aPP*bPP-c_w*aMM*bPP-c_e*aPP*bMM+c_w*aMM*bMM)/d4,  (c_e*aPP*bPe-c_e*aPP*bMe)/d4+(c_s*aPP*bMM-c_n*aPP*bPP)/d2,
        -c_w*aMM*bMw/d4-c_s*aMs*bMM/d2,  (c_s*aMs*bMM-c_s*aPs*bMM)/d2+(c_e*aPP*bMM-c_w*aMM*bMM)/d4,  c_e*aPP*bMe/d4+c_s*aPs*bMM/d2,
        0,  -4*c_n*bPP/dy2,  0,
        -c_w*aMM/dx2,  (4*c_n*bPP+4*c_s*bMM)/dy2+(c_e*aPP+c_w*aMM)/dx2,  -c_e*aPP/dx2,
        0,  -4*c_s*bMM/dy2,  0,
      };

      /* i indices */
      const int I[] = {
        i-1,  i,  i+1,
        i-1,  i,  i+1,
        i-1,  i,  i+1,
        i-1,  i,  i+1,
        i-1,  i,  i+1,
        i-1,  i,  i+1,
      };

      /* j indices */
      const int J[] = {
        j+1,  j+1,  j+1,
        j,  j,  j,
        j-1,  j-1,  j-1,
        j+1,  j+1,  j+1,
        j,  j,  j,
        j-1,  j-1,  j-1,
      };

      /* component indices */
      const int C[] = {
        0,  0,  0,
        0,  0,  0,
        0,  0,  0,
        1,  1,  1,
        1,  1,  1,
        1,  1,  1,
      };
      /* end Maxima-generated code */

      /* Dragging ice experiences friction at the bed determined by the
       *    IceBasalResistancePlasticLaw::drag() methods.  These may be a plastic,
       *    pseudo-plastic, or linear friction law.  Dragging is done implicitly
       *    (i.e. on left side of SSA eqns).  */
      double beta = 0.0;
      if (include_basal_shear) {
        if (grounded_ice(M_ij)) {
          beta = m_basal_sliding_law->drag(tauc(i,j), vel(i,j).u, vel(i,j).v);
        } else if (ice_free_land(M_ij)) {
          // apply drag even in this case, to help with margins; note ice free
          // areas already have a strength extension
          beta = beta_ice_free_bedrock;
        }
        if (sub_gl) {
          // reduce the basal drag at grid cells that are partially grounded:
          if (icy(M_ij)) {
            beta = grounded_fraction(i,j) * m_basal_sliding_law->drag(tauc(i,j), vel(i,j).u, vel(i,j).v);
          }
        }
      }

      // add beta to diagonal entries
      eq1[4]  += beta;
      eq2[13] += beta;

      // check diagonal entries:
      const double eps = 1e-16;
      if (fabs(eq1[4]) < eps) {
        if (replace_zero_diagonal_entries) {
          eq1[4] = beta_ice_free_bedrock;
        } else {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION, "first  (X) equation in the SSAFD system:"
                                        " zero diagonal entry at a regular (not Dirichlet B.C.)"
                                        " location: i = %d, j = %d\n", i, j);
        }
      }
      if (fabs(eq2[13]) < eps) {
        if (replace_zero_diagonal_entries) {
          eq2[13] = beta_ice_free_bedrock;
        } else {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION, "second (Y) equation in the SSAFD system:"
                                        " zero diagonal entry at a regular (not Dirichlet B.C.)"
                                        " location: i = %d, j = %d\n", i, j);
        }
      }

      row.i = i;
      row.j = j;
      for (int m = 0; m < sten; m++) {
        col[m].i = I[m];
        col[m].j = J[m];
        col[m].c = C[m];
      }

      // set coefficients of the first equation:
      row.c = 0;
      ierr = MatSetValuesStencil(A, 1, &row, sten, col, eq1, INSERT_VALUES);
      PISM_CHK(ierr, "MatSetValuesStencil");

      // set coefficients of the second equation:
      row.c = 1;
      ierr = MatSetValuesStencil(A, 1, &row, sten, col, eq2, INSERT_VALUES);
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
#if (PISM_DEBUG==1)
  ierr = MatSetOption(A,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);
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
option `-ssa_rtol`. The outer loop also terminates when a maximum
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
  // be done once.
  {
    assemble_rhs(inputs);
    compute_hardav_staggered(inputs);
  }

  for (unsigned int k = 0; k < 3; ++k) {
    try {
      if (k == 0) {
        // default strategy
        picard_iteration(inputs, m_config->get_double("stress_balance.ssa.epsilon"), 1.0);

        break;
      } else if (k == 1) {
        // try underrelaxing the iteration
        const double underrelax = m_config->get_double("stress_balance.ssa.fd.nuH_iter_failure_underrelaxation");
        m_log->message(1,
                   "  re-trying with effective viscosity under-relaxation (parameter = %.2f) ...\n",
                   underrelax);
        picard_iteration(inputs, m_config->get_double("stress_balance.ssa.epsilon"), underrelax);

        break;
      } else if (k == 2) {
        // try over-regularization
        picard_strategy_regularization(inputs);

        break;
      } else {
        // if we reached this, then all strategies above failed
        write_system_petsc("all_strategies_failed");
        throw RuntimeError(PISM_ERROR_LOCATION, "all SSAFD strategies failed");
      }
    } catch (PicardFailure &f) {
      // proceed to the next strategy
    }
  }

  // Post-process velocities if the user asked for it:
  if (m_config->get_boolean("stress_balance.ssa.fd.brutal_sliding")) {
    const double brutal_sliding_scaleFactor = m_config->get_double("stress_balance.ssa.fd.brutal_sliding_scale");
    m_velocity.scale(brutal_sliding_scaleFactor);

    m_velocity.update_ghosts();
  }
}

void SSAFD::picard_iteration(const Inputs &inputs,
                             double nuH_regularization,
                             double nuH_iter_failure_underrelax) {

  if (m_default_pc_failure_count < m_default_pc_failure_max_count) {
    // Give BJACOBI another shot if we haven't tried it enough yet

    try {
      pc_setup_bjacobi();
      picard_manager(inputs, nuH_regularization,
                     nuH_iter_failure_underrelax);

    } catch (KSPFailure &f) {

      m_default_pc_failure_count += 1;

      m_log->message(1,
                 "  re-trying using the Additive Schwarz preconditioner...\n");

      pc_setup_asm();

      m_velocity.copy_from(m_velocity_old);

      picard_manager(inputs, nuH_regularization,
                     nuH_iter_failure_underrelax);
    }

  } else {
    // otherwise use ASM
    pc_setup_asm();

    picard_manager(inputs, nuH_regularization,
                   nuH_iter_failure_underrelax);
  }
}

//! \brief Manages the Picard iteration loop.
void SSAFD::picard_manager(const Inputs &inputs,
                           double nuH_regularization,
                           double nuH_iter_failure_underrelax) {
  PetscErrorCode ierr;
  double   nuH_norm, nuH_norm_change;
  // ksp_iterations should be a PetscInt because it is used in the
  // KSPGetIterationNumber() call below
  PetscInt    ksp_iterations, ksp_iterations_total = 0, outer_iterations;
  KSPConvergedReason  reason;

  unsigned int max_iterations = static_cast<int>(m_config->get_double("stress_balance.ssa.fd.max_iterations"));
  double ssa_relative_tolerance = m_config->get_double("stress_balance.ssa.fd.relative_convergence");
  char tempstr[100] = "";
  bool verbose = m_log->get_threshold() >= 2,
    very_verbose = m_log->get_threshold() > 2;

  // set the initial guess:
  m_velocity_global.copy_from(m_velocity);

  m_stdout_ssa.clear();

  bool use_cfbc = m_config->get_boolean("stress_balance.calving_front_stress_bc");

  if (use_cfbc == true) {
    compute_nuH_staggered_cfbc(*inputs.geometry, nuH_regularization, m_nuH);
  } else {
    compute_nuH_staggered(*inputs.geometry, nuH_regularization, m_nuH);
  }
  update_nuH_viewers();

  // outer loop
  for (unsigned int k = 0; k < max_iterations; ++k) {

    if (very_verbose) {
      snprintf(tempstr, 100, "  %2d:", k);
      m_stdout_ssa += tempstr;
    }

    // in preparation of measuring change of effective viscosity:
    m_nuH_old.copy_from(m_nuH);

    // assemble (or re-assemble) matrix, which depends on updated viscosity
    assemble_matrix(inputs, true, m_A);

    if (very_verbose) {

      m_stdout_ssa += "A:";
    }

    // Call PETSc to solve linear system by iterative method; "inner iteration":
    ierr = KSPSetOperators(m_KSP, m_A, m_A);
    PISM_CHK(ierr, "KSPSetOperator");

    ierr = KSPSolve(m_KSP, m_b.vec(), m_velocity_global.vec());
    PISM_CHK(ierr, "KSPSolve");

    // Check if diverged; report to standard out about iteration
    ierr = KSPGetConvergedReason(m_KSP, &reason);
    PISM_CHK(ierr, "KSPGetConvergedReason");

    if (reason < 0) {
      // KSP diverged
      m_log->message(1,
                 "PISM WARNING:  KSPSolve() reports 'diverged'; reason = %d = '%s'\n",
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
      snprintf(tempstr, 100, "S:%d,%d: ", (int)ksp_iterations, reason);
      m_stdout_ssa += tempstr;
    }

    // limit ice speed
    {
      auto max_speed = m_config->get_double("stress_balance.ssa.fd.max_speed", "m second-1");
      int high_speed_counter = 0;

      IceModelVec::AccessList list{&m_velocity_global};

      for (Points p(*m_grid); p; p.next()) {
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
    if (use_cfbc == true) {
      compute_nuH_staggered_cfbc(*inputs.geometry, nuH_regularization, m_nuH);
    } else {
      compute_nuH_staggered(*inputs.geometry, nuH_regularization, m_nuH);
    }

    if (nuH_iter_failure_underrelax != 1.0) {
      m_nuH.scale(nuH_iter_failure_underrelax);
      m_nuH.add(1.0 - nuH_iter_failure_underrelax, m_nuH_old);
    }
    compute_nuH_norm(nuH_norm, nuH_norm_change);

    update_nuH_viewers();

    if (very_verbose) {
      snprintf(tempstr, 100, "|nu|_2, |Delta nu|_2/|nu|_2 = %10.3e %10.3e\n",
               nuH_norm, nuH_norm_change/nuH_norm);

      m_stdout_ssa += tempstr;

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
    snprintf(tempstr, 100, "... =%5d outer iterations, ~%3.1f KSP iterations each\n",
             (int)outer_iterations, ((double) ksp_iterations_total) / outer_iterations);

    m_stdout_ssa += tempstr;
  } else if (verbose) {
    // at default verbosity, just record last nuH_norm_change and iterations
    snprintf(tempstr, 100, "%5d outer iterations, ~%3.1f KSP iterations each\n",
             (int)outer_iterations, ((double) ksp_iterations_total) / outer_iterations);

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
  double nuH_regularization = m_config->get_double("stress_balance.ssa.epsilon");
  unsigned int k = 0, max_tries = 5;

  if (nuH_regularization <= 0.0) {
    throw PicardFailure("will not attempt over-regularization of nuH\n"
                        "because nuH_regularization == 0.0.");
  }

  while (k < max_tries) {
    m_velocity.copy_from(m_velocity_old);
    m_log->message(1,
               "  re-trying with nuH_regularization multiplied by %8.2f...\n",
               DEFAULT_EPSILON_MULTIPLIER_SSA);

    nuH_regularization *= DEFAULT_EPSILON_MULTIPLIER_SSA;

    try {
      // 1.0 is the under-relaxation parameter
      picard_iteration(inputs, nuH_regularization, 1.0);
      // if this call succeeded, stop over-regularizing
      break;
    }
    catch (PicardFailure &f) {
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

  const double area = m_grid->cell_area();
  const NormType MY_NORM = NORM_1;

  // Test for change in nu
  m_nuH_old.add(-1, m_nuH);

  std::vector<double>
    nuNorm   = m_nuH.norm_all(MY_NORM),
    nuChange = m_nuH_old.norm_all(MY_NORM);

  nuChange[0] *= area;
  nuChange[1] *= area;
  nuNorm[0]   *= area;
  nuNorm[1]   *= area;

  norm_change = sqrt(PetscSqr(nuChange[0]) + PetscSqr(nuChange[1]));
  norm = sqrt(PetscSqr(nuNorm[0]) + PetscSqr(nuNorm[1]));
}

//! \brief Computes vertically-averaged ice hardness on the staggered grid.
void SSAFD::compute_hardav_staggered(const Inputs &inputs) {
  const IceModelVec2S
    &thickness = inputs.geometry->ice_thickness;

  const IceModelVec3 &enthalpy = *inputs.enthalpy;

  const double
    *E_ij     = NULL,
    *E_offset = NULL;

  std::vector<double> E(m_grid->Mz());

  IceModelVec::AccessList list{&thickness, &enthalpy, &m_hardness, &m_mask};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      E_ij = enthalpy.get_column(i,j);
      for (int o=0; o<2; o++) {
        const int oi = 1-o, oj=o;
        double H;

        if (m_mask.icy(i,j) && m_mask.icy(i+oi,j+oj)) {
          H = 0.5 * (thickness(i,j) + thickness(i+oi,j+oj));
        } else if (m_mask.icy(i,j)) {
          H = thickness(i,j);
        }  else {
          H = thickness(i+oi,j+oj);
        }

        if (H == 0) {
          m_hardness(i,j,o) = -1e6; // an obviously impossible value
          continue;
        }

        E_offset = enthalpy.get_column(i+oi,j+oj);
        // build a column of enthalpy values a the current location:
        for (unsigned int k = 0; k < m_grid->Mz(); ++k) {
          E[k] = 0.5 * (E_ij[k] + E_offset[k]);
        }

        m_hardness(i,j,o) = rheology::averaged_hardness(*m_flow_law,
                                                        H, m_grid->kBelowHeight(H),
                                                        &(m_grid->z()[0]), &E[0]);
      } // o
    } // loop over points
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
void SSAFD::fracture_induced_softening(const IceModelVec2S *fracture_density) {
  if (not fracture_density) {
    return;
  }

  const double
    epsilon = m_config->get_double("fracture_density.softening_lower_limit"),
    n_glen  = m_flow_law->exponent();

  IceModelVec::AccessList list{&m_hardness, fracture_density};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    for (int o=0; o<2; o++) {
      const int oi = 1-o, oj=o;

      const double
        // fracture density on the staggered grid:
        phi       = 0.5 * ((*fracture_density)(i,j) + (*fracture_density)(i+oi,j+oj)),
        // the line below implements equation (6) in the paper
        softening = pow((1.0-(1.0-epsilon)*phi), -n_glen);

      m_hardness(i,j,o) *= pow(softening,-1.0/n_glen);
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
void SSAFD::compute_nuH_staggered(const Geometry &geometry,
                                  double nuH_regularization,
                                  IceModelVec2Stag &result) {

  const IceModelVec2V &uv = m_velocity; // shortcut

  IceModelVec::AccessList list{&result, &uv, &m_hardness, &geometry.ice_thickness};

  double ssa_enhancement_factor = m_flow_law->enhancement_factor(),
    n_glen = m_flow_law->exponent(),
    nu_enhancement_scaling = 1.0 / pow(ssa_enhancement_factor, 1.0/n_glen);

  const double dx = m_grid->dx(), dy = m_grid->dy();

  for (int o=0; o<2; ++o) {
    const int oi = 1 - o, oj=o;
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double H = 0.5 * (geometry.ice_thickness(i,j) + geometry.ice_thickness(i+oi,j+oj));

      if (H < strength_extension->get_min_thickness()) {
        result(i,j,o) = strength_extension->get_notional_strength();
        continue;
      }

      double u_x, u_y, v_x, v_y;
      // Check the offset to determine how to differentiate velocity
      if (o == 0) {
        u_x = (uv(i+1,j).u - uv(i,j).u) / dx;
        u_y = (uv(i,j+1).u + uv(i+1,j+1).u - uv(i,j-1).u - uv(i+1,j-1).u) / (4*dy);
        v_x = (uv(i+1,j).v - uv(i,j).v) / dx;
        v_y = (uv(i,j+1).v + uv(i+1,j+1).v - uv(i,j-1).v - uv(i+1,j-1).v) / (4*dy);
      } else {
        u_x = (uv(i+1,j).u + uv(i+1,j+1).u - uv(i-1,j).u - uv(i-1,j+1).u) / (4*dx);
        u_y = (uv(i,j+1).u - uv(i,j).u) / dy;
        v_x = (uv(i+1,j).v + uv(i+1,j+1).v - uv(i-1,j).v - uv(i-1,j+1).v) / (4*dx);
        v_y = (uv(i,j+1).v - uv(i,j).v) / dy;
      }

      double nu = 0.0;
      m_flow_law->effective_viscosity(m_hardness(i,j,o),
                                      secondInvariant_2D(Vector2(u_x, v_x),
                                                         Vector2(u_y, v_y)),
                                      &nu, NULL);

      result(i,j,o) = nu * H;

      // include the SSA enhancement factor; in most cases ssa_enhancement_factor is 1
      result(i,j,o) *= nu_enhancement_scaling;

      // We ensure that nuH is bounded below by a positive constant.
      result(i,j,o) += nuH_regularization;

    } // i,j-loop
  } // o-loop


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
 * m_work storage scheme:
 *
 * m_work(i,j,0) - u_x on the i-offset
 * m_work(i,j,1) - v_x on the i-offset
 * m_work(i,j,2) - i-offset weight
 * m_work(i,j,3) - u_y on the j-offset
 * m_work(i,j,4) - v_y on the j-offset
 * m_work(i,j,5) - j-offset weight
 *
 * @return 0 on success
 */
void SSAFD::compute_nuH_staggered_cfbc(const Geometry &geometry,
                                       double nuH_regularization,
                                       IceModelVec2Stag &result) {

  const IceModelVec2S &thickness = geometry.ice_thickness;

  const IceModelVec2V &uv = m_velocity; // shortcut
  double ssa_enhancement_factor = m_flow_law->enhancement_factor(),
    n_glen = m_flow_law->exponent(),
    nu_enhancement_scaling = 1.0 / pow(ssa_enhancement_factor, 1.0/n_glen);

  const unsigned int U_X = 0, V_X = 1, W_I = 2, U_Y = 3, V_Y = 4, W_J = 5;

  const double dx = m_grid->dx(), dy = m_grid->dy();

  IceModelVec::AccessList list{&m_mask, &m_work, &m_velocity};

  assert(m_velocity.stencil_width() >= 2);
  assert(m_mask.stencil_width()    >= 2);
  assert(m_work.stencil_width()     >= 1);

  for (PointsWithGhosts p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // x-derivative, i-offset
    {
      if (m_mask.icy(i,j) && m_mask.icy(i+1,j)) {
        m_work(i,j,U_X) = (uv(i+1,j).u - uv(i,j).u) / dx; // u_x
        m_work(i,j,V_X) = (uv(i+1,j).v - uv(i,j).v) / dx; // v_x
        m_work(i,j,W_I) = 1.0;
      } else {
        m_work(i,j,U_X) = 0.0;
        m_work(i,j,V_X) = 0.0;
        m_work(i,j,W_I) = 0.0;
      }
    }

    // y-derivative, j-offset
    {
      if (m_mask.icy(i,j) && m_mask.icy(i,j+1)) {
        m_work(i,j,U_Y) = (uv(i,j+1).u - uv(i,j).u) / dy; // u_y
        m_work(i,j,V_Y) = (uv(i,j+1).v - uv(i,j).v) / dy; // v_y
        m_work(i,j,W_J) = 1.0;
      } else {
        m_work(i,j,U_Y) = 0.0;
        m_work(i,j,V_Y) = 0.0;
        m_work(i,j,W_J) = 0.0;
      }
    }
  }

  list.add({&result, &m_hardness, &thickness});

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double u_x, u_y, v_x, v_y, H, nu, W;
    // i-offset
    {
      if (m_mask.icy(i,j) && m_mask.icy(i+1,j)) {
        H = 0.5 * (thickness(i,j) + thickness(i+1,j));
      }
      else if (m_mask.icy(i,j)) {
        H = thickness(i,j);
      } else {
        H = thickness(i+1,j);
      }

      if (H >= strength_extension->get_min_thickness()) {
        u_x = m_work(i,j,U_X);
        v_x = m_work(i,j,V_X);

        W = m_work(i,j,W_J) + m_work(i,j-1,W_J) + m_work(i+1,j-1,W_J) + m_work(i+1,j,W_J);
        if (W > 0) {
          u_y = 1.0/W * (m_work(i,j,U_Y) + m_work(i,j-1,U_Y) +
                         m_work(i+1,j-1,U_Y) + m_work(i+1,j,U_Y));
          v_y = 1.0/W * (m_work(i,j,V_Y) + m_work(i,j-1,V_Y) +
                         m_work(i+1,j-1,V_Y) + m_work(i+1,j,V_Y));
        } else {
          u_y = 0.0;
          v_y = 0.0;
        }

        m_flow_law->effective_viscosity(m_hardness(i,j,0),
                                        secondInvariant_2D(Vector2(u_x, v_x),
                                                           Vector2(u_y, v_y)),
                                        &nu, NULL);
        result(i,j,0) = nu * H;
      } else {
        result(i,j,0) = strength_extension->get_notional_strength();
      }
    }

    // j-offset
    {
      if (m_mask.icy(i,j) && m_mask.icy(i,j+1)) {
        H = 0.5 * (thickness(i,j) + thickness(i,j+1));
      } else if (m_mask.icy(i,j)) {
        H = thickness(i,j);
      } else {
        H = thickness(i,j+1);
      }

      if (H >= strength_extension->get_min_thickness()) {
        u_y = m_work(i,j,U_Y);
        v_y = m_work(i,j,V_Y);

        W = m_work(i,j,W_I) + m_work(i-1,j,W_I) + m_work(i-1,j+1,W_I) + m_work(i,j+1,W_I);
        if (W > 0.0) {
          u_x = 1.0/W * (m_work(i,j,U_X) + m_work(i-1,j,U_X) +
                         m_work(i-1,j+1,U_X) + m_work(i,j+1,U_X));
          v_x = 1.0/W * (m_work(i,j,V_X) + m_work(i-1,j,V_X) +
                         m_work(i-1,j+1,V_X) + m_work(i,j+1,V_X));
        } else {
          u_x = 0.0;
          v_x = 0.0;
        }

        m_flow_law->effective_viscosity(m_hardness(i,j,1),
                                        secondInvariant_2D(Vector2(u_x, v_x),
                                                           Vector2(u_y, v_y)),
                                        &nu, NULL);
        result(i,j,1) = nu * H;
      } else {
        result(i,j,1) = strength_extension->get_notional_strength();
      }
    }

    // adjustments:
    for (unsigned int o = 0; o < 2; ++o) {
      // include the SSA enhancement factor; in most cases ssa_enhancement_factor is 1
      result(i,j,o) *= nu_enhancement_scaling;

      // We ensure that nuH is bounded below by a positive constant.
      result(i,j,o) += nuH_regularization;
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

  IceModelVec2S tmp;
  tmp.create(m_grid, "nuH", WITHOUT_GHOSTS);
  tmp.set_attrs("temporary",
                "log10 of (viscosity * thickness)",
                "Pa s m", "");

  IceModelVec::AccessList list{&m_nuH, &tmp};

  for (Points p(*m_grid); p; p.next()) {
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

  tmp.view(m_nuh_viewer, petsc::Viewer::Ptr());
}

void SSAFD::set_diagonal_matrix_entry(Mat A, int i, int j,
                                      double value) {
  PetscErrorCode ierr;
  MatStencil row, col;
  row.i = i;
  row.j = j;
  col.i = i;
  col.j = j;

  row.c = 0;
  col.c = 0;

  ierr = MatSetValuesStencil(A, 1, &row, 1, &col, &value, INSERT_VALUES);
  PISM_CHK(ierr, "MatSetValuesStencil");

  row.c = 1;
  col.c = 1;

  ierr = MatSetValuesStencil(A, 1, &row, 1, &col, &value, INSERT_VALUES);
  PISM_CHK(ierr, "MatSetValuesStencil");
}

//! \brief Checks if a cell is near or at the ice front.
/*!
 * You need to create IceModelVec::AccessList object and add mask to it.
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

  StarStencil<int> M = m_mask.int_star(i, j);
  const int
    // "diagonal" neighbors
    M_ne = m_mask.as_int(i + 1,j + 1),
    M_se = m_mask.as_int(i + 1,j - 1),
    M_nw = m_mask.as_int(i - 1,j + 1),
    M_sw = m_mask.as_int(i - 1,j - 1);

  if (ssa_dirichlet_bc) {
    return icy(M.ij) &&
      (ice_free(M.e) || ice_free(M.w) || ice_free(M.n) || ice_free(M.s) ||
       ice_free(M_ne) || ice_free(M_se) || ice_free(M_nw) || ice_free(M_sw));
  } else {
    return icy(M.ij) &&
      (ice_free_ocean(M.e) || ice_free_ocean(M.w) ||
       ice_free_ocean(M.n) || ice_free_ocean(M.s) ||
       ice_free_ocean(M_ne) || ice_free_ocean(M_se) ||
       ice_free_ocean(M_nw) || ice_free_ocean(M_sw));
  }
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

  ierr = VecView(m_b.vec(), viewer);
  PISM_CHK(ierr, "VecView");
}

SSAFD_nuH::SSAFD_nuH(const SSAFD *m)
  : Diag<SSAFD>(m) {

  // set metadata:
  m_vars = {SpatialVariableMetadata(m_sys, "nuH[0]"),
            SpatialVariableMetadata(m_sys, "nuH[1]")};

  set_attrs("ice thickness times effective viscosity, i-offset", "",
            "Pa s m", "kPa s m", 0);
  set_attrs("ice thickness times effective viscosity, j-offset", "",
            "Pa s m", "kPa s m", 1);
}

IceModelVec::Ptr SSAFD_nuH::compute_impl() const {

  IceModelVec2Stag::Ptr result(new IceModelVec2Stag);
  result->create(m_grid, "nuH", WITH_GHOSTS);
  result->metadata(0) = m_vars[0];
  result->metadata(1) = m_vars[1];

  result->copy_from(model->integrated_viscosity());

  return result;
}

DiagnosticList SSAFD::diagnostics_impl() const {
  DiagnosticList result = SSA::diagnostics_impl();

  result["nuH"] = Diagnostic::Ptr(new SSAFD_nuH(this));

  return result;
}

const IceModelVec2Stag & SSAFD::integrated_viscosity() const {
  return m_nuH;
}


} // end of namespace stressbalance
} // end of namespace pism
