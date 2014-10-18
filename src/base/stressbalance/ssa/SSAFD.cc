// Copyright (C) 2004--2014 Constantine Khroulev, Ed Bueler and Jed Brown
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

#include "SSAFD.hh"
#include "Mask.hh"
#include "basal_resistance.hh"
#include "pism_options.hh"
#include "flowlaws.hh"
#include "PISMVars.hh"

#include <assert.h>
#include <stdexcept>

namespace pism {

using namespace mask;

SSA* SSAFDFactory(IceGrid &g, EnthalpyConverter &ec, const Config &c) {
  return new SSAFD(g,ec,c);
}

SSAFD::SSAFD(IceGrid &g, EnthalpyConverter &e, const Config &c)
  : SSA(g,e,c) {
  PetscErrorCode ierr = allocate_fd();
  if (ierr != 0) {
    throw std::runtime_error("SSAFD allocation failed");
  }
}

SSAFD::~SSAFD() {
  PetscErrorCode ierr = deallocate_fd();
  if (ierr != 0) {
    PetscPrintf(grid.com, "FATAL ERROR: SSAFD de-allocation failed.\n");
    abort();
  }
}


PetscErrorCode SSAFD::pc_setup_bjacobi() {
  PetscErrorCode ierr;
  PC pc;

  ierr = KSPSetType(m_KSP, KSPGMRES); CHKERRQ(ierr);
#if PETSC_VERSION_LT(3,5,0)
  ierr = KSPSetOperators(m_KSP, m_A, m_A, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
#else
  ierr = KSPSetOperators(m_KSP, m_A, m_A); CHKERRQ(ierr);
#endif

  // Get the PC from the KSP solver:
  ierr = KSPGetPC(m_KSP, &pc); CHKERRQ(ierr);

  // Set the PC type:
  ierr = PCSetType(pc, PCBJACOBI); CHKERRQ(ierr);

  // Process options:
  ierr = KSPSetFromOptions(m_KSP); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode SSAFD::pc_setup_asm() {
  PetscErrorCode ierr;
  PC pc, sub_pc;

  // Set parameters equivalent to
  // -ksp_type gmres -ksp_norm_type unpreconditioned -ksp_pc_side right -pc_type asm -sub_pc_type lu

  ierr = KSPSetType(m_KSP, KSPGMRES); CHKERRQ(ierr);
#if PETSC_VERSION_LT(3,5,0)
  ierr = KSPSetOperators(m_KSP, m_A, m_A, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
#else
  ierr = KSPSetOperators(m_KSP, m_A, m_A); CHKERRQ(ierr);
#endif
    
  // Switch to using the "unpreconditioned" norm.
  ierr = KSPSetNormType(m_KSP, KSP_NORM_UNPRECONDITIONED); CHKERRQ(ierr);

  // Switch to "right" preconditioning.
  ierr = KSPSetPCSide(m_KSP, PC_RIGHT); CHKERRQ(ierr);

  // Get the PC from the KSP solver:
  ierr = KSPGetPC(m_KSP, &pc); CHKERRQ(ierr);
  
  // Set the PC type:
  ierr = PCSetType(pc, PCASM); CHKERRQ(ierr);

  // Set the sub-KSP object to "preonly"
  KSP *sub_ksp;
  ierr = PCSetUp(pc); CHKERRQ(ierr);
  ierr = PCASMGetSubKSP(pc, NULL, NULL, &sub_ksp); CHKERRQ(ierr);

  ierr = KSPSetType(*sub_ksp, KSPPREONLY); CHKERRQ(ierr);

  // Set the PC of the sub-KSP to "LU".
  ierr = KSPGetPC(*sub_ksp, &sub_pc); CHKERRQ(ierr);

  ierr = PCSetType(sub_pc, PCLU); CHKERRQ(ierr);
    
  // Let the user override all this:
  // Process options:
  ierr = KSPSetFromOptions(m_KSP); CHKERRQ(ierr);

  return 0;
}

//! \brief Allocate objects specific to the SSAFD object.
/*!
Because the FD implementation of the SSA uses Picard iteration, a PETSc KSP
and Mat are used directly.  In particular we set up \f$A\f$
(Mat m_A) and a \f$b\f$ (= Vec m_b) and iteratively solve
linear systems
  \f[ A x = b \f]
where \f$x\f$ (= Vec SSAX).  A PETSc SNES object is never created.
 */
PetscErrorCode SSAFD::allocate_fd() {
  PetscErrorCode ierr;

  fracture_density = NULL;
  m_melange_back_pressure = NULL;

#if PETSC_VERSION_LT(3,5,0)
  ierr = DMCreateMatrix(*m_da, MATAIJ, &m_A); CHKERRQ(ierr);
#else
  ierr = DMSetMatType(*m_da, MATAIJ); CHKERRQ(ierr);
  ierr = DMCreateMatrix(*m_da, &m_A); CHKERRQ(ierr);
#endif

  ierr = KSPCreate(grid.com, &m_KSP); CHKERRQ(ierr);
  ierr = KSPSetOptionsPrefix(m_KSP, "ssafd_"); CHKERRQ(ierr);

  // Use non-zero initial guess (i.e. SSA velocities from the last
  // solve() call).
  ierr = KSPSetInitialGuessNonzero(m_KSP, PETSC_TRUE); CHKERRQ(ierr);

  ierr = m_b.create(grid, "right_hand_side", WITHOUT_GHOSTS); CHKERRQ(ierr);

  ierr = m_velocity_old.create(grid, "velocity_old", WITH_GHOSTS); CHKERRQ(ierr);
  ierr = m_velocity_old.set_attrs("internal",
                                  "old SSA velocity field; used for re-trying with a different epsilon",
                                  "m s-1", ""); CHKERRQ(ierr);
  
  const double power = 1.0 / flow_law->exponent();
  char unitstr[TEMPORARY_STRING_LENGTH];
  snprintf(unitstr, sizeof(unitstr), "Pa s%f", power);
  ierr = hardness.create(grid, "hardness", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = hardness.set_attrs("diagnostic",
                            "vertically-averaged ice hardness",
                            unitstr, ""); CHKERRQ(ierr);

  ierr = nuH.create(grid, "nuH", WITH_GHOSTS); CHKERRQ(ierr);
  ierr = nuH.set_attrs("internal",
                       "ice thickness times effective viscosity",
                       "Pa s m", ""); CHKERRQ(ierr);

  ierr = nuH_old.create(grid, "nuH_old", WITH_GHOSTS); CHKERRQ(ierr);
  ierr = nuH_old.set_attrs("internal",
                           "ice thickness times effective viscosity (before an update)",
                           "Pa s m", ""); CHKERRQ(ierr);

  ierr = m_work.create(grid, "m_work", WITH_GHOSTS,
                       2, /* stencil width */
                       6  /* dof */); CHKERRQ(ierr);
  ierr = m_work.set_attrs("internal",
                          "temporary storage used to compute nuH",
                          "", ""); CHKERRQ(ierr);

  m_scaling = 1.0e9;  // comparable to typical beta for an ice stream;

  // The nuH viewer:
  view_nuh = false;
  nuh_viewer_size = 300;
  nuh_viewer = NULL;

  dump_system_matlab = false;

  return 0;
}

//! \brief De-allocate SSAFD internal objects.
PetscErrorCode SSAFD::deallocate_fd() {
  PetscErrorCode ierr;

  if (m_KSP != NULL) {
    ierr = KSPDestroy(&m_KSP); CHKERRQ(ierr);
  }

  if (m_A != NULL) {
    ierr = MatDestroy(&m_A); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode SSAFD::init(Vars &vars) {
  PetscErrorCode ierr;
  ierr = SSA::init(vars); CHKERRQ(ierr);

  // The FD solver does not support direct specification of a driving stress;
  // a surface elevation must be explicitly given.
  if (surface == NULL) {
    throw RuntimeError("The finite difference SSA solver requires a surface elevation.\n"
                       "An explicit driving stress was specified instead and cannot be used.");
  }

  ierr = verbPrintf(2,grid.com,
                    "  [using the KSP-based finite difference implementation]\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "SSAFD options", ""); CHKERRQ(ierr);
  {
    bool flag;
    ierr = OptionsInt("-ssa_nuh_viewer_size", "nuH viewer size",
                          nuh_viewer_size, flag); CHKERRQ(ierr);
    ierr = OptionsIsSet("-ssa_view_nuh", "Enable the SSAFD nuH runtime viewer",
                            view_nuh); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (config.get_flag("calving_front_stress_boundary_condition")) {
    ierr = verbPrintf(2,grid.com,
      "  using PISM-PIK calving-front stress boundary condition ...\n"); CHKERRQ(ierr);
  }

  // option to save linear system in Matlab-readable ASCII format at end of each
  // numerical solution of SSA equations; can be given with or without filename prefix
  // (i.e. "-ssa_matlab " or "-ssa_matlab foo" are both legal; in former case get
  // "pism_SSA_[year].m" if "pism_SSA" is default prefix, and in latter case get "foo_[year].m")
  std::string tempPrefix;
  ierr = OptionsIsSet("-ssafd_matlab", "Save linear system in Matlab-readable ASCII format",
                          dump_system_matlab); CHKERRQ(ierr);

  m_default_pc_failure_count     = 0;
  m_default_pc_failure_max_count = 5;

  if (config.get_flag("do_fracture_density")) {
    fracture_density = vars.get_2d_scalar("fracture_density");
  }

  return 0;
}

PetscErrorCode SSAFD::update(bool fast, IceModelVec2S& melange_back_pressure) {
  PetscErrorCode ierr;

  m_melange_back_pressure = &melange_back_pressure;

  ierr = SSA::update(fast, melange_back_pressure); CHKERRQ(ierr);

  return 0;
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
come from known velocity values.  The fields vel_bc and bc_locations are used for
this.
 */
PetscErrorCode SSAFD::assemble_rhs() {
  PetscErrorCode ierr;
  const double dx = grid.dx, dy = grid.dy;

  const double ice_free_default_velocity = 0.0;

  const double standard_gravity = config.get("standard_gravity"),
    rho_ocean = config.get("sea_water_density"),
    rho_ice = config.get("ice_density");
  const bool use_cfbc = config.get_flag("calving_front_stress_boundary_condition");

  // FIXME: bedrock_boundary is a misleading name
  bool bedrock_boundary = config.get_flag("ssa_dirichlet_bc");

  ierr = m_b.set(0.0); CHKERRQ(ierr);

  // get driving stress components
  ierr = compute_driving_stress(taud); CHKERRQ(ierr);

  IceModelVec::AccessList list;
  list.add(taud);
  list.add(m_b);

  if (m_vel_bc && bc_locations) {
    list.add(*m_vel_bc);
    list.add(*bc_locations);
  }

  if (use_cfbc) {
    list.add(*thickness);
    list.add(*bed);
    list.add(*mask);
  }

  if (use_cfbc && m_melange_back_pressure != NULL) {
    list.add(*m_melange_back_pressure);
  }

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_vel_bc != NULL &&
        bc_locations->as_int(i, j) == 1) {
      m_b(i, j).u = m_scaling * (*m_vel_bc)(i, j).u;
      m_b(i, j).v = m_scaling * (*m_vel_bc)(i, j).v;
      continue;
    }

    if (use_cfbc) {
      double H_ij = (*thickness)(i,j);
      int M_ij = mask->as_int(i,j),
        M_e = mask->as_int(i + 1,j),
        M_w = mask->as_int(i - 1,j),
        M_n = mask->as_int(i,j + 1),
        M_s = mask->as_int(i,j - 1);

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
          if (ice_free_ocean(M_e)) aPP = 0;
          if (ice_free_ocean(M_w)) aMM = 0;
          if (ice_free_ocean(M_n)) bPP = 0;
          if (ice_free_ocean(M_s)) bMM = 0;
        } else {
          if (ice_free(M_e)) aPP = 0;
          if (ice_free(M_w)) aMM = 0;
          if (ice_free(M_n)) bPP = 0;
          if (ice_free(M_s)) bMM = 0;
        }

        const double H_ij2 = H_ij*H_ij,
          b     = (*bed)(i,j),
          rho_g = rho_ice * standard_gravity;

        double ocean_pressure;
        // this is not really the ocean_pressure, but the difference
        // between ocean_pressure and isotrop.normal stresses
        // (=pressure) from within the ice

        if (ocean(M_ij)) {
          // floating shelf
          ocean_pressure = 0.5 * rho_g * (1.0 - (rho_ice / rho_ocean)) * H_ij2;
        } else {
          // grounded terminus
          if (b >= sea_level) {
            ocean_pressure = 0.5 * rho_g * H_ij2;
          } else {
            ocean_pressure = 0.5 * rho_g * (H_ij2 - (rho_ocean / rho_ice) * pow(sea_level - b, 2.0));
          }
        }

        if (m_melange_back_pressure != NULL) {
          double lambda = (*m_melange_back_pressure)(i, j);

          // adjust the "pressure imbalance term" using the provided
          // "melange back pressure fraction".
          ocean_pressure *= (1.0 - lambda);
        }

        // Note that if the current cell is "marginal" but not a CFBC
        // location, the following two lines are equaivalent to the "usual
        // case" below.
        m_b(i, j).u = taud(i,j).u + (aMM - aPP) * ocean_pressure / dx;
        m_b(i, j).v = taud(i,j).v + (bMM - bPP) * ocean_pressure / dy;

        continue;
      } // end of "if (is_marginal(i, j))"

        // If we reached this point, then CFBC are enabled, but we are in the
        // interior of a sheet or shelf. See "usual case" below.

    }   // end of "if (use_cfbc)"

    // usual case: use already computed driving stress
    m_b(i, j).u = taud(i, j).u;
    m_b(i, j).v = taud(i, j).v;
  }

  return 0;
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
PetscErrorCode SSAFD::assemble_matrix(bool include_basal_shear, Mat A) {
  PetscErrorCode  ierr;
  int zero_pivot_flag = 0;

  const double   dx=grid.dx, dy=grid.dy;
  const double   beta_ice_free_bedrock = config.get("beta_ice_free_bedrock");
  const bool use_cfbc = config.get_flag("calving_front_stress_boundary_condition");

  // FIXME: bedrock_boundary is a misleading name
  const bool bedrock_boundary = config.get_flag("ssa_dirichlet_bc");

  // shortcut:
  IceModelVec2V &vel = m_velocity;

  ierr = MatZeroEntries(A); CHKERRQ(ierr);

  IceModelVec::AccessList list;
  list.add(nuH);
  list.add(*tauc);
  list.add(vel);
  list.add(*mask);

  if (m_vel_bc && bc_locations) {
    list.add(*bc_locations);
  }
  
  const bool sub_gl = config.get_flag("sub_groundingline");
  if (sub_gl) {
    list.add(*gl_mask);
  }

  // handles friction of the ice cell along ice-free bedrock margins when bedrock higher than ice surface (in simplified setups)
  bool nu_bedrock_set=config.get_flag("nu_bedrock_set");
  if (nu_bedrock_set) {
    list.add(*thickness);
    list.add(*bed);
    list.add(*surface);
  }
  double nu_bedrock=config.get("nu_bedrock");
  double HminFrozen=0.0;

  /* matrix assembly loop */
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // Handle the easy case: provided Dirichlet boundary conditions
    if (m_vel_bc && bc_locations && bc_locations->as_int(i,j) == 1) {
      // set diagonal entry to one (scaled); RHS entry will be known velocity;
      ierr = set_diagonal_matrix_entry(A, i, j, m_scaling); CHKERRQ(ierr);
      continue;
    }

    /* Provide shorthand for the following staggered coefficients  nu H:
     *      c_n
     *  c_w     c_e
     *      c_s
     */
    // const
    double c_w = nuH(i-1,j,0);
    double c_e = nuH(i,j,0);
    double c_s = nuH(i,j-1,1);
    double c_n = nuH(i,j,1);

    if (nu_bedrock_set) {
      // if option is set, the viscosity at ice-bedrock boundary layer will
      // be prescribed and is a temperature-independent free (user determined) parameter

      // direct neighbors
      int  M_e = mask->as_int(i + 1,j),
        M_w = mask->as_int(i - 1,j),
        M_n = mask->as_int(i,j + 1),
        M_s = mask->as_int(i,j - 1);

      if ((*thickness)(i,j) > HminFrozen) {
        if ((*bed)(i-1,j) > (*surface)(i,j) && ice_free_land(M_w)) {
          c_w = nu_bedrock * 0.5 * ((*thickness)(i,j)+(*thickness)(i-1,j));
        }
        if ((*bed)(i+1,j) > (*surface)(i,j) && ice_free_land(M_e)) {
          c_e = nu_bedrock * 0.5 * ((*thickness)(i,j)+(*thickness)(i+1,j));
        }
        if ((*bed)(i,j+1) > (*surface)(i,j) && ice_free_land(M_n)) {
          c_n = nu_bedrock * 0.5 * ((*thickness)(i,j)+(*thickness)(i,j+1));
        }
        if ((*bed)(i,j-1) > (*surface)(i,j) && ice_free_land(M_s)) {
          c_s = nu_bedrock * 0.5 * ((*thickness)(i,j)+(*thickness)(i+1,j));
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

    int M_ij = mask->as_int(i,j);

    if (use_cfbc) {
      int
        // direct neighbors
        M_e = mask->as_int(i + 1,j),
        M_w = mask->as_int(i - 1,j),
        M_n = mask->as_int(i,j + 1),
        M_s = mask->as_int(i,j - 1),
        // "diagonal" neighbors
        M_ne = mask->as_int(i + 1,j + 1),
        M_se = mask->as_int(i + 1,j - 1),
        M_nw = mask->as_int(i - 1,j + 1),
        M_sw = mask->as_int(i - 1,j - 1);

      // Note: this sets velocities at both ice-free ocean and ice-free
      // bedrock to zero. This means that we need to set boundary conditions
      // at both ice/ice-free-ocean and ice/ice-free-bedrock interfaces below
      // to be consistent.
      if (ice_free(M_ij)) {
        ierr = set_diagonal_matrix_entry(A, i, j, m_scaling); CHKERRQ(ierr);
        continue;
      }

      if (is_marginal(i, j, bedrock_boundary)) {
        // If at least one of the following four conditions is "true", we're
        // at a CFBC location.
        if (bedrock_boundary) {

          if (ice_free_ocean(M_e)) aPP = 0;
          if (ice_free_ocean(M_w)) aMM = 0;
          if (ice_free_ocean(M_n)) bPP = 0;
          if (ice_free_ocean(M_s)) bMM = 0;

          // decide whether to use centered or one-sided differences
          if (ice_free_ocean(M_n) || ice_free_ocean(M_ne)) aPn = 0;
          if (ice_free_ocean(M_e) || ice_free_ocean(M_ne)) bPe = 0;
          if (ice_free_ocean(M_e) || ice_free_ocean(M_se)) bMe = 0;
          if (ice_free_ocean(M_s) || ice_free_ocean(M_se)) aPs = 0;
          if (ice_free_ocean(M_s) || ice_free_ocean(M_sw)) aMs = 0;
          if (ice_free_ocean(M_w) || ice_free_ocean(M_sw)) bMw = 0;
          if (ice_free_ocean(M_w) || ice_free_ocean(M_nw)) bPw = 0;
          if (ice_free_ocean(M_n) || ice_free_ocean(M_nw)) aMn = 0;
        } else {
          if (ice_free(M_e)) aPP = 0;
          if (ice_free(M_w)) aMM = 0;
          if (ice_free(M_n)) bPP = 0;
          if (ice_free(M_s)) bMM = 0;

          // decide whether to use centered or one-sided differences
          if (ice_free(M_n) || ice_free(M_ne)) aPn = 0;
          if (ice_free(M_e) || ice_free(M_ne)) bPe = 0;
          if (ice_free(M_e) || ice_free(M_se)) bMe = 0;
          if (ice_free(M_s) || ice_free(M_se)) aPs = 0;
          if (ice_free(M_s) || ice_free(M_sw)) aMs = 0;
          if (ice_free(M_w) || ice_free(M_sw)) bMw = 0;
          if (ice_free(M_w) || ice_free(M_nw)) bPw = 0;
          if (ice_free(M_n) || ice_free(M_nw)) aMn = 0;                           }
      }
    } // end of "if (use_cfbc)"

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
        beta = basal_sliding_law->drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);
      } else if (ice_free_land(M_ij)) {
        // apply drag even in this case, to help with margins; note ice free
        // areas already have a strength extension
        beta = beta_ice_free_bedrock;
      }
      if (sub_gl) {
        // reduce the basal drag at grid cells that are partially grounded:
        if (icy(M_ij)) {
          beta = (*gl_mask)(i,j) * basal_sliding_law->drag((*tauc)(i,j), vel(i,j).u, vel(i,j).v);
        }
      }
    }

    // add beta to diagonal entries
    eq1[4]  += beta;
    eq2[13] += beta;

    // check diagonal entries:
    const double eps = 1e-16;
    if (fabs(eq1[4]) < eps) {
      fprintf(stderr, "PISM ERROR: first  (X) equation in the SSAFD system: zero diagonal entry at a regular (not Dirichlet B.C.) location: i = %d, j = %d\n", i, j);
      zero_pivot_flag = 1;
    }
    if (fabs(eq2[13]) < eps) {
      fprintf(stderr, "PISM ERROR: second (Y) equation in the SSAFD system: zero diagonal entry at a regular (not Dirichlet B.C.) location: i = %d, j = %d\n", i, j);
      zero_pivot_flag = 1;
    }

    // build equations: NOTE TRANSPOSE
    row.j = i; row.i = j;
    for (int m = 0; m < sten; m++) {
      col[m].j = I[m]; col[m].i = J[m]; col[m].c = C[m];
    }

    // set coefficients of the first equation:
    row.c = 0;
    ierr = MatSetValuesStencil(A, 1, &row, sten, col, eq1, INSERT_VALUES); CHKERRQ(ierr);

    // set coefficients of the second equation:
    row.c = 1;
    ierr = MatSetValuesStencil(A, 1, &row, sten, col, eq2, INSERT_VALUES); CHKERRQ(ierr);
  } // loop over points

  ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
#if (PISM_DEBUG==1)
  ierr = MatSetOption(A,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE); CHKERRQ(ierr);
#endif

  int zero_pivot_flag_global = 0;
  MPI_Allreduce(&zero_pivot_flag, &zero_pivot_flag_global, 1, MPI_INT, MPI_MAX, grid.com);
  if (zero_pivot_flag_global != 0) {
    fprintf(stderr, "PISM ERROR: zero pivot detected.\n");
    return 1;
  }

  return 0;
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
PetscErrorCode SSAFD::solve() {
  PetscErrorCode ierr;

  // Store away old SSA velocity (it might be needed in case a solver
  // fails).
  ierr = m_velocity.copy_to(m_velocity_old); CHKERRQ(ierr);

  {
    // These computations do not depend on the solution, so they need
    // to be done once.
    ierr = assemble_rhs(); CHKERRQ(ierr);
    ierr = compute_hardav_staggered(); CHKERRQ(ierr);
  }

  // Try with default settings:
  {
    if (m_default_pc_failure_count < m_default_pc_failure_max_count) {
      ierr = pc_setup_bjacobi(); CHKERRQ(ierr);
    }

    ierr = picard_iteration(static_cast<int>(config.get("max_iterations_ssafd")),
                            config.get("ssafd_relative_convergence"),
                            config.get("epsilon_ssa"),
                            1.0);
  }

  if (ierr == 1) {
    // if KSP diverged (i.e. linear solver failure) then try different linear solver
    ierr = strategy_2_asm();
  }

  if (ierr == 2) {
    // if effective viscosity (i.e. outer) iteration failure then try underrelaxing the iteration
    const double underrelax = config.get("ssafd_nuH_iter_failure_underrelaxation");
    ierr = verbPrintf(1, grid.com,
                      "  re-trying with effective viscosity under-relaxation (parameter = %.2f) ...\n",
                      underrelax); CHKERRQ(ierr);
    ierr = picard_iteration(static_cast<int>(config.get("max_iterations_ssafd")),
                            config.get("ssafd_relative_convergence"),
                            config.get("epsilon_ssa"),
                            underrelax);
  }

  if (ierr != 0) {
    // try again the old strategy ... it never worked that well ...
    m_default_pc_failure_count += 1;
    ierr = strategy_1_regularization();
  }

  if (ierr != 0) {
    PetscPrintf(grid.com,
                "PISM ERROR: all SSAFD strategies failed.\n");
    write_system_petsc("allstrategiesfailed");
    return ierr;
  }

  if (dump_system_matlab) {
    ierr = write_system_matlab(""); CHKERRQ(ierr);
  }
  
  if (config.get_flag("brutal_sliding")) {
    const double brutal_sliding_scaleFactor = config.get("brutal_sliding_scale");
    ierr = m_velocity.scale(brutal_sliding_scaleFactor); CHKERRQ(ierr);

    ierr = m_velocity.update_ghosts(); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Manages the Picard iteration loop.
/*!
 *
 * Returns 0 on success, 1 if the KSP solver failed to converge, and 2 if the
 * Picard iteration failed to converge.
 *
 */
PetscErrorCode SSAFD::picard_iteration(unsigned int max_iterations,
                                       double ssa_relative_tolerance,
                                       double nuH_regularization,
                                       double nuH_iter_failure_underrelax) {
  PetscErrorCode ierr;
  double   nuH_norm, nuH_norm_change;
  PetscInt    ksp_iterations, ksp_iterations_total = 0, outer_iterations;
  KSPConvergedReason  reason;

  char tempstr[100] = "";
  bool verbose = getVerbosityLevel() >= 2,
    very_verbose = getVerbosityLevel() > 2;

  // set the initial guess:
  ierr = m_velocity.copy_to(m_velocity_global); CHKERRQ(ierr);

  stdout_ssa.clear();

  bool use_cfbc = config.get_flag("calving_front_stress_boundary_condition");

  if (use_cfbc == true) {
    ierr = compute_nuH_staggered_cfbc(nuH, nuH_regularization); CHKERRQ(ierr);
  } else {
    ierr = compute_nuH_staggered(nuH, nuH_regularization); CHKERRQ(ierr);
  }
  ierr = update_nuH_viewers(); CHKERRQ(ierr);

  // outer loop
  for (unsigned int k = 0; k < max_iterations; ++k) {

    if (very_verbose) {
      snprintf(tempstr, 100, "  %2d:", k);
      stdout_ssa += tempstr;
    }

    // in preparation of measuring change of effective viscosity:
    ierr = nuH.copy_to(nuH_old); CHKERRQ(ierr);

    // assemble (or re-assemble) matrix, which depends on updated viscosity
    ierr = assemble_matrix(true, m_A); CHKERRQ(ierr);

    if (very_verbose)
      stdout_ssa += "A:";

    // Call PETSc to solve linear system by iterative method; "inner iteration":
#if PETSC_VERSION_LT(3,5,0)
    ierr = KSPSetOperators(m_KSP, m_A, m_A, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
#else
    ierr = KSPSetOperators(m_KSP, m_A, m_A); CHKERRQ(ierr);
#endif
    ierr = KSPSolve(m_KSP, m_b.get_vec(), m_velocity_global.get_vec()); CHKERRQ(ierr); // SOLVE

    // Check if diverged; report to standard out about iteration
    ierr = KSPGetConvergedReason(m_KSP, &reason); CHKERRQ(ierr);

    if (reason < 0) {
      // KSP diverged
      ierr = verbPrintf(1, grid.com,
                        "PISM WARNING:  KSPSolve() reports 'diverged'; reason = %d = '%s'\n",
                        reason, KSPConvergedReasons[reason]); CHKERRQ(ierr);

      ierr = write_system_petsc("kspdivergederror"); CHKERRQ(ierr);

      // Tell the caller that we failed. (The caller might try again,
      // though.)
      return 1;
    }

    // report on KSP success; the "inner" iteration is done
    ierr = KSPGetIterationNumber(m_KSP, &ksp_iterations); CHKERRQ(ierr);
    ksp_iterations_total += ksp_iterations;

    if (very_verbose) {
      snprintf(tempstr, 100, "S:%d,%d: ", (int)ksp_iterations, reason);
      stdout_ssa += tempstr;
    }

    // Communicate so that we have stencil width for evaluation of effective
    // viscosity on next "outer" iteration (and geometry etc. if done):
    // Note that copy_from() updates ghosts of m_velocity.
    ierr = m_velocity_global.copy_to(m_velocity); CHKERRQ(ierr);

    // update viscosity and check for viscosity convergence
    if (use_cfbc == true) {
      ierr = compute_nuH_staggered_cfbc(nuH, nuH_regularization); CHKERRQ(ierr);
    } else {
      ierr = compute_nuH_staggered(nuH, nuH_regularization); CHKERRQ(ierr);
    }
    if (nuH_iter_failure_underrelax != 1.0) {
      ierr = nuH.scale(nuH_iter_failure_underrelax); CHKERRQ(ierr);
      ierr = nuH.add(1.0 - nuH_iter_failure_underrelax, nuH_old); CHKERRQ(ierr);
    }
    ierr = compute_nuH_norm(nuH_norm, nuH_norm_change); CHKERRQ(ierr);

    ierr = update_nuH_viewers(); CHKERRQ(ierr);

    if (very_verbose) {
      snprintf(tempstr, 100, "|nu|_2, |Delta nu|_2/|nu|_2 = %10.3e %10.3e\n",
               nuH_norm, nuH_norm_change/nuH_norm);

      stdout_ssa += tempstr;

      // assume that high verbosity shows interest in immediate
      // feedback about SSA iterations
      ierr = verbPrintf(2, grid.com, stdout_ssa.c_str()); CHKERRQ(ierr);

      stdout_ssa.clear();
    }

    outer_iterations = k + 1;

    if (nuH_norm == 0 || nuH_norm_change / nuH_norm < ssa_relative_tolerance)
      goto done;
    
  } // outer loop (k)

  // If we're here, it means that we exceeded max_iterations and still
  // failed.

  ierr = verbPrintf(1, grid.com,
                    "PISM WARNING: Effective viscosity not converged after %d outer iterations\n"
                    "  with nuH_regularization=%8.2e.\n",
                    max_iterations, nuH_regularization); CHKERRQ(ierr);

  return 2;

 done:

  if (very_verbose) {
    snprintf(tempstr, 100, "... =%5d outer iterations, ~%3.1f KSP iterations each\n",
             (int)outer_iterations, ((double) ksp_iterations_total) / outer_iterations);

    stdout_ssa += tempstr;
  } else if (verbose) {
    // at default verbosity, just record last nuH_norm_change and iterations
    snprintf(tempstr, 100, "%5d outer iterations, ~%3.1f KSP iterations each\n",
             (int)outer_iterations, ((double) ksp_iterations_total) / outer_iterations);

    stdout_ssa += tempstr;
  }

  if (verbose)
    stdout_ssa = "  SSA: " + stdout_ssa;

  return 0;
}

//! Old SSAFD recovery strategy: increase the SSA regularization parameter.
PetscErrorCode SSAFD::strategy_1_regularization() {
  PetscErrorCode ierr;
  // this has no units; epsilon goes up by this ratio when previous value failed
  const double DEFAULT_EPSILON_MULTIPLIER_SSA = 4.0;
  double nuH_regularization = config.get("epsilon_ssa");
  unsigned int k = 0, max_tries = 5;

  if (nuH_regularization <= 0.0) {
      ierr = verbPrintf(1, grid.com,
                        "PISM WARNING: will not attempt over-regularization of nuH\n"
                        "  because nuH_regularization == 0.0.\n"); CHKERRQ(ierr);
    return 1;
  }
  
  while (k < max_tries) {
    ierr = m_velocity.copy_from(m_velocity_old); CHKERRQ(ierr);
    ierr = verbPrintf(1, grid.com,
                      "  re-trying with nuH_regularization multiplied by %8.2f...\n",
                      DEFAULT_EPSILON_MULTIPLIER_SSA); CHKERRQ(ierr);

    nuH_regularization *= DEFAULT_EPSILON_MULTIPLIER_SSA;


    ierr = picard_iteration(static_cast<int>(config.get("max_iterations_ssafd")),
                            config.get("ssafd_relative_convergence"),
                            nuH_regularization,
                            1.0);
    if (ierr == 0)
      break;

    k += 1;

    if (ierr != 0 && k == max_tries)
      return ierr;
  }

  return 0;
}

//! New SSAFD strategy: switch to direct solves on sub-domains.
PetscErrorCode SSAFD::strategy_2_asm() {
  PetscErrorCode ierr;

  ierr = pc_setup_asm(); CHKERRQ(ierr);

  ierr = verbPrintf(1, grid.com,
                    "  re-trying using the Additive Schwarz preconditioner...\n"); CHKERRQ(ierr);
  
  ierr = m_velocity.copy_from(m_velocity_old); CHKERRQ(ierr);

  ierr = picard_iteration(static_cast<int>(config.get("max_iterations_ssafd")),
                          config.get("ssafd_relative_convergence"),
                          config.get("epsilon_ssa"),
                          1.0);

  return ierr;
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
PetscErrorCode SSAFD::compute_nuH_norm(double &norm, double &norm_change) {
  PetscErrorCode ierr;

  std::vector<double> nuNorm, nuChange;

  const double area = grid.dx * grid.dy;
#define MY_NORM     NORM_1

  // Test for change in nu
  ierr = nuH_old.add(-1, nuH); CHKERRQ(ierr);

  ierr = nuH_old.norm_all(MY_NORM, nuChange); CHKERRQ(ierr);
  ierr =     nuH.norm_all(MY_NORM, nuNorm);   CHKERRQ(ierr);

  nuChange[0] *= area;
  nuChange[1] *= area;
  nuNorm[0] *= area;
  nuNorm[1] *= area;

  norm_change = sqrt(PetscSqr(nuChange[0]) + PetscSqr(nuChange[1]));
  norm = sqrt(PetscSqr(nuNorm[0]) + PetscSqr(nuNorm[1]));

  return 0;
}

//! \brief Computes vertically-averaged ice hardness on the staggered grid.
PetscErrorCode SSAFD::compute_hardav_staggered() {
  PetscErrorCode ierr;
  double *E_ij, *E_offset;

  std::vector<double> E(grid.Mz);

  IceModelVec::AccessList list;
  list.add(*thickness);
  list.add(*enthalpy);
  list.add(hardness);
  list.add(*mask);

  MaskQuery m(*mask);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    ierr = enthalpy->getInternalColumn(i,j,&E_ij); CHKERRQ(ierr);
    for (int o=0; o<2; o++) {
      const int oi = 1-o, oj=o;
      double H;

      if (m.icy(i,j) && m.icy(i+oi,j+oj))
        H = 0.5 * ((*thickness)(i,j) + (*thickness)(i+oi,j+oj));
      else if (m.icy(i,j))
        H = (*thickness)(i,j);
      else
        H = (*thickness)(i+oi,j+oj);

      if (H == 0) {
        hardness(i,j,o) = -1e6; // an obviously impossible value
        continue;
      }

      ierr = enthalpy->getInternalColumn(i+oi,j+oj,&E_offset); CHKERRQ(ierr);
      // build a column of enthalpy values a the current location:
      for (unsigned int k = 0; k < grid.Mz; ++k) {
        E[k] = 0.5 * (E_ij[k] + E_offset[k]);
      }

      hardness(i,j,o) = flow_law->averaged_hardness(H, grid.kBelowHeight(H),
                                                    &grid.zlevels[0], &E[0]); CHKERRQ(ierr);
    } // o
  }     // loop over points


  ierr = fracture_induced_softening(); CHKERRQ(ierr);

  return 0;
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
PetscErrorCode SSAFD::fracture_induced_softening() {
  if (config.get_flag("do_fracture_density") == false)
    return 0;

  const double
    epsilon = config.get("fracture_density_softening_lower_limit"),
    n_glen  = flow_law->exponent();

  IceModelVec::AccessList list;
  list.add(hardness);
  list.add(*fracture_density);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    for (int o=0; o<2; o++) {
      const int oi = 1-o, oj=o;

      const double
        // fracture density on the staggered grid:
        phi       = 0.5 * ((*fracture_density)(i,j) + (*fracture_density)(i+oi,j+oj)),
        // the line below implements equation (6) in the paper
        softening = pow((1.0-(1.0-epsilon)*phi), -n_glen);

      hardness(i,j,o) *= pow(softening,-1.0/n_glen);
    }
  }

  return 0;
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
PetscErrorCode SSAFD::compute_nuH_staggered(IceModelVec2Stag &result,
                                            double nuH_regularization) {
  PetscErrorCode ierr;

  IceModelVec2V &uv = m_velocity; // shortcut

  IceModelVec::AccessList list;
  list.add(result);
  list.add(uv);
  list.add(hardness);
  list.add(*thickness);

  double ssa_enhancement_factor = flow_law->enhancement_factor(),
    n_glen = flow_law->exponent(),
    nu_enhancement_scaling = 1.0 / pow(ssa_enhancement_factor, 1.0/n_glen);

  const double dx = grid.dx, dy = grid.dy;

  for (int o=0; o<2; ++o) {
    const int oi = 1 - o, oj=o;
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double H = 0.5 * ((*thickness)(i,j) + (*thickness)(i+oi,j+oj));

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

      double nu;
      flow_law->effective_viscosity(hardness(i,j,o),
                                    secondInvariant_2D(u_x, u_y, v_x, v_y),
                                    &nu, NULL);

      result(i,j,o) = nu * H;

      // include the SSA enhancement factor; in most cases ssa_enhancement_factor is 1
      result(i,j,o) *= nu_enhancement_scaling;

      // We ensure that nuH is bounded below by a positive constant.
      result(i,j,o) += nuH_regularization;

    }
  } // o


  // Some communication
  ierr = result.update_ghosts(); CHKERRQ(ierr);
  return 0;
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
PetscErrorCode SSAFD::compute_nuH_staggered_cfbc(IceModelVec2Stag &result,
                                                 double nuH_regularization) {

  PetscErrorCode ierr;
  IceModelVec2V &uv = m_velocity; // shortcut
  double ssa_enhancement_factor = flow_law->enhancement_factor(),
    n_glen = flow_law->exponent(),
    nu_enhancement_scaling = 1.0 / pow(ssa_enhancement_factor, 1.0/n_glen);

  const unsigned int U_X = 0, V_X = 1, W_I = 2, U_Y = 3, V_Y = 4, W_J = 5;

  const double dx = grid.dx, dy = grid.dy;

  MaskQuery m(*mask);

  IceModelVec::AccessList list;
  list.add(*mask);
  list.add(m_work);
  list.add(m_velocity);

  assert(m_velocity.get_stencil_width() >= 2);
  assert(mask->get_stencil_width()      >= 2);
  assert(m_work.get_stencil_width()     >= 1);

  for (PointsWithGhosts p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // x-derivative, i-offset
    {
      if (m.icy(i,j) && m.icy(i+1,j)) {
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
      if (m.icy(i,j) && m.icy(i,j+1)) {
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

  list.add(result);
  list.add(hardness);
  list.add(*thickness);
 
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double u_x, u_y, v_x, v_y, H, nu, W;
    // i-offset
    {
      if (m.icy(i,j) && m.icy(i+1,j))
        H = 0.5 * ((*thickness)(i,j) + (*thickness)(i+1,j));
      else if (m.icy(i,j))
        H = (*thickness)(i,j);
      else
        H = (*thickness)(i+1,j);

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

        flow_law->effective_viscosity(hardness(i,j,0),
                                      secondInvariant_2D(u_x, u_y, v_x, v_y),
                                      &nu, NULL);
        result(i,j,0) = nu * H;
      } else {
        result(i,j,0) = strength_extension->get_notional_strength();
      }
    }
      
    // j-offset
    {
      if (m.icy(i,j) && m.icy(i,j+1))
        H = 0.5 * ((*thickness)(i,j) + (*thickness)(i,j+1));
      else if (m.icy(i,j))
        H = (*thickness)(i,j);
      else
        H = (*thickness)(i,j+1);

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

        flow_law->effective_viscosity(hardness(i,j,1),
                                      secondInvariant_2D(u_x, u_y, v_x, v_y),
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
  ierr = result.update_ghosts(); CHKERRQ(ierr);

  return 0;
}

//! Update the nuH viewer, which shows log10(nu H).
PetscErrorCode SSAFD::update_nuH_viewers() {
  PetscErrorCode ierr;

  if (!view_nuh) return 0;

  IceModelVec2S tmp;
  ierr = tmp.create(grid, "nuH", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = tmp.set_attrs("temporary",
                       "log10 of (viscosity * thickness)",
                       "Pa s m", ""); CHKERRQ(ierr);

  IceModelVec::AccessList list;
  list.add(nuH);
  list.add(tmp);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double avg_nuH = 0.5 * (nuH(i,j,0) + nuH(i,j,1));
    if (avg_nuH > 1.0e14) {
      tmp(i,j) = log10(avg_nuH);
    } else {
      tmp(i,j) = 14.0;
    }
  }

  if (nuh_viewer == NULL) {
    ierr = grid.create_viewer(nuh_viewer_size, "nuH", nuh_viewer); CHKERRQ(ierr);
  }

  ierr = tmp.view(nuh_viewer, NULL); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode SSAFD::set_diagonal_matrix_entry(Mat A, int i, int j,
                                                double value) {
  PetscErrorCode ierr;
  MatStencil row, col;
  row.j = i; row.i = j;
  col.j = i; col.i = j;

  row.c = 0; col.c = 0;
  ierr = MatSetValuesStencil(A, 1, &row, 1, &col, &value, INSERT_VALUES); CHKERRQ(ierr);

  row.c = 1; col.c = 1;
  ierr = MatSetValuesStencil(A, 1, &row, 1, &col, &value, INSERT_VALUES); CHKERRQ(ierr);

  return 0;
}

//! \brief Checks if a cell is near or at the ice front.
/*!
 * You need to call create IceModelVec::AccessList object and add mask to it.
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

  const int M_ij = mask->as_int(i,j),
    // direct neighbors
    M_e = mask->as_int(i + 1,j),
    M_w = mask->as_int(i - 1,j),
    M_n = mask->as_int(i,j + 1),
    M_s = mask->as_int(i,j - 1),
    // "diagonal" neighbors
    M_ne = mask->as_int(i + 1,j + 1),
    M_se = mask->as_int(i + 1,j - 1),
    M_nw = mask->as_int(i - 1,j + 1),
    M_sw = mask->as_int(i - 1,j - 1);

  if (ssa_dirichlet_bc) {
    return icy(M_ij) &&
      (ice_free(M_e) || ice_free(M_w) || ice_free(M_n) || ice_free(M_s) ||
       ice_free(M_ne) || ice_free(M_se) || ice_free(M_nw) || ice_free(M_sw));}
  else {
    return icy(M_ij) &&
      (ice_free_ocean(M_e) || ice_free_ocean(M_w) ||
       ice_free_ocean(M_n) || ice_free_ocean(M_s) ||
       ice_free_ocean(M_ne) || ice_free_ocean(M_se) ||
       ice_free_ocean(M_nw) || ice_free_ocean(M_sw));
  }
}

PetscErrorCode SSAFD::write_system_petsc(const std::string &namepart) {
  PetscErrorCode ierr;

  // write a file with a fixed filename; avoid zillions of files
  std::string filename = "SSAFD_" + namepart + ".petsc";
  ierr = verbPrintf(1, grid.com,
                    "  writing linear system to PETSc binary file %s ...\n", filename.c_str());
                    CHKERRQ(ierr);
  PetscViewer viewer;
  ierr = PetscViewerBinaryOpen(grid.com, filename.c_str(), FILE_MODE_WRITE,
                               &viewer); CHKERRQ(ierr);
  ierr = MatView(m_A, viewer); CHKERRQ(ierr);
  ierr = VecView(m_b.get_vec(), viewer); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  return 0;
}

//! \brief Write the SSA system to an .m (MATLAB) file (for debugging).
PetscErrorCode SSAFD::write_system_matlab(const std::string &namepart) {
  PetscErrorCode ierr;
  PetscViewer    viewer;
  std::string    prefix = "SSAFD_" + namepart, file_name;
  char           yearappend[PETSC_MAX_PATH_LEN];

  IceModelVec2S component;
  ierr = component.create(grid, "temp_storage", WITHOUT_GHOSTS); CHKERRQ(ierr);

  bool flag;
  ierr = OptionsString("-ssafd_matlab",
                           "Save the linear system to an ASCII .m file. Sets the file prefix.",
                           prefix, flag); CHKERRQ(ierr);

  snprintf(yearappend, PETSC_MAX_PATH_LEN, "_y%.0f.m",
           grid.convert(grid.time->current(), "seconds", "years"));
  file_name = prefix + std::string(yearappend);

  ierr = verbPrintf(2, grid.com,
                    "writing Matlab-readable file for SSAFD system A xsoln = rhs to file `%s' ...\n",
                    file_name.c_str()); CHKERRQ(ierr);
  ierr = PetscViewerCreate(grid.com, &viewer); CHKERRQ(ierr);
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII); CHKERRQ(ierr);
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(viewer, file_name.c_str()); CHKERRQ(ierr);

  // get the command which started the run
  std::string cmdstr = pism_args_string();

  // save linear system; gives system A xsoln = rhs at last (nonlinear) iteration of SSA
  ierr = PetscViewerASCIIPrintf(viewer,
    "%% A PISM linear system report for the finite difference implementation\n"
    "%% of the SSA stress balance from this run:\n"
    "%%   '%s'\n"
    "%% Writes matrix A (sparse), and vectors uv and rhs, for the linear\n"
    "%% system which was solved at the last step of the nonlinear iteration:\n"
    "%%    A * uv = rhs.\n"
    "%% Also writes the year, the coordinates x,y, their gridded versions\n"
    "%% xx,yy, and the thickness (thk) and surface elevation (usurf).\n"
    "%% Also writes i-offsetvalues of vertically-integrated viscosity\n"
    "%% (nuH_0 = nu * H), and j-offset version of same thing (nuH_1 = nu * H);\n"
    "%% these are on the staggered grid.\n",
    cmdstr.c_str());  CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"\n\necho off\n");  CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) m_A,"A"); CHKERRQ(ierr);
  ierr = MatView(m_A, viewer); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"clear zzz\n\n");  CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject) m_b.get_vec(), "rhs"); CHKERRQ(ierr);
  ierr = VecView(m_b.get_vec(), viewer); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) m_velocity_global.get_vec(), "uv"); CHKERRQ(ierr);
  ierr = VecView(m_velocity_global.get_vec(), viewer); CHKERRQ(ierr);

  // save coordinates (for viewing, primarily)
  ierr = PetscViewerASCIIPrintf(viewer,"\nyear=%10.6f;\n",
                                grid.convert(grid.time->current(), "seconds", "years"));
  CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,
            "x=%12.3f + (0:%d)*%12.3f;\n"
            "y=%12.3f + (0:%d)*%12.3f;\n",
            -grid.Lx,grid.Mx-1,grid.dx,-grid.Ly,grid.My-1,grid.dy); CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"[xx,yy]=meshgrid(x,y);\n");  CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(viewer,"echo on\n");  CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);

  return 0;
}

SSAFD_nuH::SSAFD_nuH(SSAFD *m, IceGrid &g, Vars &my_vars)
  : Diag<SSAFD>(m, g, my_vars) {

  // set metadata:
  dof = 2;
  vars.resize(dof, NCSpatialVariable(g.get_unit_system()));
  vars[0].init_2d("nuH[0]", grid);
  vars[1].init_2d("nuH[1]", grid);

  set_attrs("ice thickness times effective viscosity, i-offset", "",
            "Pa s m", "kPa s m", 0);
  set_attrs("ice thickness times effective viscosity, j-offset", "",
            "Pa s m", "kPa s m", 1);
}

PetscErrorCode SSAFD_nuH::compute(IceModelVec* &output) {
  PetscErrorCode ierr;

  IceModelVec2Stag *result = new IceModelVec2Stag;
  ierr = result->create(grid, "nuH", WITH_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->metadata(1) = vars[1];
  result->write_in_glaciological_units = true;

  ierr = model->nuH.copy_to(*result); CHKERRQ(ierr);

  output = result;
  return 0;
}

void SSAFD::get_diagnostics(std::map<std::string, Diagnostic*> &dict,
                            std::map<std::string, TSDiagnostic*> &ts_dict) {
  SSA::get_diagnostics(dict, ts_dict);

  dict["nuH"] = new SSAFD_nuH(this, grid, *variables);
}

} // end of namespace pism
