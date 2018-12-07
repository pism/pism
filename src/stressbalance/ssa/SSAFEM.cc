// Copyright (C) 2009--2018 Jed Brown and Ed Bueler and Constantine Khroulev and David Maxwell
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

#include "pism/util/IceGrid.hh"
#include "SSAFEM.hh"
#include "pism/util/FETools.hh"
#include "pism/util/Mask.hh"
#include "pism/basalstrength/basal_resistance.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/Vars.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/geometry/Geometry.hh"

#include "pism/util/node_types.hh"

namespace pism {
namespace stressbalance {

/** The Q1 finite element SSA solver.
 *
 *
 *
 */
SSAFEM::SSAFEM(IceGrid::ConstPtr g)
  : SSA(g),
    m_bc_mask(NULL),
    m_bc_values(NULL),
    m_gc(*m_config),
    m_element_index(*g),
    m_element(*g),
    m_quadrature(g->dx(), g->dy(), 1.0) {

  const double ice_density = m_config->get_double("constants.ice.density");
  m_alpha = 1 - ice_density / m_config->get_double("constants.sea_water.density");
  m_rho_g = ice_density * m_config->get_double("constants.standard_gravity");

  m_driving_stress_x = NULL;
  m_driving_stress_y = NULL;

  PetscErrorCode ierr;

  m_dirichletScale = 1.0;
  m_beta_ice_free_bedrock = m_config->get_double("basal_resistance.beta_ice_free_bedrock");

  ierr = SNESCreate(m_grid->com, m_snes.rawptr());
  PISM_CHK(ierr, "SNESCreate");

  // Set the SNES callbacks to call into our compute_local_function and compute_local_jacobian.
  m_callback_data.da = *m_da;
  m_callback_data.ssa = this;

  ierr = DMDASNESSetFunctionLocal(*m_da, INSERT_VALUES,
                                  (DMDASNESFunction)function_callback,
                                  &m_callback_data);
  PISM_CHK(ierr, "DMDASNESSetFunctionLocal");

  ierr = DMDASNESSetJacobianLocal(*m_da,
                                  (DMDASNESJacobian)jacobian_callback,
                                  &m_callback_data);
  PISM_CHK(ierr, "DMDASNESSetJacobianLocal");

  ierr = DMSetMatType(*m_da, "baij");
  PISM_CHK(ierr, "DMSetMatType");

  ierr = DMSetApplicationContext(*m_da, &m_callback_data);
  PISM_CHK(ierr, "DMSetApplicationContext");

  ierr = SNESSetDM(m_snes, *m_da);
  PISM_CHK(ierr, "SNESSetDM");

  // Default of maximum 200 iterations; possibly overridden by command line options
  int snes_max_it = 200;
  ierr = SNESSetTolerances(m_snes, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT,
                           snes_max_it, PETSC_DEFAULT);
  PISM_CHK(ierr, "SNESSetTolerances");

  ierr = SNESSetFromOptions(m_snes);
  PISM_CHK(ierr, "SNESSetFromOptions");

  // Allocate m_coefficients, which contains coefficient data at the nodes of all the elements.
  m_coefficients.create(m_grid, "ssa_coefficients", WITH_GHOSTS, 1);

  m_node_type.create(m_grid, "node_type", WITH_GHOSTS, 1);
  m_node_type.set_attrs("internal", // intent
                        "node types: interior, boundary, exterior", // long name
                        "", ""); // no units or standard name

  // ElementMap::nodal_values() expects a ghosted IceModelVec2S. Ghosts if this field are never
  // assigned to and not communocated, though.
  m_boundary_integral.create(m_grid, "boundary_integral", WITH_GHOSTS, 1);
  m_boundary_integral.set_attrs("internal", // intent
                                "residual contribution from lateral boundaries", // long name
                                "", ""); // no units or standard name
}

SSA* SSAFEMFactory(IceGrid::ConstPtr g) {
  return new SSAFEM(g);
}

SSAFEM::~SSAFEM() {
  // empty
}

// Initialize the solver, called once by the client before use.
void SSAFEM::init_impl() {

  SSA::init_impl();

  // Use explicit driving stress if provided.
  if (m_grid->variables().is_available("ssa_driving_stress_x") and
      m_grid->variables().is_available("ssa_driving_stress_y")) {
    m_driving_stress_x = m_grid->variables().get_2d_scalar("ssa_driving_stress_x");
    m_driving_stress_y = m_grid->variables().get_2d_scalar("ssa_driving_stress_y");
  }

  m_log->message(2,
                 "  [using the SNES-based finite element method implementation]\n");

  // process command-line options
  {
    m_dirichletScale = 1.0e9;
    m_dirichletScale = options::Real("-ssa_fe_dirichlet_scale",
                                     "Enforce Dirichlet conditions with this additional scaling",
                                     m_dirichletScale);

  }

  // On restart, SSA::init() reads the SSA velocity from a PISM output file
  // into IceModelVec2V "velocity". We use that field as an initial guess.
  // If we are not restarting from a PISM file, "velocity" is identically zero,
  // and the call below clears m_velocity_global.

  m_velocity_global.copy_from(m_velocity);
}

/**  Solve the SSA system of equations.

 The FEM solver computes values of the various coefficients (most notably: the vertically-averaged
 ice hardness) at each of the element nodes *only once*.

 When running in an ice-model context, at each time step, SSA::update() is called, which calls
 SSAFEM::solve(). Since coefficients have generally changed between time steps, we need to recompute
 coefficients. On the other hand, in the context of inversion, coefficients will not change between
 iteration and there is no need to recompute their values.

 So there are two different solve methods, SSAFEM::solve() and SSAFEM::solve_nocache(). The only
 difference is that SSAFEM::solve() recomputes the cached values of the coefficients before calling
 SSAFEM::solve_nocache().
 */
void SSAFEM::solve(const Inputs &inputs) {

  TerminationReason::Ptr reason = solve_with_reason(inputs);
  if (reason->failed()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "SSAFEM solve failed to converge (SNES reason %s)",
                                  reason->description().c_str());
  } else if (m_log->get_threshold() > 2) {
    m_stdout_ssa += "SSAFEM converged (SNES reason " + reason->description() + ")";
  }
}

TerminationReason::Ptr SSAFEM::solve_with_reason(const Inputs &inputs) {

  // Set up the system to solve.
  cache_inputs(inputs);

  return solve_nocache();
}

//! Solve the SSA without first recomputing the values of coefficients at quad
//! points.  See the disccusion of SSAFEM::solve for more discussion.
TerminationReason::Ptr SSAFEM::solve_nocache() {
  PetscErrorCode ierr;

  m_epsilon_ssa = m_config->get_double("stress_balance.ssa.epsilon");

  options::String filename("-ssa_view", "");
  if (filename.is_set()) {
    petsc::Viewer viewer;
    ierr = PetscViewerASCIIOpen(m_grid->com, filename->c_str(), viewer.rawptr());
    PISM_CHK(ierr, "PetscViewerASCIIOpen");

    ierr = PetscViewerASCIIPrintf(viewer, "SNES before SSASolve_FE\n");
    PISM_CHK(ierr, "PetscViewerASCIIPrintf");

    ierr = SNESView(m_snes, viewer);
    PISM_CHK(ierr, "SNESView");

    ierr = PetscViewerASCIIPrintf(viewer, "solution vector before SSASolve_FE\n");
    PISM_CHK(ierr, "PetscViewerASCIIPrintf");

    ierr = VecView(m_velocity_global.vec(), viewer);
    PISM_CHK(ierr, "VecView");
  }

  m_stdout_ssa.clear();
  if (m_log->get_threshold() >= 2) {
    m_stdout_ssa = "  SSA: ";
  }

  // Solve:
  ierr = SNESSolve(m_snes, NULL, m_velocity_global.vec());
  PISM_CHK(ierr, "SNESSolve");

  // See if it worked.
  SNESConvergedReason snes_reason;
  ierr = SNESGetConvergedReason(m_snes, &snes_reason); PISM_CHK(ierr, "SNESGetConvergedReason");

  TerminationReason::Ptr reason(new SNESTerminationReason(snes_reason));
  if (not reason->failed()) {

    // Extract the solution back from SSAX to velocity and communicate.
    m_velocity.copy_from(m_velocity_global);
    m_velocity.update_ghosts();

    bool view_solution = options::Bool("-ssa_view_solution", "view solution of the SSA system");
    if (view_solution) {
      petsc::Viewer viewer;
      ierr = PetscViewerASCIIOpen(m_grid->com, filename->c_str(), viewer.rawptr());
      PISM_CHK(ierr, "PetscViewerASCIIOpen");

      ierr = PetscViewerASCIIPrintf(viewer, "solution vector after SSASolve\n");
      PISM_CHK(ierr, "PetscViewerASCIIPrintf");

      ierr = VecView(m_velocity_global.vec(), viewer);
      PISM_CHK(ierr, "VecView");
    }

  }

  return reason;
}

//! Initialize stored data from the coefficients in the SSA.  Called by SSAFEM::solve.
/**
   This method is should be called after SSAFEM::init and whenever any geometry or temperature
   related coefficients have changed. The method stores the values of the coefficients the nodes of
   each element so that these are available to the residual and Jacobian evaluation methods.

   We store the vertical average of the ice hardness to avoid re-doing this computation every
   iteration.

   In addition to coefficients at element nodes we store "node types" used to identify interior
   elements, exterior elements, and boundary faces.
*/
void SSAFEM::cache_inputs(const Inputs &inputs) {

  // Hold on to pointers to the B.C. mask and values: they are needed in SNES callbacks and
  // inputs.bc_{mask,values} are not available there.
  m_bc_mask   = inputs.bc_mask;
  m_bc_values = inputs.bc_values;

  const std::vector<double> &z = m_grid->z();

  IceModelVec::AccessList list{&m_coefficients,
      inputs.enthalpy,
      &inputs.geometry->ice_thickness,
      &inputs.geometry->bed_elevation,
      &inputs.geometry->sea_level_elevation,
      inputs.basal_yield_stress};

  bool use_explicit_driving_stress = (m_driving_stress_x != NULL) && (m_driving_stress_y != NULL);
  if (use_explicit_driving_stress) {
    list.add({m_driving_stress_x, m_driving_stress_y});
  }

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double thickness = inputs.geometry->ice_thickness(i, j);

      Vector2 tau_d;
      if (use_explicit_driving_stress) {
        tau_d.u = (*m_driving_stress_x)(i, j);
        tau_d.v = (*m_driving_stress_y)(i, j);
      } else {
	// tau_d above is set to zero by the Vector2
	// constructor, but is not used
      }

      const double *enthalpy = inputs.enthalpy->get_column(i, j);
      double hardness = rheology::averaged_hardness(*m_flow_law, thickness,
                                                    m_grid->kBelowHeight(thickness),
                                                    &z[0], enthalpy);

      Coefficients c;
      c.thickness      = thickness;
      c.bed            = inputs.geometry->bed_elevation(i, j);
      c.sea_level      = inputs.geometry->sea_level_elevation(i, j);
      c.tauc           = (*inputs.basal_yield_stress)(i, j);
      c.hardness       = hardness;
      c.driving_stress = tau_d;

      m_coefficients(i, j) = c;
    } // loop over owned grid points
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_coefficients.update_ghosts();

  const bool use_cfbc = m_config->get_boolean("stress_balance.calving_front_stress_bc");
  if (use_cfbc) {
    // Note: the call below uses ghosts of inputs.geometry->ice_thickness.
    compute_node_types(inputs.geometry->ice_thickness,
                       m_config->get_double("stress_balance.ice_free_thickness_standard"),
                       m_node_type);
  } else {
    m_node_type.set(NODE_INTERIOR);
  }

  cache_residual_cfbc(inputs);

}

//! Compute quadrature point values of various coefficients given a quadrature `Q` and nodal values.
void SSAFEM::quad_point_values(const fem::Quadrature &Q,
                               const Coefficients *x,
                               int *mask,
                               double *thickness,
                               double *tauc,
                               double *hardness) const {
  const fem::Germs *test = Q.test_function_values();
  const unsigned int n = Q.n();

  for (unsigned int q = 0; q < n; q++) {
    double
      bed       = 0.0,
      sea_level = 0.0;

    thickness[q] = 0.0;
    tauc[q]      = 0.0;
    hardness[q]  = 0.0;

    for (unsigned int k = 0; k < fem::q1::n_chi; k++) {
      const fem::Germ &psi  = test[q][k];

      thickness[q] += psi.val * x[k].thickness;
      bed          += psi.val * x[k].bed;
      sea_level    += psi.val * x[k].sea_level;
      tauc[q]      += psi.val * x[k].tauc;
      hardness[q]  += psi.val * x[k].hardness;
    }

    mask[q] = m_gc.mask(sea_level, bed, thickness[q]);
  }
}

//! Compute gravitational driving stress at quadrature points.
//! Uses explicitly-provided nodal values.
void SSAFEM::explicit_driving_stress(const fem::Quadrature &Q,
                                     const Coefficients *x,
                                     Vector2 *result) const {
  const fem::Germs *test = Q.test_function_values();
  const unsigned int n = Q.n();

  for (unsigned int q = 0; q < n; q++) {
    result[q] = 0.0;

    for (unsigned int k = 0; k < fem::q1::n_chi; k++) {
      const fem::Germ &psi  = test[q][k];
      result[q]  += psi.val * x[k].driving_stress;
    }
  }
}

/** Compute the the driving stress at quadrature points.
   The surface elevation @f$ h @f$ is defined by
   @f{equation}{
   h =
   \begin{cases}
   b + H, & \mathrm{grounded},\\
   z_{sl} + \alpha H, & \mathrm{shelf},
   \end{cases}
   @f}
   where @f$ b @f$ is the bed elevation, @f$ H @f$ is the ice thickness, and @f$ z_{sl} @f$ is the
   sea level, and @f$ \alpha @f$ is the ratio of densities of ice and sea water.

   Note that if @f$ b @f$, @f$ H @f$, and @f$ z_{sl} @f$ are bilinear (or @f$ C^{1} @f$ in
   general), @f$ h @f$ is continuous but not continuously differentiable. In
   other words, the gradient of @f$ h @f$ is not continuous, so near the
   grounding line we cannot compute the gravitational driving stress
   @f$ \tau_{d} = - \rho g H \nabla h @f$ using the @f$ Q_1 @f$ or @f$ P_1 @f$ FE basis
   expansion of @f$ h @f$.

   We differentiate the equation above instead, getting
   @f{equation}{
   \nabla h =
   \begin{cases}
   \nabla b + \nabla H, & \mathrm{grounded},\\
   \nabla z_{sl} + \alpha \nabla H, & \mathrm{shelf}.
   \end{cases}
   @f}

   and use basis expansions of @f$ \nabla b @f$ and @f$ \nabla H @f$.

   Overall, we have

   @f{align*}{
   \tau_{d} &= \rho g H \nabla h\\
   &=
   - \rho g H
   \begin{cases}
   \nabla b + \nabla H, & \mathrm{grounded},\\
   \alpha \nabla H, & \mathrm{shelf},
   \end{cases}
   @f}

   because @f$ z = z_{sl} @f$ defines the geoid surface and so *its gradient
   does not contribute to the driving stress*.
*/
void SSAFEM::driving_stress(const fem::Quadrature &Q,
                            const Coefficients *x,
                            Vector2 *result) const {
  const fem::Germs *test = Q.test_function_values();
  const unsigned int n = Q.n();

  for (unsigned int q = 0; q < n; q++) {
    double
      sea_level = 0.0,
      b         = 0.0,
      b_x       = 0.0,
      b_y       = 0.0,
      H         = 0.0,
      H_x       = 0.0,
      H_y       = 0.0;

    result[q] = 0.0;

    for (unsigned int k = 0; k < fem::q1::n_chi; k++) {
      const fem::Germ &psi  = test[q][k];

      b   += psi.val * x[k].bed;
      b_x += psi.dx * x[k].bed;
      b_y += psi.dy * x[k].bed;

      H   += psi.val * x[k].thickness;
      H_x += psi.dx * x[k].thickness;
      H_y += psi.dy * x[k].thickness;

      sea_level += psi.val * x[k].sea_level;
    }

    const int M = m_gc.mask(sea_level, b, H);
    const bool grounded = mask::grounded(M);

    const double
      pressure = m_rho_g * H,
      h_x = grounded ? b_x + H_x : m_alpha * H_x,
      h_y = grounded ? b_y + H_y : m_alpha * H_y;

    result[q].u = - pressure * h_x;
    result[q].v = - pressure * h_y;
  }
}


/** @brief Compute the "(regularized effective viscosity) x (ice thickness)" and effective viscous
 *  bed strength from the current solution, at a single quadrature point.
 *
 * @param[in] thickness ice thickness
 * @param[in] hardness ice hardness
 * @param[in] mask cell type mask
 * @param[in] tauc basal yield stress
 * @param[in] U the value of the solution
 * @param[in] U_x x-derivatives of velocity components
 * @param[in] U_y y-derivatives of velocity components
 * @param[out] nuH product of the ice viscosity and thickness @f$ \nu H @f$
 * @param[out] dnuH derivative of @f$ \nu H @f$ with respect to the
 *                  second invariant @f$ \gamma @f$. Set to NULL if
 *                  not desired.
 * @param[out] beta basal drag coefficient @f$ \beta @f$
 * @param[out] dbeta derivative of @f$ \beta @f$ with respect to the
 *                   second invariant @f$ \gamma @f$. Set to NULL if
 *                   not desired.
 */
void SSAFEM::PointwiseNuHAndBeta(double thickness,
                                 double hardness,
                                 int mask,
                                 double tauc,
                                 const Vector2 &U,
                                 const Vector2 &U_x,
                                 const Vector2 &U_y,
                                 double *nuH, double *dnuH,
                                 double *beta, double *dbeta) {

  if (thickness < strength_extension->get_min_thickness()) {
    *nuH = strength_extension->get_notional_strength();
    if (dnuH) {
      *dnuH = 0;
    }
  } else {
    m_flow_law->effective_viscosity(hardness, secondInvariant_2D(U_x, U_y),
                                    nuH, dnuH);

    *nuH  = m_epsilon_ssa + *nuH * thickness;

    if (dnuH) {
      *dnuH *= thickness;
    }
  }

  if (mask::grounded_ice(mask)) {
    m_basal_sliding_law->drag_with_derivative(tauc, U.u, U.v, beta, dbeta);
  } else {
    *beta = 0;

    if (mask::ice_free_land(mask)) {
      *beta = m_beta_ice_free_bedrock;
    }

    if (dbeta) {
      *dbeta = 0;
    }
  }
}


//! Compute and cache residual contributions from the integral over the lateral boundary.
/**

   This method computes FIXME.

 */
void SSAFEM::cache_residual_cfbc(const Inputs &inputs) {

  fem::BoundaryQuadrature2 bq(m_grid->dx(), m_grid->dy(), 1.0);

  const unsigned int Nk = fem::q1::n_chi;
  const unsigned int Nq = bq.n();
  const unsigned int n_sides = fem::q1::n_sides;
  using mask::ocean;

  const Vector2 *outward_normal = fem::q1::outward_normals();

  const bool
    use_cfbc          = m_config->get_boolean("stress_balance.calving_front_stress_bc"),
    is_dry_simulation = m_config->get_boolean("ocean.always_grounded");

  const double
    ice_density      = m_config->get_double("constants.ice.density"),
    ocean_density    = m_config->get_double("constants.sea_water.density"),
    standard_gravity = m_config->get_double("constants.standard_gravity");

  // Reset the boundary integral so that all values are overwritten.
  m_boundary_integral.set(0.0);

  if (not use_cfbc) {
    // If CFBC is not used then we're done.
    return;
  }

  IceModelVec::AccessList list{&m_node_type,
      &inputs.geometry->ice_thickness,
      &inputs.geometry->bed_elevation,
      &inputs.geometry->sea_level_elevation,
      &m_boundary_integral};

  // Iterate over the elements.
  const int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  ParallelSection loop(m_grid->com);
  try {
    for (int j = ys; j < ys + ym; j++) {
      for (int i = xs; i < xs + xm; i++) {
        // Initialize the map from global to element degrees of freedom.
        m_element.reset(i, j);

        int node_type[Nk];
        m_element.nodal_values(m_node_type, node_type);

        // an element is "interior" if all its nodes are interior or boundary
        const bool interior_element = (node_type[0] < NODE_EXTERIOR and
                                       node_type[1] < NODE_EXTERIOR and
                                       node_type[2] < NODE_EXTERIOR and
                                       node_type[3] < NODE_EXTERIOR);

        // residual contributions at element nodes
        std::vector<Vector2> I(Nk);

        if (not interior_element) {
          // not an interior element: the contribution is zero
          continue;
        }

        double H_nodal[Nk];
        m_element.nodal_values(inputs.geometry->ice_thickness, H_nodal);

        double b_nodal[Nk];
        m_element.nodal_values(inputs.geometry->bed_elevation, b_nodal);

        double sl_nodal[Nk];
        m_element.nodal_values(inputs.geometry->sea_level_elevation, sl_nodal);

        // storage for test function values psi[0] for the first "end" of a side, psi[1] for the
        // second
        double psi[2] = {0.0, 0.0};

        // loop over element sides
        for (unsigned int side = 0; side < n_sides; ++side) {

          // nodes incident to the current side
          const int
            n0 = fem::q1::incident_nodes[side][0],
            n1 = fem::q1::incident_nodes[side][1];

          if (not (node_type[n0] == NODE_BOUNDARY and node_type[n1] == NODE_BOUNDARY)) {
            // not a boundary side; skip it
            continue;
          }

          // in our case (i.e. uniform spacing in x and y directions) W is the same at all
          // quadrature points along a side.
          const double W = bq.weights(side);

          for (unsigned int q = 0; q < Nq; ++q) {

            // test functions at nodes incident to the current side, evaluated at the quadrature point
            // q
            psi[0] = bq.germ(side, q, n0).val;
            psi[1] = bq.germ(side, q, n1).val;

            // Compute ice thickness and bed elevation at a quadrature point. This uses a 1D basis
            // expansion on the side.
            const double
              H         = H_nodal[n0]  * psi[0] + H_nodal[n1]  * psi[1],
              bed       = b_nodal[n0]  * psi[0] + b_nodal[n1]  * psi[1],
              sea_level = sl_nodal[n0] * psi[0] + sl_nodal[n1] * psi[1];

            const bool floating = ocean(m_gc.mask(sea_level, bed, H));

            // ocean pressure difference at a quadrature point
            const double dP = ocean_pressure_difference(floating, is_dry_simulation,
                                                        H, bed, sea_level,
                                                        ice_density, ocean_density,
                                                        standard_gravity);

            // This integral contributes to the residual at 2 nodes (the ones incident to the
            // current side). This is is written in a way that allows *adding* (... += ...) the
            // boundary contribution in the residual computation.
            I[n0] += W * (- psi[0] * dP) * outward_normal[side];
            I[n1] += W * (- psi[1] * dP) * outward_normal[side];
          } // q-loop

        } // loop over element sides

        m_element.add_contribution(&I[0], m_boundary_integral);

      } // i-loop
    } // j-loop
  } catch (...) {
    loop.failed();
  }
  loop.check();
}


//! Implements the callback for computing the residual.
/*!
 * Compute the residual \f[r_{ij}= G(x, \psi_{ij}) \f] where \f$G\f$
 * is the weak form of the SSA, \f$x\f$ is the current approximate
 * solution, and the \f$\psi_{ij}\f$ are test functions.
 *
 * The weak form of the SSA system is
 */
void SSAFEM::compute_local_function(Vector2 const *const *const velocity_global,
                                    Vector2 **residual_global) {

  const bool use_explicit_driving_stress = (m_driving_stress_x != NULL) && (m_driving_stress_y != NULL);

  const bool use_cfbc = m_config->get_boolean("stress_balance.calving_front_stress_bc");

  const unsigned int Nk = fem::q1::n_chi;
  const unsigned int Nq_max = fem::MAX_QUADRATURE_SIZE;

  IceModelVec::AccessList list{&m_node_type, &m_coefficients, &m_boundary_integral};

  // Set the boundary contribution of the residual. This is computed at the nodes, so we don't want
  // to set it using ElementMap::add_contribution() because that would lead to
  // double-counting. Also note that without CFBC m_boundary_integral is exactly zero.
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    residual_global[j][i] = m_boundary_integral(i, j);
  }

  // Start access to Dirichlet data if present.
  fem::DirichletData_Vector dirichlet_data(m_bc_mask, m_bc_values, m_dirichletScale);

  // Storage for the current solution and its derivatives at quadrature points.
  Vector2 U[Nq_max], U_x[Nq_max], U_y[Nq_max];

  // Iterate over the elements.
  const int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  ParallelSection loop(m_grid->com);
  try {
    for (int j = ys; j < ys + ym; j++) {
      for (int i = xs; i < xs + xm; i++) {
        // Initialize the map from global to element degrees of freedom.
        m_element.reset(i, j);
        int node_type[Nk];
        m_element.nodal_values(m_node_type, node_type);

        // an element is "interior" if all its nodes are interior or boundary
        const bool interior_element = (node_type[0] < NODE_EXTERIOR and
                                       node_type[1] < NODE_EXTERIOR and
                                       node_type[2] < NODE_EXTERIOR and
                                       node_type[3] < NODE_EXTERIOR);

        if (use_cfbc and (not interior_element)) {
          // an exterior element in the CFBC case
          continue;
        }

        // Note: without CFBC all elements are "interior".

        fem::Quadrature &Q = m_quadrature;

        // Number of quadrature points.
        const unsigned int Nq = Q.n();

        // An Nq by Nk array of test function values.
        const fem::Germs *test = Q.test_function_values();

        // Jacobian times weights for quadrature.
        const double* W = Q.weights();

        // Storage for the solution and residuals at element nodes.
        Vector2 residual[Nk];

        int    mask[Nq_max];
        double thickness[Nq_max];
        double tauc[Nq_max];
        double hardness[Nq_max];
        Vector2 tau_d[Nq_max];

        {
          Coefficients coeffs[Nk];
          m_element.nodal_values(m_coefficients, coeffs);

          quad_point_values(Q, coeffs, mask, thickness, tauc, hardness);

          if (use_explicit_driving_stress) {
            explicit_driving_stress(Q, coeffs, tau_d);
          } else {
            driving_stress(Q, coeffs, tau_d);
          }
        }

        {
          // Obtain the value of the solution at the nodes adjacent to the element.
          Vector2 velocity_nodal[Nk];
          m_element.nodal_values(velocity_global, velocity_nodal);

          // These values now need to be adjusted if some nodes in the element have Dirichlet data.
          if (dirichlet_data) {
            // Set elements of velocity_nodal that correspond to Dirichlet nodes to prescribed
            // values.
            dirichlet_data.enforce(m_element, velocity_nodal);
            // mark Dirichlet nodes in m_element so that they are not touched by
            // add_contribution() below
            dirichlet_data.constrain(m_element);
          }

          // Compute the solution values and its gradient at the quadrature points.
          quadrature_point_values(Q, velocity_nodal, // input
                                  U, U_x, U_y);   // outputs
        }

        // Zero out the element-local residual in preparation for updating it.
        for (unsigned int k = 0; k < Nk; k++) {
          residual[k].u = 0;
          residual[k].v = 0;
        }

        // loop over quadrature points:
        for (unsigned int q = 0; q < Nq; q++) {

          double eta = 0.0, beta = 0.0;
          PointwiseNuHAndBeta(thickness[q], hardness[q], mask[q], tauc[q],
                              U[q], U_x[q], U_y[q], // inputs
                              &eta, NULL, &beta, NULL);              // outputs

          // The next few lines compute the actual residual for the element.
          const Vector2 tau_b = U[q] * (- beta); // basal shear stress

          const double
            jw           = W[q],
            u_x          = U_x[q].u,
            v_y          = U_y[q].v,
            u_y_plus_v_x = U_y[q].u + U_x[q].v;

          // Loop over test functions.
          for (unsigned int k = 0; k < Nk; k++) {
            const fem::Germ &psi = test[q][k];

            residual[k].u += jw * (eta * (psi.dx * (4.0 * u_x + 2.0 * v_y) + psi.dy * u_y_plus_v_x)
                                   - psi.val * (tau_b.u + tau_d[q].u));
            residual[k].v += jw * (eta * (psi.dx * u_y_plus_v_x + psi.dy * (2.0 * u_x + 4.0 * v_y))
                                   - psi.val * (tau_b.v + tau_d[q].v));
          } // k (test functions)
        }   // q (quadrature points)

        m_element.add_contribution(residual, residual_global);
      } // i-loop
    } // j-loop
  } catch (...) {
    loop.failed();
  }
  loop.check();

  // Until now we have not touched rows in the residual corresponding to Dirichlet data.
  // We fix this now.
  if (dirichlet_data) {
    dirichlet_data.fix_residual(velocity_global, residual_global);
  }

  if (use_cfbc) {
    // Prescribe homogeneous Dirichlet B.C. at ice-free nodes. This uses the fact that m_node_type
    // can be used as a "Dirichlet B.C. mask", i.e. ice-free nodes (and only ice-free nodes) are
    // marked with ones.
    fem::DirichletData_Vector dirichlet_ice_free(&m_node_type, NULL, m_dirichletScale);
    dirichlet_ice_free.fix_residual_homogeneous(residual_global);
  }

  monitor_function(velocity_global, residual_global);
}

void SSAFEM::monitor_function(Vector2 const *const *const velocity_global,
                              Vector2 const *const *const residual_global) {
  PetscErrorCode ierr;
  bool monitorFunction = options::Bool("-ssa_monitor_function", "monitor the SSA residual");
  if (not monitorFunction) {
    return;
  }

  ierr = PetscPrintf(m_grid->com,
                     "SSA Solution and Function values (pointwise residuals)\n");
  PISM_CHK(ierr, "PetscPrintf");

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      ierr = PetscSynchronizedPrintf(m_grid->com,
                                     "[%2d, %2d] u=(%12.10e, %12.10e)  f=(%12.4e, %12.4e)\n",
                                     i, j,
                                     velocity_global[j][i].u, velocity_global[j][i].v,
                                     residual_global[j][i].u, residual_global[j][i].v);
      PISM_CHK(ierr, "PetscSynchronizedPrintf");
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  ierr = PetscSynchronizedFlush(m_grid->com, NULL);
  PISM_CHK(ierr, "PetscSynchronizedFlush");
}


//! Implements the callback for computing the Jacobian.
/*!
  Compute the Jacobian

  @f[ J_{ij}{kl} \frac{d r_{ij}}{d x_{kl}}= G(x, \psi_{ij}) @f]

  where \f$G\f$ is the weak form of the SSA, \f$x\f$ is the current
  approximate solution, and the \f$\psi_{ij}\f$ are test functions.

*/
void SSAFEM::compute_local_jacobian(Vector2 const *const *const velocity_global, Mat Jac) {

  const unsigned int Nk     = fem::q1::n_chi;
  const unsigned int Nq_max = fem::MAX_QUADRATURE_SIZE;

  const bool use_cfbc = m_config->get_boolean("stress_balance.calving_front_stress_bc");

  // Zero out the Jacobian in preparation for updating it.
  PetscErrorCode ierr = MatZeroEntries(Jac);
  PISM_CHK(ierr, "MatZeroEntries");

  IceModelVec::AccessList list{&m_node_type, &m_coefficients};

  // Start access to Dirichlet data if present.
  fem::DirichletData_Vector dirichlet_data(m_bc_mask, m_bc_values, m_dirichletScale);

  // Storage for the current solution at quadrature points.
  Vector2 U[Nq_max], U_x[Nq_max], U_y[Nq_max];

  // Loop through all the elements.
  int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  ParallelSection loop(m_grid->com);
  try {
    for (int j = ys; j < ys + ym; j++) {
      for (int i = xs; i < xs + xm; i++) {
        // Initialize the map from global to element degrees of freedom.
        m_element.reset(i, j);

        int node_type[Nk];
        m_element.nodal_values(m_node_type, node_type);
        // an element is "interior" if all its nodes are interior or boundary
        const bool interior_element = (node_type[0] < NODE_EXTERIOR and
                                       node_type[1] < NODE_EXTERIOR and
                                       node_type[2] < NODE_EXTERIOR and
                                       node_type[3] < NODE_EXTERIOR);

        if (use_cfbc and (not interior_element)) {
          // an exterior element in the CFBC case
          continue;
        }

        fem::Quadrature &Q = m_quadrature;

        // Number of quadrature points.
        const unsigned int Nq = Q.n();

        // Jacobian times weights for quadrature.
        const double* W = Q.weights();

        // Values of the finite element test functions at the quadrature points.
        // This is an Nq by Nk array of function germs
        const fem::Germs *test = Q.test_function_values();

        int    mask[Nq_max];
        double thickness[Nq_max];
        double tauc[Nq_max];
        double hardness[Nq_max];

        {
          Coefficients coeffs[Nk];
          m_element.nodal_values(m_coefficients, coeffs);

          quad_point_values(Q, coeffs,
                            mask, thickness, tauc, hardness);
        }

        {
          // Values of the solution at the nodes of the current element.
          Vector2 velocity_nodal[Nk];
          // Obtain the value of the solution at the adjacent nodes to the element.
          m_element.nodal_values(velocity_global, velocity_nodal);

          // These values now need to be adjusted if some nodes in the element have
          // Dirichlet data.
          if (dirichlet_data) {
            dirichlet_data.enforce(m_element, velocity_nodal);
            dirichlet_data.constrain(m_element);
          }
          // Compute the values of the solution at the quadrature points.
          quadrature_point_values(Q, velocity_nodal, U, U_x, U_y);
        }

        // Element-local Jacobian matrix (there are Nk vector valued degrees
        // of freedom per element, for a total of (2*Nk)*(2*Nk) = 16
        // entries in the local Jacobian.
        double K[2*Nk][2*Nk];
        // Build the element-local Jacobian.
        ierr = PetscMemzero(K, sizeof(K));
        PISM_CHK(ierr, "PetscMemzero");

        for (unsigned int q = 0; q < Nq; q++) {
          const double
            jw           = W[q],
            u            = U[q].u,
            v            = U[q].v,
            u_x          = U_x[q].u,
            v_y          = U_y[q].v,
            u_y_plus_v_x = U_y[q].u + U_x[q].v;

          double eta = 0.0, deta = 0.0, beta = 0.0, dbeta = 0.0;
          PointwiseNuHAndBeta(thickness[q], hardness[q], mask[q], tauc[q],
                              U[q], U_x[q], U_y[q],
                              &eta, &deta, &beta, &dbeta);

          for (unsigned int l = 0; l < Nk; l++) { // Trial functions

            // Current trial function and its derivatives:
            const fem::Germ &phi = test[q][l];

            // Derivatives of \gamma with respect to u_l and v_l:
            const double
              gamma_u = (2.0 * u_x + v_y) * phi.dx + 0.5 * u_y_plus_v_x * phi.dy,
              gamma_v = 0.5 * u_y_plus_v_x * phi.dx + (u_x + 2.0 * v_y) * phi.dy;

            // Derivatives of \eta = \nu*H with respect to u_l and v_l:
            const double
              eta_u = deta * gamma_u,
              eta_v = deta * gamma_v;

            // Derivatives of the basal shear stress term (\tau_b):
            const double
              taub_xu = -dbeta * u * u * phi.val - beta * phi.val,  // x-component, derivative with respect to u_l
              taub_xv = -dbeta * u * v * phi.val,                   // x-component, derivative with respect to u_l
              taub_yu = -dbeta * v * u * phi.val,                   // y-component, derivative with respect to v_l
              taub_yv = -dbeta * v * v * phi.val - beta * phi.val;  // y-component, derivative with respect to v_l

            for (unsigned int k = 0; k < Nk; k++) {   // Test functions

              // Current test function and its derivatives:

              const fem::Germ &psi = test[q][k];

              if (eta == 0) {
                ierr = PetscPrintf(PETSC_COMM_SELF, "eta=0 i %d j %d q %d k %d\n", i, j, q, k);
                PISM_CHK(ierr, "PetscPrintf");
              }

              // u-u coupling
              K[k*2 + 0][l*2 + 0] += jw * (eta_u * (psi.dx * (4 * u_x + 2 * v_y) + psi.dy * u_y_plus_v_x)
                                           + eta * (4 * psi.dx * phi.dx + psi.dy * phi.dy) - psi.val * taub_xu);
              // u-v coupling
              K[k*2 + 0][l*2 + 1] += jw * (eta_v * (psi.dx * (4 * u_x + 2 * v_y) + psi.dy * u_y_plus_v_x)
                                           + eta * (2 * psi.dx * phi.dy + psi.dy * phi.dx) - psi.val * taub_xv);
              // v-u coupling
              K[k*2 + 1][l*2 + 0] += jw * (eta_u * (psi.dx * u_y_plus_v_x + psi.dy * (2 * u_x + 4 * v_y))
                                           + eta * (psi.dx * phi.dy + 2 * psi.dy * phi.dx) - psi.val * taub_yu);
              // v-v coupling
              K[k*2 + 1][l*2 + 1] += jw * (eta_v * (psi.dx * u_y_plus_v_x + psi.dy * (2 * u_x + 4 * v_y))
                                           + eta * (psi.dx * phi.dx + 4 * psi.dy * phi.dy) - psi.val * taub_yv);

            } // l
          } // k
        } // q
        m_element.add_contribution(&K[0][0], Jac);
      } // j
    } // i
  } catch (...) {
    loop.failed();
  }
  loop.check();


  // Until now, the rows and columns correspoinding to Dirichlet data
  // have not been set. We now put an identity block in for these
  // unknowns. Note that because we have takes steps to not touching
  // these columns previously, the symmetry of the Jacobian matrix is
  // preserved.
  if (dirichlet_data) {
    dirichlet_data.fix_jacobian(Jac);
  }

  if (use_cfbc) {
    // Prescribe homogeneous Dirichlet B.C. at ice-free nodes. This uses the fact that m_node_type
    // can be used as a "Dirichlet B.C. mask", i.e. ice-free nodes (and only ice-free nodes) are
    // marked with ones.
    fem::DirichletData_Vector dirichlet_ice_free(&m_node_type, NULL, m_dirichletScale);
    dirichlet_ice_free.fix_jacobian(Jac);
  }

  ierr = MatAssemblyBegin(Jac, MAT_FINAL_ASSEMBLY);
  PISM_CHK(ierr, "MatAssemblyBegin");

  ierr = MatAssemblyEnd(Jac, MAT_FINAL_ASSEMBLY);
  PISM_CHK(ierr, "MatAssemblyEnd");

  ierr = MatSetOption(Jac, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
  PISM_CHK(ierr, "MatSetOption");

  ierr = MatSetOption(Jac, MAT_SYMMETRIC, PETSC_TRUE);
  PISM_CHK(ierr, "MatSetOption");

  monitor_jacobian(Jac);
}

void SSAFEM::monitor_jacobian(Mat Jac) {
  PetscErrorCode ierr;
  bool mon_jac = options::Bool("-ssa_monitor_jacobian", "monitor the SSA Jacobian");

  if (not mon_jac) {
    return;
  }

  // iter has to be a PetscInt because it is used in the
  // SNESGetIterationNumber() call below.
  PetscInt iter = 0;
  ierr = SNESGetIterationNumber(m_snes, &iter);
  PISM_CHK(ierr, "SNESGetIterationNumber");

  char file_name[PETSC_MAX_PATH_LEN];
  snprintf(file_name, PETSC_MAX_PATH_LEN, "PISM_SSAFEM_J%d.m", (int)iter);

  m_log->message(2,
                 "writing Matlab-readable file for SSAFEM system A xsoln = rhs to file `%s' ...\n",
                 file_name);

  petsc::Viewer viewer(m_grid->com);

  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);
  PISM_CHK(ierr, "PetscViewerSetType");

  ierr = PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  PISM_CHK(ierr, "PetscViewerPushFormat");

  ierr = PetscViewerFileSetName(viewer, file_name);
  PISM_CHK(ierr, "PetscViewerFileSetName");

  ierr = PetscObjectSetName((PetscObject) Jac, "A");
  PISM_CHK(ierr, "PetscObjectSetName");

  ierr = MatView(Jac, viewer);
  PISM_CHK(ierr, "MatView");

  ierr = PetscViewerPopFormat(viewer);
  PISM_CHK(ierr, "PetscViewerPopFormat");
}

//!
PetscErrorCode SSAFEM::function_callback(DMDALocalInfo *info,
                                         Vector2 const *const *const velocity,
                                         Vector2 **residual,
                                         CallbackData *fe) {
  try {
    (void) info;
    fe->ssa->compute_local_function(velocity, residual);
  } catch (...) {
    MPI_Comm com = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)fe->da, &com); CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

PetscErrorCode SSAFEM::jacobian_callback(DMDALocalInfo *info,
                                         Vector2 const *const *const velocity,
                                         Mat A, Mat J, CallbackData *fe) {
  try {
    (void) A;
    (void) info;
    fe->ssa->compute_local_jacobian(velocity, J);
  } catch (...) {
    MPI_Comm com = MPI_COMM_SELF;
    PetscErrorCode ierr = PetscObjectGetComm((PetscObject)fe->da, &com); CHKERRQ(ierr);
    handle_fatal_errors(com);
    SETERRQ(com, 1, "A PISM callback failed");
  }
  return 0;
}

} // end of namespace stressbalance
} // end of namespace pism
