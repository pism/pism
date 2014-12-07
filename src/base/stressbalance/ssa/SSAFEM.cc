// Copyright (C) 2009--2014 Jed Brown and Ed Bueler and Constantine Khroulev and David Maxwell
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

#include "SSAFEM.hh"
#include "FETools.hh"
#include "Mask.hh"
#include "basal_resistance.hh"
#include "flowlaws.hh"
#include "pism_options.hh"
#include <stdexcept>

namespace pism {

/** The Q1 finite element SSA solver.
 *
 *
 *
 */
SSAFEM::SSAFEM(IceGrid &g, EnthalpyConverter &e, const Config &c)
  : SSA(g, e, c), m_element_index(g), m_quadrature(grid, 1.0), m_quadrature_vector(grid, 1.0) {
  PetscErrorCode ierr = allocate_fem();
  if (ierr != 0) {
    throw std::runtime_error("SSAFEM allocation failed");
  }
}

#if PETSC_VERSION_LT(3,5,0)
// FIXME: we should never define PETSc objects!
typedef PetscErrorCode (*DMDASNESJacobianLocal)(DMDALocalInfo*, void*, Mat, Mat, MatStructure*, void*);
typedef PetscErrorCode (*DMDASNESFunctionLocal)(DMDALocalInfo*, void*, void*, void*);
#endif

SSA* SSAFEMFactory(IceGrid &g, EnthalpyConverter &ec, const Config &c) {
  return new SSAFEM(g, ec, c);
}

SSAFEM::~SSAFEM() {
  deallocate_fem();
}

//! \brief Allocating SSAFEM-specific objects; called by the constructor.
PetscErrorCode SSAFEM::allocate_fem() {
  PetscErrorCode ierr;

  m_dirichletScale = 1.0;
  m_ocean_rho = config.get("sea_water_density");
  m_earth_grav = config.get("standard_gravity");
  m_beta_ice_free_bedrock = config.get("beta_ice_free_bedrock");

  ierr = SNESCreate(grid.com, &m_snes); CHKERRQ(ierr);

  // Set the SNES callbacks to call into our compute_local_function and compute_local_jacobian
  // methods via SSAFEFunction and SSAFEJ
  m_callback_data.da = *m_da;
  m_callback_data.ssa = this;
#if PETSC_VERSION_LT(3,5,0)
  ierr = DMDASNESSetFunctionLocal(*m_da, INSERT_VALUES,
                                  (DMDASNESFunctionLocal)SSAFEFunction, &m_callback_data); CHKERRQ(ierr);
  ierr = DMDASNESSetJacobianLocal(*m_da,
                                  (DMDASNESJacobianLocal)SSAFEJacobian, &m_callback_data); CHKERRQ(ierr);
#else
  ierr = DMDASNESSetFunctionLocal(*m_da, INSERT_VALUES,
                                  (DMDASNESFunction)SSAFEFunction, &m_callback_data); CHKERRQ(ierr);
  ierr = DMDASNESSetJacobianLocal(*m_da,
                                  (DMDASNESJacobian)SSAFEJacobian, &m_callback_data); CHKERRQ(ierr);
#endif

  ierr = DMSetMatType(*m_da, "baij"); CHKERRQ(ierr);
  ierr = DMSetApplicationContext(*m_da, &m_callback_data); CHKERRQ(ierr);

  ierr = SNESSetDM(m_snes, *m_da); CHKERRQ(ierr);

  // Default of maximum 200 iterations; possibly overridded by commandline
  int snes_max_it = 200;
  ierr = SNESSetTolerances(m_snes, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT,
                           snes_max_it, PETSC_DEFAULT);
  // ierr = SNESSetOptionsPrefix(snes, ((PetscObject)this)->prefix); CHKERRQ(ierr);

  ierr = SNESSetFromOptions(m_snes); CHKERRQ(ierr);

  // Allocate m_coefficients, which contains coefficient data at the quadrature points of all the elements.
  // There are nElement elements, and FEQuadrature::Nq quadrature points.
  int nElements = m_element_index.element_count();
  m_coefficients = new SSACoefficients[FEQuadrature::Nq*nElements];

  return 0;
}

//! Undo the allocations of SSAFEM::allocate_fem; called by the destructor.
PetscErrorCode SSAFEM::deallocate_fem() {
  PetscErrorCode ierr;

  ierr = SNESDestroy(&m_snes); CHKERRQ(ierr);
  delete[] m_coefficients;

  return 0;
}

// Initialize the solver, called once by the client before use.
void SSAFEM::init(Vars &vars) {

  SSA::init(vars);
  verbPrintf(2, grid.com,
             "  [using the SNES-based finite element method implementation]\n");

  setFromOptions();

  // On restart, SSA::init() reads the SSA velocity from a PISM output file
  // into IceModelVec2V "velocity". We use that field as an initial guess.
  // If we are not restarting from a PISM file, "velocity" is identically zero,
  // and the call below clears SSAX.

  m_velocity.copy_to(m_velocity_global);

  // Store coefficient data at the quadrature points.
  cacheQuadPtValues();
}

//! Opportunity to modify behaviour based on command-line options.
/*! Called from SSAFEM::init */
void SSAFEM::setFromOptions() {
  PetscErrorCode ierr;

  bool flag = false;
  m_dirichletScale = 1.0e9;
  ierr = OptionsReal("-ssa_fe_dirichlet_scale",
                     "Enforce Dirichlet conditions with this additional scaling",
                     m_dirichletScale, flag);
}

//! Solve the SSA.  The FEM solver exchanges time for memory by computing
//! the value of the various coefficients at each of the quadrature points
//! only once.  When running in an ice-model context, at each time step,
//! SSA::update is called, which calls SSAFEM::solve.  Since coefficients
//! have generally changed between timesteps, we need to recompute coefficeints
//! at the quad points. On the other hand, in the context of inversion,
//! coefficients will not change between iteration and there is no need to
//! recompute the values at the quad points.  So there are two different solve
//! methods, SSAFEM::solve() and SSAFEM::solve_nocache().  The only difference
//! is that SSAFEM::solve() recomputes the cached values of the coefficients at
//! quadrature points before calling SSAFEM::solve_nocache().
void SSAFEM::solve() {

  // Set up the system to solve (store coefficient data at the quadrature points):
  cacheQuadPtValues();

  TerminationReason::Ptr reason;
  solve_nocache(reason);
  if (reason->failed()) {
    throw RuntimeError::formatted("SSAFEM solve failed to converge (SNES reason %s)",
                                  reason->description().c_str());
  } else if (getVerbosityLevel() > 2) {
    stdout_ssa += "SSAFEM converged (SNES reason " + reason->description() + ")";
  }
}

void SSAFEM::solve(TerminationReason::Ptr &reason) {

  // Set up the system to solve (store coefficient data at the quadrature points):
  cacheQuadPtValues();

  solve_nocache(reason);
}

//! Solve the SSA without first recomputing the values of coefficients at quad
//! points.  See the disccusion of SSAFEM::solve for more discussion.
void SSAFEM::solve_nocache(TerminationReason::Ptr &reason) {
  PetscErrorCode ierr;
  PetscViewer    viewer;
  std::string filename;
  bool flag = false;

  m_epsilon_ssa = config.get("epsilon_ssa");

  OptionsString("-ssa_view", "", filename, flag);
  if (flag) {
    ierr = PetscViewerASCIIOpen(grid.com, filename.c_str(), &viewer);
    PISM_PETSC_CHK(ierr, "PetscViewerASCIIOpen");

    ierr = PetscViewerASCIIPrintf(viewer, "SNES before SSASolve_FE\n");
    PISM_PETSC_CHK(ierr, "PetscViewerASCIIPrintf");

    ierr = SNESView(m_snes, viewer); PISM_PETSC_CHK(ierr, "SNESView");
    ierr = PetscViewerASCIIPrintf(viewer, "solution vector before SSASolve_FE\n");
    PISM_PETSC_CHK(ierr, "PetscViewerASCIIPrintf");

    ierr = VecView(m_velocity_global.get_vec(), viewer);
    PISM_PETSC_CHK(ierr, "VecView");
    ierr = PetscViewerDestroy(&viewer);
    PISM_PETSC_CHK(ierr, "PetscViewerDestroy");
  }

  stdout_ssa.clear();
  if (getVerbosityLevel() >= 2) {
    stdout_ssa = "  SSA: ";
  }

  // Solve:
  ierr = SNESSolve(m_snes, NULL, m_velocity_global.get_vec()); PISM_PETSC_CHK(ierr, "SNESSolve");

  // See if it worked.
  SNESConvergedReason snes_reason;
  ierr = SNESGetConvergedReason(m_snes, &snes_reason); PISM_PETSC_CHK(ierr, "SNESGetConvergedReason");
  reason.reset(new SNESTerminationReason(snes_reason));
  if (reason->failed()) {
    return;
  }

  // Extract the solution back from SSAX to velocity and communicate.
  m_velocity_global.copy_to(m_velocity);
  m_velocity.update_ghosts();

  OptionsIsSet("-ssa_view_solution", flag);
  if (flag) {
    ierr = PetscViewerASCIIOpen(grid.com, filename.c_str(), &viewer);
    PISM_PETSC_CHK(ierr, "PetscViewerASCIIOpen");
    ierr = PetscViewerASCIIPrintf(viewer, "solution vector after SSASolve\n");
    PISM_PETSC_CHK(ierr, "PetscViewerASCIIPrintf");

    ierr = VecView(m_velocity_global.get_vec(), viewer);
    PISM_PETSC_CHK(ierr, "VecView");
    ierr = PetscViewerDestroy(&viewer);
    PISM_PETSC_CHK(ierr, "PetscViewerDestroy");
  }
}

//! Initialize stored data from the coefficients in the SSA.  Called by SSAFEM::solve.
/* This method is should be called after SSAFEM::init and whenever
any geometry or temperature related coefficients have changed. The method
stores the values of the coefficients at the quadrature points of each
element so that these interpolated values do not need to be computed
during each outer iteration of the nonlinear solve.*/
void SSAFEM::cacheQuadPtValues() {
  double
    *Enth_e[4],
    *Enth_q[4];

  double ice_density = config.get("ice_density");

  for (unsigned int q=0; q<FEQuadrature::Nq; q++) {
    Enth_q[q] = new double[grid.Mz()];
  }

  GeometryCalculator gc(sea_level, config);

  IceModelVec::AccessList list;
  list.add(*enthalpy);
  bool driving_stress_explicit;
  if ((driving_stress_x != NULL) && (driving_stress_y != NULL)) {
    driving_stress_explicit = true;
    list.add(*driving_stress_x);
    list.add(*driving_stress_y);
  } else {
    // The class SSA ensures in this case that 'surface' is available
    driving_stress_explicit = false;
    list.add(*surface);
  }

  list.add(*thickness);
  list.add(*bed);
  list.add(*tauc);

  int xs = m_element_index.xs, xm = m_element_index.xm,
    ys   = m_element_index.ys, ym = m_element_index.ym;

  for (int i=xs; i<xs+xm; i++) {
    for (int j=ys; j<ys+ym; j++) {
      double hq[FEQuadrature::Nq], hxq[FEQuadrature::Nq], hyq[FEQuadrature::Nq];
      double ds_xq[FEQuadrature::Nq], ds_yq[FEQuadrature::Nq];
      if (driving_stress_explicit) {
        m_quadrature.computeTrialFunctionValues(i, j, m_dofmap, *driving_stress_x, ds_xq);
        m_quadrature.computeTrialFunctionValues(i, j, m_dofmap, *driving_stress_y, ds_yq);
      } else {
        m_quadrature.computeTrialFunctionValues(i, j, m_dofmap, *surface, hq, hxq, hyq);
      }

      double Hq[FEQuadrature::Nq], bq[FEQuadrature::Nq], taucq[FEQuadrature::Nq];
      m_quadrature.computeTrialFunctionValues(i, j, m_dofmap, *thickness, Hq);
      m_quadrature.computeTrialFunctionValues(i, j, m_dofmap, *bed, bq);
      m_quadrature.computeTrialFunctionValues(i, j, m_dofmap, *tauc, taucq);

      const int ij = m_element_index.flatten(i, j);
      SSACoefficients *coefficients = &m_coefficients[4*ij];
      for (unsigned int q = 0; q < FEQuadrature::Nq; q++) {
        coefficients[q].H  = Hq[q];
        coefficients[q].b  = bq[q];
        coefficients[q].tauc = taucq[q];
        if (driving_stress_explicit) {
          coefficients[q].driving_stress.u = ds_xq[q];
          coefficients[q].driving_stress.v = ds_yq[q];
        } else {
          coefficients[q].driving_stress.u = -ice_density*m_earth_grav*Hq[q]*hxq[q];
          coefficients[q].driving_stress.v = -ice_density*m_earth_grav*Hq[q]*hyq[q];
        }

        coefficients[q].mask = gc.mask(coefficients[q].b, coefficients[q].H);
      }

      // In the following, we obtain the averaged hardness value from enthalpy by
      // interpolating enthalpy in each column over a quadrature point and then
      // taking the average over the column.  A faster approach would be to take
      // the column average over each element nodes and then interpolate to the
      // quadrature points. Does this make a difference?

      // Obtain the values of enthalpy at each vertical level at each of the vertices
      // of the current element.
      enthalpy->getInternalColumn(i, j, &Enth_e[0]);
      enthalpy->getInternalColumn(i+1, j, &Enth_e[1]);
      enthalpy->getInternalColumn(i+1, j+1, &Enth_e[2]);
      enthalpy->getInternalColumn(i, j+1, &Enth_e[3]);

      // We now want to interpolate to the quadrature points at each of the
      // vertical levels.  It would be nice to use quadrature::computeTestFunctionValues,
      // but the way we have just obtained the values at the element vertices
      // using getInternalColumn doesn't make this straightforward.  So we compute the values
      // by hand.
      const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();
      for (unsigned int k = 0; k < grid.Mz(); k++) {
        Enth_q[0][k] = Enth_q[1][k] = Enth_q[2][k] = Enth_q[3][k] = 0;
        for (unsigned int q = 0; q < FEQuadrature::Nq; q++) {
          for (unsigned int p = 0; p < FEQuadrature::Nk; p++) {
            Enth_q[q][k] += test[q][p].val * Enth_e[p][k];
          }
        }
      }

      // Now, for each column over a quadrature point, find the averaged_hardness.
      for (unsigned int q = 0; q < FEQuadrature::Nq; q++) {
        // Evaluate column integrals in flow law at every quadrature point's column
        coefficients[q].B = flow_law->averaged_hardness(coefficients[q].H, grid.kBelowHeight(coefficients[q].H),
                                                        &(grid.z()[0]), Enth_q[q]);
      }

    } // j-loop
  } // i-loop

  for (unsigned int q = 0; q < FEQuadrature::Nq; q++) {
    delete [] Enth_q[q];
  }
}

/** @brief Compute the "(effective viscosity) x (ice thickness)"
 *  and effective viscous bed strength from the current solution, at a
 *  single quadrature point.
 *
 * @param[in] coefficients SSA coefficients at the current quadrature point
 * @param[in] u the value of the solution
 * @param[in] Du the value of the symmetric gradient of the solution
 * @param[out] nuH product of the ice viscosity and thickness @f$ \nu H @f$
 * @param[out] dnuH derivative of @f$ \nu H @f$ with respect to the
 *                  second invariant @f$ \gamma @f$. Set to NULL if
 *                  not desired.
 * @param[out] beta basal drag coefficient @f$ \beta @f$
 * @param[out] dbeta derivative of @f$ \beta @f$ with respect to the
 *                   second invariant @f$ \gamma @f$. Set to NULL if
 *                   not desired.
 *
 * @return 0 on success
 */
void SSAFEM::PointwiseNuHAndBeta(const SSACoefficients &coefficients,
                                           const Vector2 &u, const double Du[],
                                           double *nuH, double *dnuH,
                                           double *beta, double *dbeta) {

  if (coefficients.H < strength_extension->get_min_thickness()) {
    *nuH = strength_extension->get_notional_strength();
    if (dnuH) {
      *dnuH = 0;
    }
  } else {
    flow_law->effective_viscosity(coefficients.B, secondInvariantDu_2D(Du),
                                  nuH, dnuH);

    *nuH  = m_epsilon_ssa + *nuH * coefficients.H;

    if (dnuH) {
      *dnuH *= coefficients.H;
    }
  }

  if (mask::grounded_ice(coefficients.mask)) {
    basal_sliding_law->drag_with_derivative(coefficients.tauc, u.u, u.v, beta, dbeta);
  } else {
    *beta = 0;

    if (mask::ice_free_land(coefficients.mask)) {
      *beta = m_beta_ice_free_bedrock;
    }

    if (dbeta) {
      *dbeta = 0;
    }
  }
}

//! Implements the callback for computing the SNES local function.
/*!
 * Compute the residual \f[r_{ij}= G(x, \psi_{ij}) \f] where \f$G\f$
 * is the weak form of the SSA, \f$x\f$ is the current approximate
 * solution, and the \f$\psi_{ij}\f$ are test functions.
 *
 * The weak form of the SSA system is
 */
void SSAFEM::compute_local_function(DMDALocalInfo *info,
                                              const Vector2 **velocity_global,
                                              Vector2 **residual_global) {

  (void) info; // Avoid compiler warning.

  // Zero out the portion of the function we are responsible for computing.
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    residual_global[i][j].u = 0.0;
    residual_global[i][j].v = 0.0;
  }

  // Start access to Dirichlet data if present.
  DirichletData_Vector dirichlet_data;
  dirichlet_data.init(bc_locations, m_vel_bc, m_dirichletScale);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  // Storage for the current solution at quadrature points.
  Vector2 u[FEQuadrature::Nq];
  double Du[FEQuadrature::Nq][3];

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  // Iterate over the elements.
  int xs = m_element_index.xs,
    xm   = m_element_index.xm,
    ys   = m_element_index.ys,
    ym   = m_element_index.ym;

  for (int i = xs; i < xs + xm; i++) {
    for (int j = ys; j < ys + ym; j++) {
      // Storage for element-local solution and residuals.
      Vector2 velocity[FEQuadrature::Nk], residual[FEQuadrature::Nk];

      // Index into coefficient storage in m_coefficients
      const int ij = m_element_index.flatten(i, j);

      // Coefficients and weights for this quadrature point.
      const SSACoefficients *coefficients = &m_coefficients[ij*FEQuadrature::Nq];

      // Initialize the map from global to local degrees of freedom for this element.
      m_dofmap.reset(i, j, grid);

      // Obtain the value of the solution at the nodes adjacent to the element.
      m_dofmap.extractLocalDOFs(velocity_global, velocity);

      // These values now need to be adjusted if some nodes in the element have
      // Dirichlet data.
      if (dirichlet_data) {
        dirichlet_data.update(m_dofmap, velocity);
        dirichlet_data.constrain(m_dofmap);
      }

      // Zero out the element-local residual in prep for updating it.
      for (unsigned int k = 0; k < FEQuadrature::Nk; k++) {
        residual[k].u = 0;
        residual[k].v = 0;
      }

      // Compute the solution values and symmetric gradient at the quadrature points.
      m_quadrature_vector.computeTrialFunctionValues(velocity, u, Du);

      // loop over quadrature points on this element:
      for (unsigned int q = 0; q < FEQuadrature::Nq; q++) {

        // Symmetric gradient at the quadrature point.
        const double *Duq = Du[q];

        double eta = 0.0, beta = 0.0;
        PointwiseNuHAndBeta(coefficients[q], u[q], Duq,
                            &eta, NULL, &beta, NULL);

        // The next few lines compute the actual residual for the element.
        const Vector2
          tau_b = u[q] * (- beta), // basal shear stress
          tau_d = coefficients[q].driving_stress; // gravitational driving stress

        const double
          U_x          = Duq[0],
          V_y          = Duq[1],
          U_y_plus_V_x = 2.0 * Duq[2];

        // Loop over test functions.
        for (unsigned int k = 0; k < FEQuadrature::Nk; k++) {
          const FEFunctionGerm &psi = test[q][k];
          residual[k].u += JxW[q] * (eta * (psi.dx * (4.0 * U_x + 2.0 * V_y) + psi.dy * U_y_plus_V_x)
                                     - psi.val * (tau_b.u + tau_d.u));
          residual[k].v += JxW[q] * (eta * (psi.dx * U_y_plus_V_x + psi.dy * (2.0 * U_x + 4.0 * V_y))
                                     - psi.val * (tau_b.v + tau_d.v));
        } // k
      } // q

      m_dofmap.addLocalResidualBlock(residual, residual_global);
    } // j-loop
  } // i-loop

  // Until now we have not touched rows in the residual corresponding to Dirichlet data.
  // We fix this now.
  if (dirichlet_data) {
    dirichlet_data.fix_residual(velocity_global, residual_global);
  }

  dirichlet_data.finish();

  monitor_function(velocity_global, residual_global);
}

void SSAFEM::monitor_function(const Vector2 **velocity_global,
                              Vector2 **residual_global) {
  PetscErrorCode ierr;
  bool monitorFunction = false;
  OptionsIsSet("-ssa_monitor_function", monitorFunction);
  if (monitorFunction == false) {
    return;
  }

  ierr = PetscPrintf(grid.com,
                     "SSA Solution and Function values (pointwise residuals)\n");
  PISM_PETSC_CHK(ierr, "PetscPrintf");

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    ierr = PetscSynchronizedPrintf(grid.com,
                                   "[%2d, %2d] u=(%12.10e, %12.10e)  f=(%12.4e, %12.4e)\n",
                                   i, j,
                                   velocity_global[i][j].u, velocity_global[i][j].v,
                                   residual_global[i][j].u, residual_global[i][j].v);
    
    PISM_PETSC_CHK(ierr, "PetscSynchronizedPrintf");
  }

#if PETSC_VERSION_LT(3,5,0)
  ierr = PetscSynchronizedFlush(grid.com);
  PISM_PETSC_CHK(ierr, "PetscSynchronizedFlush");
#else
  ierr = PetscSynchronizedFlush(grid.com, NULL);
  PISM_PETSC_CHK(ierr, "PetscSynchronizedFlush");
#endif
}


//! Implements the callback for computing the SNES local Jacobian.
/*!
  Compute the Jacobian

  @f[ J_{ij}{kl} \frac{d r_{ij}}{d x_{kl}}= G(x, \psi_{ij}) @f]

  where \f$G\f$ is the weak form of the SSA, \f$x\f$ is the current
approximate solution, and the \f$\psi_{ij}\f$ are test functions.

*/
void SSAFEM::compute_local_jacobian(DMDALocalInfo *info,
                                              const Vector2 **velocity_global, Mat Jac) {
  PetscErrorCode ierr;

  // Avoid compiler warning.
  (void) info;

  // Zero out the Jacobian in preparation for updating it.
  ierr = MatZeroEntries(Jac); PISM_PETSC_CHK(ierr, "MatZeroEntries");

  // Start access to Dirichlet data if present.
  DirichletData_Vector dirichlet_data;
  dirichlet_data.init(bc_locations, m_vel_bc, m_dirichletScale);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  // Storage for the current solution at quadrature points.
  Vector2 u[FEQuadrature::Nq];
  double Du[FEQuadrature::Nq][3];

  // Values of the finite element test functions at the quadrature points.
  // This is an Nq by Nk array of function germs (Nq=#of quad pts, Nk=#of test functions).
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  // Loop through all the elements.
  int
    xs = m_element_index.xs,
    xm = m_element_index.xm,
    ys = m_element_index.ys,
    ym = m_element_index.ym;

  for (int i = xs; i < xs + xm; i++) {
    for (int j = ys; j < ys + ym; j++) {
      // Values of the solution at the nodes of the current element.
      Vector2 velocity[FEQuadrature::Nk];

      // Element-local Jacobian matrix (there are FEQuadrature::Nk vector valued degrees
      // of freedom per element, for a total of (2*FEQuadrature::Nk)*(2*FEQuadrature::Nk) = 16
      // entries in the local Jacobian.
      double K[2*FEQuadrature::Nk][2*FEQuadrature::Nk];

      // Index into the coefficient storage array.
      const int ij = m_element_index.flatten(i, j);

      // Coefficients at quadrature points in the current element:
      const SSACoefficients *coefficients = &m_coefficients[ij*FEQuadrature::Nq];

      // Initialize the map from global to local degrees of freedom for this element.
      m_dofmap.reset(i, j, grid);

      // Obtain the value of the solution at the adjacent nodes to the element.
      m_dofmap.extractLocalDOFs(velocity_global, velocity);

      // These values now need to be adjusted if some nodes in the element have
      // Dirichlet data.
      if (dirichlet_data) {
        dirichlet_data.update(m_dofmap, velocity);
        dirichlet_data.constrain(m_dofmap);
      }

      // Compute the values of the solution at the quadrature points.
      m_quadrature_vector.computeTrialFunctionValues(velocity, u, Du);

      // Build the element-local Jacobian.
      ierr = PetscMemzero(K, sizeof(K));
      PISM_PETSC_CHK(ierr, "PetscMemzero");
      for (unsigned int q = 0; q < FEQuadrature::Nq; q++) {
        const double
          jw           = JxW[q],
          U            = u[q].u,
          V            = u[q].v,
          U_x          = Du[q][0],
          V_y          = Du[q][1],
          U_y_plus_V_x = 2.0 * Du[q][2]; // u_y + v_x is twice the symmetric gradient

        double eta = 0.0, deta = 0.0, beta = 0.0, dbeta = 0.0;
        PointwiseNuHAndBeta(coefficients[q], u[q], Du[q],
                            &eta, &deta, &beta, &dbeta);

        for (unsigned int l = 0; l < FEQuadrature::Nk; l++) { // Trial functions

          // Current trial function and its derivatives:
          const double
            phi   = test[q][l].val,
            phi_x = test[q][l].dx,
            phi_y = test[q][l].dy;

          // Derivatives of \gamma with respect to u_l and v_l:
          const double
            gamma_u = (2.0 * U_x + V_y) * phi_x + Du[q][2] * phi_y,
            gamma_v = Du[q][2] * phi_x + (U_x + 2.0 * V_y) * phi_y;

          // Derivatives if \eta (\nu*H) with respect to u_l and v_l:
          const double
            eta_u = deta * gamma_u,
            eta_v = deta * gamma_v;

          // Derivatives of the basal shear stress term (\tau_b):
          const double
            taub_xu = -dbeta * U * U * phi - beta * phi, // x-component, derivative with respect to u_l
            taub_xv = -dbeta * U * V * phi,              // x-component, derivative with respect to u_l
            taub_yu = -dbeta * V * U * phi,              // y-component, derivative with respect to v_l
            taub_yv = -dbeta * V * V * phi - beta * phi; // y-component, derivative with respect to v_l

          for (unsigned int k = 0; k < FEQuadrature::Nk; k++) {   // Test functions

            // Current test function and its derivatives:
            const double
              psi   = test[q][k].val,
              psi_x = test[q][k].dx,
              psi_y = test[q][k].dy;

            if (eta == 0) {
              PetscPrintf(PETSC_COMM_SELF, "eta=0 i %d j %d q %d k %d\n", i, j, q, k);
            }

            // u-u coupling
            K[k*2 + 0][l*2 + 0] += jw * (eta_u * (psi_x * (4 * U_x + 2 * V_y) + psi_y * U_y_plus_V_x)
                                         + eta * (4 * psi_x * phi_x + psi_y * phi_y) - psi * taub_xu);
            // u-v coupling
            K[k*2 + 0][l*2 + 1] += jw * (eta_v * (psi_x * (4 * U_x + 2 * V_y) + psi_y * U_y_plus_V_x)
                                         + eta * (2 * psi_x * phi_y + psi_y * phi_x) - psi * taub_xv);
            // v-u coupling
            K[k*2 + 1][l*2 + 0] += jw * (eta_u * (psi_x * U_y_plus_V_x + psi_y * (2 * U_x + 4 * V_y))
                                         + eta * (psi_x * phi_y + 2 * psi_y * phi_x) - psi * taub_yu);
            // v-v coupling
            K[k*2 + 1][l*2 + 1] += jw * (eta_v * (psi_x * U_y_plus_V_x + psi_y * (2 * U_x + 4 * V_y))
                                         + eta * (psi_x * phi_x + 4 * psi_y * phi_y) - psi * taub_yv);

          } // l
        } // k
      } // q
      m_dofmap.addLocalJacobianBlock(&K[0][0], Jac);
    } // j
  } // i

  // Until now, the rows and columns correspoinding to Dirichlet data
  // have not been set. We now put an identity block in for these
  // unknowns. Note that because we have takes steps to not touching
  // these columns previously, the symmetry of the Jacobian matrix is
  // preserved.
  if (dirichlet_data) {
    dirichlet_data.fix_jacobian(Jac);
  }

  dirichlet_data.finish();

  ierr = MatAssemblyBegin(Jac, MAT_FINAL_ASSEMBLY); PISM_PETSC_CHK(ierr, "MatAssemblyBegin");
  ierr = MatAssemblyEnd(Jac, MAT_FINAL_ASSEMBLY); PISM_PETSC_CHK(ierr, "MatAssemblyEnd");

  ierr = MatSetOption(Jac, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE); PISM_PETSC_CHK(ierr, "MatSetOption");
  ierr = MatSetOption(Jac, MAT_SYMMETRIC, PETSC_TRUE); PISM_PETSC_CHK(ierr, "MatSetOption");

  monitor_jacobian(Jac);
}

void SSAFEM::monitor_jacobian(Mat Jac) {
  PetscErrorCode ierr;
  bool mon_jac = false;
  OptionsIsSet("-ssa_monitor_jacobian", mon_jac);

  if (mon_jac == false) {
    return;
  }

  PetscViewer viewer;

  char file_name[PETSC_MAX_PATH_LEN];
  PetscInt iter = 0;
  ierr = SNESGetIterationNumber(m_snes, &iter); PISM_PETSC_CHK(ierr, "SNESGetIterationNumber");

  snprintf(file_name, PETSC_MAX_PATH_LEN, "PISM_SSAFEM_J%d.m", (int)iter);

  verbPrintf(2, grid.com,
             "writing Matlab-readable file for SSAFEM system A xsoln = rhs to file `%s' ...\n",
             file_name);

  ierr = PetscViewerCreate(grid.com, &viewer);
  PISM_PETSC_CHK(ierr, "PetscViewerCreate");
  ierr = PetscViewerSetType(viewer, PETSCVIEWERASCII);
  PISM_PETSC_CHK(ierr, "PetscViewerSetType");
  ierr = PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
  PISM_PETSC_CHK(ierr, "PetscViewerSetFormat");
  ierr = PetscViewerFileSetName(viewer, file_name);
  PISM_PETSC_CHK(ierr, "PetscViewerFileSetName");

  ierr = PetscObjectSetName((PetscObject) Jac, "A");
  PISM_PETSC_CHK(ierr, "PetscObjectSetName");
  ierr = MatView(Jac, viewer); PISM_PETSC_CHK(ierr, "MatView");
}

//!
PetscErrorCode SSAFEFunction(DMDALocalInfo *info,
			     const Vector2 **velocity, Vector2 **residual,
			     SSAFEM_SNESCallbackData *fe) {
  fe->ssa->compute_local_function(info, velocity, residual);

  return 0;
}

#if PETSC_VERSION_LT(3,5,0)
PetscErrorCode SSAFEJacobian(DMDALocalInfo *info, const Vector2 **velocity,
			     Mat A, Mat J, MatStructure *str, SSAFEM_SNESCallbackData *fe) {

  (void) A;

  fe->ssa->compute_local_jacobian(info, velocity, J);

  *str = SAME_NONZERO_PATTERN;

  return 0;
}
#else
PetscErrorCode SSAFEJacobian(DMDALocalInfo *info, const Vector2 **velocity,
			     Mat A, Mat J, SSAFEM_SNESCallbackData *fe) {

  (void) A;

  fe->ssa->compute_local_jacobian(info, velocity, J);

  return 0;
}
#endif
} // end of namespace pism
