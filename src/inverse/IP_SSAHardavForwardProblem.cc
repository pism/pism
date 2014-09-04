// Copyright (C) 2013, 2014  David Maxwell
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

#include "IP_SSAHardavForwardProblem.hh"
#include <assert.h>
#include "PISMVars.hh"
#include "Mask.hh"
#include "basal_resistance.hh"

namespace pism {


IP_SSAHardavForwardProblem::IP_SSAHardavForwardProblem(IceGrid &g, EnthalpyConverter &e,
                                                       IPDesignVariableParameterization &tp,
                                                       const Config &c)
  : SSAFEM(g, e, c),
    m_grid(grid), m_zeta(NULL),
    m_fixed_design_locations(NULL),
    m_design_param(tp),
    m_element_index(m_grid),
    m_quadrature(grid, 1.0),
    m_rebuild_J_state(true) {
  PetscErrorCode ierr = this->construct();
  CHKERRCONTINUE(ierr);
  assert(ierr == 0);
}

IP_SSAHardavForwardProblem::~IP_SSAHardavForwardProblem() {
  PetscErrorCode ierr = this->destruct();
  CHKERRCONTINUE(ierr);
  assert(ierr == 0);
}

PetscErrorCode IP_SSAHardavForwardProblem::construct() {
  PetscErrorCode ierr;
  int stencilWidth = 1;

  ierr = m_dzeta_local.create(m_grid, "d_zeta_local", WITH_GHOSTS, stencilWidth); CHKERRQ(ierr);
  ierr = m_hardav.create(m_grid, "hardav", WITH_GHOSTS, stencilWidth); CHKERRQ(ierr);

  ierr = m_du_global.create(m_grid, "linearization work vector (sans ghosts)", WITHOUT_GHOSTS, stencilWidth); CHKERRQ(ierr);
  ierr = m_du_local.create(m_grid, "linearization work vector (with ghosts)", WITH_GHOSTS, stencilWidth); CHKERRQ(ierr);

  ierr = DMCreateMatrix(SSADA->get(), "baij", &m_J_state); CHKERRQ(ierr);

  ierr = KSPCreate(m_grid.com, &m_ksp); CHKERRQ(ierr);
  double ksp_rtol = 1e-12;
  ierr = KSPSetTolerances(m_ksp, ksp_rtol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(m_ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCBJACOBI); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(m_ksp); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IP_SSAHardavForwardProblem::destruct() {
  PetscErrorCode ierr;
  ierr = MatDestroy(&m_J_state); CHKERRQ(ierr);
  ierr = KSPDestroy(&m_ksp); CHKERRQ(ierr);
  return 0;
}

//! Sets the current value of of the design paramter \f$\zeta\f$.
/*! This method sets \f$\zeta\f$ but does not solve the %SSA.
It it intended for inverse methods that simultaneously compute
the pair \f$u\f$ and \f$\zeta\f$ without ever solving the %SSA
directly.  Use this method in conjuction with
\ref assemble_jacobian_state and \ref apply_jacobian_design and their friends.
The vector \f$\zeta\f$ is not copied; a reference to the IceModelVec is
kept.
*/
PetscErrorCode IP_SSAHardavForwardProblem::set_design(IceModelVec2S &new_zeta)
{
  PetscErrorCode ierr;

  m_zeta = &new_zeta;

  // Convert zeta to hardav.
  m_design_param.convertToDesignVariable(*m_zeta, m_hardav);

  // Cache hardav at the quadrature points in m_coefficients.
  double hardav_q[FEQuadrature::Nq];
  ierr = m_hardav.begin_access(); CHKERRQ(ierr);
  int xs = m_element_index.xs, xm = m_element_index.xm,
    ys = m_element_index.ys, ym = m_element_index.ym;
  for (int i = xs; i < xs + xm; i++) {
    for (int j = ys; j < ys + ym; j++) {
      m_quadrature.computeTrialFunctionValues(i, j, m_dofmap, m_hardav, hardav_q);
      const int ij = m_element_index.flatten(i, j);
      SSACoefficients *coefficients = &m_coefficients[ij*FEQuadrature::Nq];
      for (unsigned int q = 0; q < FEQuadrature::Nq; q++) {
        coefficients[q].B = hardav_q[q];
      }
    }
  }
  ierr = m_hardav.end_access(); CHKERRQ(ierr);

  // Flag the state jacobian as needing rebuilding.
  m_rebuild_J_state = true;

  return 0;
}

//! Sets the current value of the design variable \f$\zeta\f$ and solves the %SSA to find the associated \f$u_{\rm SSA}\f$.
/* Use this method for inverse methods employing the reduced gradient. Use this method
in conjuction with apply_linearization and apply_linearization_transpose.*/
PetscErrorCode IP_SSAHardavForwardProblem::linearize_at(IceModelVec2S &zeta, TerminationReason::Ptr &reason) {

  PetscErrorCode ierr;
  ierr = this->set_design(zeta); CHKERRQ(ierr);

  ierr = this->solve_nocache(reason); CHKERRQ(ierr);

  return 0;
}

//! Computes the residual function \f$\mathcal{R}(u, \zeta)\f$ as defined in the class-level documentation.
/* The value of \f$\zeta\f$ is set prior to this call via set_design or linearize_at. The value
of the residual is returned in \a RHS.*/
PetscErrorCode IP_SSAHardavForwardProblem::assemble_residual(IceModelVec2V &u, IceModelVec2V &RHS) {
  PetscErrorCode ierr;

  Vector2 **u_a, **rhs_a;

  ierr = u.get_array(u_a); CHKERRQ(ierr);
  ierr = RHS.get_array(rhs_a); CHKERRQ(ierr);

  DMDALocalInfo *info = NULL;
  ierr = this->compute_local_function(info, const_cast<const Vector2 **>(u_a), rhs_a); CHKERRQ(ierr);

  ierr = u.end_access(); CHKERRQ(ierr);
  ierr = RHS.end_access(); CHKERRQ(ierr);

  return 0;
}

//! Computes the residual function \f$\mathcal{R}(u, \zeta)\f$ defined in the class-level documentation.
/* The return value is specified via a Vec for the benefit of certain TAO routines.  Otherwise,
the method is identical to the assemble_residual returning values as a StateVec (an IceModelVec2V).*/
PetscErrorCode IP_SSAHardavForwardProblem::assemble_residual(IceModelVec2V &u, Vec RHS) {
  PetscErrorCode ierr;

  Vector2 **u_a, **rhs_a;

  ierr = u.get_array(u_a); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(SSADA->get(), RHS, &rhs_a); CHKERRQ(ierr);

  DMDALocalInfo *info = NULL;
  ierr = this->compute_local_function(info, const_cast<const Vector2 **>(u_a), rhs_a); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(SSADA->get(), RHS, &rhs_a); CHKERRQ(ierr);
  ierr = u.end_access(); CHKERRQ(ierr);

  return 0;
}

//! Assembles the state Jacobian matrix.
/* The matrix depends on the current value of the design variable \f$\zeta\f$ and the current
value of the state variable \f$u\f$.  The specification of \f$\zeta\f$ is done earlier
with set_design or linearize_at.  The value of \f$u\f$ is specified explicitly as an argument
to this method.
  @param[in] u Current state variable value.
  @param[out] J computed state Jacobian.
*/
PetscErrorCode IP_SSAHardavForwardProblem::assemble_jacobian_state(IceModelVec2V &u, Mat Jac) {
  PetscErrorCode ierr;

  Vector2 **u_a;
  ierr = u.get_array(u_a); CHKERRQ(ierr);

  DMDALocalInfo *info = NULL;
  ierr = this->compute_local_jacobian(info, const_cast<const Vector2 **>(u_a), Jac); CHKERRQ(ierr);

  ierr = u.end_access(); CHKERRQ(ierr);

  return 0;
}

//! Applies the design Jacobian matrix to a perturbation of the design variable.
/*! The return value uses a DesignVector (IceModelVec2V), which can be ghostless. Ghosts (if present) are updated.
\overload
*/
PetscErrorCode IP_SSAHardavForwardProblem::apply_jacobian_design(IceModelVec2V &u, IceModelVec2S &dzeta, IceModelVec2V &du) {
  PetscErrorCode ierr;
  Vector2 **du_a;
  ierr = du.get_array(du_a); CHKERRQ(ierr);
  ierr = this->apply_jacobian_design(u, dzeta, du_a);
  ierr = du.end_access(); CHKERRQ(ierr);
  return 0;
}

//! Applies the design Jacobian matrix to a perturbation of the design variable.
/*! The return value is a Vec for the benefit of TAO. It is assumed to be ghostless; no communication is done.
\overload
*/
PetscErrorCode IP_SSAHardavForwardProblem::apply_jacobian_design(IceModelVec2V &u, IceModelVec2S &dzeta, Vec du) {
  PetscErrorCode ierr;
  Vector2 **du_a;
  ierr = DMDAVecGetArray(SSADA->get(), du, &du_a); CHKERRQ(ierr);
  ierr = this->apply_jacobian_design(u, dzeta, du_a);
  ierr = DMDAVecRestoreArray(SSADA->get(), du, &du_a); CHKERRQ(ierr);
  return 0;
}

//! @brief Applies the design Jacobian matrix to a perturbation of the
//! design variable.

/*! The matrix depends on the current value of the design variable
    \f$\zeta\f$ and the current value of the state variable \f$u\f$.
    The specification of \f$\zeta\f$ is done earlier with set_design
    or linearize_at. The value of \f$u\f$ is specified explicitly as
    an argument to this method.

  @param[in] u Current state variable value.

  @param[in] dzeta Perturbation of the design variable. Prefers
                   vectors with ghosts; will copy to a ghosted vector
                   if needed.

  @param[out] du_a Computed corresponding perturbation of the state
                   variable. The array \a du_a should be extracted
                   first from a Vec or an IceModelVec.

  Typically this method is called via one of its overloads.
*/
PetscErrorCode IP_SSAHardavForwardProblem::apply_jacobian_design(IceModelVec2V &u,
                                                                 IceModelVec2S &dzeta,
                                                                 Vector2 **du_a) {
  PetscErrorCode ierr;

  ierr = m_zeta->begin_access(); CHKERRQ(ierr);

  ierr = u.begin_access(); CHKERRQ(ierr);

  IceModelVec2S *dzeta_local;
  if (dzeta.has_ghosts()) {
    dzeta_local = &dzeta;
  } else {
    ierr = m_dzeta_local.copy_from(dzeta); CHKERRQ(ierr);
    dzeta_local = &m_dzeta_local;
  }
  ierr = dzeta_local->begin_access(); CHKERRQ(ierr);

  // Zero out the portion of the function we are responsible for computing.
  for (int i =grid.xs; i<grid.xs+grid.xm; i++) {
    for (int j =grid.ys; j<grid.ys+grid.ym; j++) {
      du_a[i][j].u = 0.0;
      du_a[i][j].v = 0.0;
    }
  }

  // Aliases to help with notation consistency below.
  IceModelVec2Int *m_dirichletLocations = bc_locations;
  IceModelVec2V   *m_dirichletValues    = m_vel_bc;
  double           m_dirichletWeight    = m_dirichletScale;

  Vector2 u_e[FEQuadrature::Nk];
  Vector2 u_q[FEQuadrature::Nq];
  double Du_q[FEQuadrature::Nq][3];

  Vector2 du_e[FEQuadrature::Nk];

  double dzeta_e[FEQuadrature::Nk];

  double zeta_e[FEQuadrature::Nk];

  double dB_e[FEQuadrature::Nk];
  double dB_q[FEQuadrature::Nq];

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  DirichletData_Vector dirichletBC;
  ierr = dirichletBC.init(m_dirichletLocations, m_dirichletValues,
                          m_dirichletWeight); CHKERRQ(ierr);
  DirichletData_Scalar fixedZeta;
  ierr = fixedZeta.init(m_fixed_design_locations, NULL);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  // Loop through all elements.
  int xs = m_element_index.xs,
    xm   = m_element_index.xm,
    ys   = m_element_index.ys,
    ym   = m_element_index.ym;

  for (int i =xs; i<xs+xm; i++) {
    for (int j =ys; j<ys+ym; j++) {

      // Zero out the element-local residual in prep for updating it.
      for (unsigned int k=0; k<FEQuadrature::Nk; k++) {
        du_e[k].u = 0;
        du_e[k].v = 0;
      }

      // Index into coefficient storage in m_coefficients
      const int ij = m_element_index.flatten(i, j);

      // Initialize the map from global to local degrees of freedom for this element.
      m_dofmap.reset(i, j, m_grid);

      // Obtain the value of the solution at the nodes adjacent to the element,
      // fix dirichlet values, and compute values at quad pts.
      m_dofmap.extractLocalDOFs(i, j, u, u_e);
      if (dirichletBC) {
        dirichletBC.constrain(m_dofmap);
        dirichletBC.update(m_dofmap, u_e);
      }
      m_quadrature_vector.computeTrialFunctionValues(u_e, u_q, Du_q);

      // Compute dzeta at the nodes
      m_dofmap.extractLocalDOFs(i, j, *dzeta_local, dzeta_e);
      if (fixedZeta) fixedZeta.update_homogeneous(m_dofmap, dzeta_e);

      // Compute the change in hardav with respect to zeta at the quad points.
      m_dofmap.extractLocalDOFs(i, j, *m_zeta, zeta_e);
      for (unsigned int k=0; k<FEQuadrature::Nk; k++) {
        m_design_param.toDesignVariable(zeta_e[k], NULL, dB_e + k);
        dB_e[k]*=dzeta_e[k];
      }
      m_quadrature.computeTrialFunctionValues(dB_e, dB_q);

      for (unsigned int q=0; q<FEQuadrature::Nq; q++) {
        // Symmetric gradient at the quadrature point.
        double *Duqq = Du_q[q];

        const SSACoefficients *coefficients = &m_coefficients[ij*FEQuadrature::Nq+q];

        double d_nuH = 0;
        if (coefficients->H >= strength_extension->get_min_thickness()) {
          flow_law->effective_viscosity(dB_q[q], secondInvariantDu_2D(Duqq), &d_nuH, NULL);
          d_nuH  *= (2*coefficients->H);
        }

        for (unsigned int k=0; k<FEQuadrature::Nk; k++) {
          const FEFunctionGerm &testqk = test[q][k];
          du_e[k].u += JxW[q]*d_nuH*(testqk.dx*(2*Duqq[0]+Duqq[1]) + testqk.dy*Duqq[2]);
          du_e[k].v += JxW[q]*d_nuH*(testqk.dy*(2*Duqq[1]+Duqq[0]) + testqk.dx*Duqq[2]);
        }
      } // q
      m_dofmap.addLocalResidualBlock(du_e, du_a);
    } // j
  } // i

  if (dirichletBC) dirichletBC.fix_residual_homogeneous(du_a);

  ierr = dirichletBC.finish(); CHKERRQ(ierr);
  ierr = fixedZeta.finish(); CHKERRQ(ierr);

  ierr = m_zeta->end_access(); CHKERRQ(ierr);
  ierr = u.end_access(); CHKERRQ(ierr);
  ierr = dzeta_local->end_access(); CHKERRQ(ierr);

  return 0;
}

//! Applies the transpose of the design Jacobian matrix to a perturbation of the state variable.
/*! The return value uses a StateVector (IceModelVec2S) which can be ghostless; ghosts (if present) are updated.
\overload
*/
PetscErrorCode IP_SSAHardavForwardProblem::apply_jacobian_design_transpose(IceModelVec2V &u,
                                                                           IceModelVec2V &du,
                                                                           IceModelVec2S &dzeta) {
  PetscErrorCode ierr;
  double **dzeta_a;
  ierr = dzeta.get_array(dzeta_a); CHKERRQ(ierr);
  ierr = this->apply_jacobian_design_transpose(u, du, dzeta_a); CHKERRQ(ierr);
  ierr = dzeta.end_access(); CHKERRQ(ierr);
  return 0;
}

//! Applies the transpose of the design Jacobian matrix to a perturbation of the state variable.
/*! The return value uses a Vec for the benefit of TAO.  It is assumed to be ghostless; no communication is done.
\overload */
PetscErrorCode IP_SSAHardavForwardProblem::apply_jacobian_design_transpose(IceModelVec2V &u,
                                                                           IceModelVec2V &du,
                                                                           Vec dzeta) {
  PetscErrorCode ierr;
  double **dzeta_a;
  PISMDM::Ptr da2;
  ierr = m_grid.get_dm(1, m_grid.max_stencil_width, da2); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(da2->get(), dzeta, &dzeta_a); CHKERRQ(ierr);
  ierr = this->apply_jacobian_design_transpose(u, du, dzeta_a); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArray(da2->get(), dzeta, &dzeta_a); CHKERRQ(ierr);
  return 0;
}

//! @brief Applies the transpose of the design Jacobian matrix to a
//! perturbation of the state variable.

/*! The matrix depends on the current value of the design variable
    \f$\zeta\f$ and the current value of the state variable \f$u\f$.
    The specification of \f$\zeta\f$ is done earlier with set_design
    or linearize_at. The value of \f$u\f$ is specified explicitly as
    an argument to this method.

  @param[in] u Current state variable value.

  @param[in] du Perturbation of the state variable. Prefers vectors
                with ghosts; will copy to a ghosted vector if need be.

  @param[out] dzeta_a Computed corresponding perturbation of the
                      design variable. The array \a dzeta_a should be
                      extracted first from a Vec or an IceModelVec.

  Typically this method is called via one of its overloads.
*/
PetscErrorCode IP_SSAHardavForwardProblem::apply_jacobian_design_transpose(IceModelVec2V &u,
                                                                           IceModelVec2V &du,
                                                                           double **dzeta_a) {
  PetscErrorCode ierr;

  ierr = m_zeta->begin_access(); CHKERRQ(ierr);

  ierr = u.begin_access(); CHKERRQ(ierr);

  IceModelVec2V *du_local;
  if (du.has_ghosts()) {
    du_local = &du;
  } else {
    ierr = m_du_local.copy_from(du); CHKERRQ(ierr);
    du_local = &m_du_local;
  }
  ierr = du_local->begin_access(); CHKERRQ(ierr);

  Vector2 u_e[FEQuadrature::Nk];
  Vector2 u_q[FEQuadrature::Nq];
  double Du_q[FEQuadrature::Nq][3];

  Vector2 du_e[FEQuadrature::Nk];
  Vector2 du_q[FEQuadrature::Nq];
  Vector2 du_dx_q[FEQuadrature::Nq];
  Vector2 du_dy_q[FEQuadrature::Nq];

  double dzeta_e[FEQuadrature::Nk];

  // An Nq by Nk array of test function values.
  const FEFunctionGerm (*test)[FEQuadrature::Nk] = m_quadrature.testFunctionValues();

  DirichletData_Vector dirichletBC;
  // Aliases to help with notation consistency.
  IceModelVec2Int *m_dirichletLocations = bc_locations;
  IceModelVec2V   *m_dirichletValues = m_vel_bc;
  double        m_dirichletWeight = m_dirichletScale;
  ierr = dirichletBC.init(m_dirichletLocations, m_dirichletValues,
                          m_dirichletWeight); CHKERRQ(ierr);

  // Jacobian times weights for quadrature.
  const double* JxW = m_quadrature.getWeightedJacobian();

  // Zero out the portion of the function we are responsible for computing.
  for (int i = grid.xs; i < grid.xs + grid.xm; i++) {
    for (int j = grid.ys; j < grid.ys + grid.ym; j++) {
      dzeta_a[i][j] = 0;
    }
  }

  int xs = m_element_index.xs, xm = m_element_index.xm,
           ys = m_element_index.ys, ym = m_element_index.ym;
  for (int i = xs; i < xs + xm; i++) {
    for (int j = ys; j < ys + ym; j++) {
      // Index into coefficient storage in m_coefficients
      const int ij = m_element_index.flatten(i, j);

      // Initialize the map from global to local degrees of freedom for this element.
      m_dofmap.reset(i, j, m_grid);

      // Obtain the value of the solution at the nodes adjacent to the element.
      // Compute the solution values and symmetric gradient at the quadrature points.
      m_dofmap.extractLocalDOFs(i, j, du, du_e);
      if (dirichletBC) dirichletBC.update_homogeneous(m_dofmap, du_e);
      m_quadrature_vector.computeTrialFunctionValues(du_e, du_q, du_dx_q, du_dy_q);

      m_dofmap.extractLocalDOFs(i, j, u, u_e);
      if (dirichletBC) dirichletBC.update(m_dofmap, u_e);
      m_quadrature_vector.computeTrialFunctionValues(u_e, u_q, Du_q);

      // Zero out the element - local residual in prep for updating it.
      for (unsigned int k = 0; k < FEQuadrature::Nk; k++) {
        dzeta_e[k] = 0;
      }

      for (unsigned int q = 0; q < FEQuadrature::Nq; q++) {
        // Symmetric gradient at the quadrature point.
        double *Duqq = Du_q[q];

        const SSACoefficients *coefficients = &m_coefficients[ij*FEQuadrature::Nq + q];

        // Determine "d_nuH / dB" at the quadrature point
        double d_nuH_dB = 0;
        if (coefficients->H >= strength_extension->get_min_thickness()) {
          flow_law->effective_viscosity(1., secondInvariantDu_2D(Duqq), &d_nuH_dB, NULL);
          d_nuH_dB *= (2*coefficients->H);
        }

        for (unsigned int k = 0; k < FEQuadrature::Nk; k++) {
          dzeta_e[k] += JxW[q]*d_nuH_dB*test[q][k].val*((du_dx_q[q].u*(2*Duqq[0] + Duqq[1]) +
                                                         du_dy_q[q].u*Duqq[2]) +
                                                        (du_dy_q[q].v*(2*Duqq[1] + Duqq[0]) +
                                                         du_dx_q[q].v*Duqq[2]));
        }
      } // q

      m_dofmap.addLocalResidualBlock(dzeta_e, dzeta_a);
    } // j
  } // i
  ierr = dirichletBC.finish(); CHKERRQ(ierr);

  for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
    for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
      double dB_dzeta;
      m_design_param.toDesignVariable((*m_zeta)(i, j), NULL, &dB_dzeta);
      dzeta_a[i][j] *= dB_dzeta;
    }
  }

  if (m_fixed_design_locations) {
    DirichletData_Scalar fixedZeta;
    ierr = fixedZeta.init(m_fixed_design_locations, NULL); CHKERRQ(ierr);
    fixedZeta.fix_residual_homogeneous(dzeta_a);
    ierr = fixedZeta.finish(); CHKERRQ(ierr);
  }

  ierr = m_zeta->end_access(); CHKERRQ(ierr);
  ierr = u.end_access(); CHKERRQ(ierr);
  ierr = du_local->end_access(); CHKERRQ(ierr);

  return 0;
}

/*!\brief Applies the linearization of the forward map (i.e. the reduced gradient \f$DF\f$ described in
the class-level documentation.) */
/*! As described previously,
\f[
Df = J_{\rm State}^{-1} J_{\rm Design}.
\f]
Applying the linearization then involves the solution of a linear equation.
The matrices \f$J_{\rm State}\f$ and \f$J_{\rm Design}\f$ both depend on the value of the
design variable \f$\zeta\f$ and the value of the corresponding state variable \f$u=F(\zeta)\f$.
These are established by first calling linearize_at.
  @param[in]   dzeta     Perturbation of the design variable
  @param[out]  du        Computed corresponding perturbation of the state variable; ghosts (if present) are updated.
*/
PetscErrorCode IP_SSAHardavForwardProblem::apply_linearization(IceModelVec2S &dzeta, IceModelVec2V &du) {

  PetscErrorCode ierr;

  if (m_rebuild_J_state) {
    ierr = this->assemble_jacobian_state(m_velocity, m_J_state); CHKERRQ(ierr);
    m_rebuild_J_state = false;
  }

  ierr = this->apply_jacobian_design(m_velocity, dzeta, m_du_global); CHKERRQ(ierr);
  ierr = m_du_global.scale(-1); CHKERRQ(ierr);

  // call PETSc to solve linear system by iterative method.
  ierr = KSPSetOperators(m_ksp, m_J_state, m_J_state, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_ksp, m_du_global.get_vec(), m_du_global.get_vec()); CHKERRQ(ierr); // SOLVE

  KSPConvergedReason  reason;
  ierr = KSPGetConvergedReason(m_ksp, &reason); CHKERRQ(ierr);
  if (reason < 0) {
    SETERRQ1(grid.com, 1,
      "IP_SSAHardavForwardProblem::apply_linearization solve failed to converge (KSP reason %s)\n\n", KSPConvergedReasons[reason]);
  }
  else
  {
    verbPrintf(4, grid.com, "IP_SSAHardavForwardProblem::apply_linearization converged (KSP reason %s)\n", KSPConvergedReasons[reason]);
  }

  ierr = du.copy_from(m_du_global); CHKERRQ(ierr);
  return 0;
}

/*! \brief Applies the transpose of the linearization of the forward map
 (i.e. the transpose of the reduced gradient \f$DF\f$ described in the class-level documentation.) */
/*!  As described previously,
\f[
Df = J_{\rm State}^{-1} J_{\rm Design}.
\f]
so
\f[
Df^t = J_{\rm Design}^t \; (J_{\rm State}^t)^{-1} .
\f]
Applying the transpose of the linearization then involves the solution of a linear equation.
The matrices \f$J_{\rm State}\f$ and \f$J_{\rm Design}\f$ both depend on the value of the
design variable \f$\zeta\f$ and the value of the corresponding state variable \f$u=F(\zeta)\f$.
These are established by first calling linearize_at.
  @param[in]   du     Perturbation of the state variable
  @param[out]  dzeta  Computed corresponding perturbation of the design variable; ghosts (if present) are updated.
*/
PetscErrorCode IP_SSAHardavForwardProblem::apply_linearization_transpose(IceModelVec2V &du,
                                                                         IceModelVec2S &dzeta) {

  PetscErrorCode ierr;

  if (m_rebuild_J_state) {
    ierr = this->assemble_jacobian_state(m_velocity, m_J_state); CHKERRQ(ierr);
    m_rebuild_J_state = false;
  }

  // Aliases to help with notation consistency below.
  IceModelVec2Int *m_dirichletLocations = bc_locations;
  IceModelVec2V   *m_dirichletValues    = m_vel_bc;
  double        m_dirichletWeight    = m_dirichletScale;

  ierr = m_du_global.copy_from(du); CHKERRQ(ierr);
  Vector2 **du_a;
  ierr = m_du_global.get_array(du_a); CHKERRQ(ierr);
  DirichletData_Vector dirichletBC;
  ierr = dirichletBC.init(m_dirichletLocations, m_dirichletValues, m_dirichletWeight); CHKERRQ(ierr);
  if (dirichletBC) dirichletBC.fix_residual_homogeneous(du_a);
  ierr = dirichletBC.finish(); CHKERRQ(ierr);
  ierr = m_du_global.end_access(); CHKERRQ(ierr);

  // call PETSc to solve linear system by iterative method.
  ierr = KSPSetOperators(m_ksp, m_J_state, m_J_state, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
  ierr = KSPSolve(m_ksp, m_du_global.get_vec(), m_du_global.get_vec()); CHKERRQ(ierr); // SOLVE

  KSPConvergedReason  reason;
  ierr = KSPGetConvergedReason(m_ksp, &reason); CHKERRQ(ierr);
  if (reason < 0) {
    SETERRQ1(grid.com, 1,
      "IP_SSAHardavForwardProblem::apply_linearization solve failed to converge (KSP reason %s)\n\n", KSPConvergedReasons[reason]);
  }
  else
  {
    verbPrintf(4, grid.com, "IP_SSAHardavForwardProblem::apply_linearization converged (KSP reason %s)\n", KSPConvergedReasons[reason]);
  }

  ierr = this->apply_jacobian_design_transpose(m_velocity, m_du_global, dzeta); CHKERRQ(ierr);
  ierr = dzeta.scale(-1); CHKERRQ(ierr);

  if (dzeta.has_ghosts()) {
    ierr = dzeta.update_ghosts();
  }

  return 0;
}

} // end of namespace pism
