// Copyright (C) 2026 Andy Aschwanden and Constantine Khroulev
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

#include "pism/inverse/IP_BlatterHardavForwardProblem.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Vars.hh"
#include "pism/util/error_handling.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/util/petscwrappers/DM.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/Logger.hh"

namespace pism {
namespace inverse {

IP_BlatterHardavForwardProblem::IP_BlatterHardavForwardProblem(
    std::shared_ptr<const Grid> grid,
    int Mz,
    int coarsening_factor,
    IPDesignVariableParameterization &tp)
  : Blatter(grid, Mz, coarsening_factor),
    m_zeta(NULL),
    m_hardav(m_grid, "hardav"),
    m_fixed_design_locations(NULL),
    m_design_param(tp),
    m_rebuild_J_state(true)
{
  m_surface_velocity.reset(new array::Vector(m_grid, "surface_velocity"));
  m_surface_velocity->metadata(0)
    .long_name("x-component of surface velocity from Blatter solver")
    .units("m s^-1")
    .output_units("m year^-1");
  m_surface_velocity->metadata(1)
    .long_name("y-component of surface velocity from Blatter solver")
    .units("m s^-1")
    .output_units("m year^-1");

  m_hardav.metadata(0)
    .long_name("vertically-averaged ice hardness (Glen B = A^{-1/n})")
    .units("Pa s^(1/3)");

  // Zero-initialize PETSc wrapper members. The Wrapper<T> class does NOT
  // zero-initialize m_value, so rawptr() always returns non-null (it's
  // &m_value), and get() returns garbage. We must explicitly null them.
  *m_J_picard.rawptr() = nullptr;
  *m_ksp.rawptr() = nullptr;
}

void IP_BlatterHardavForwardProblem::init() {
  init_impl();
}

//! Sets the current value of the design parameter `zeta`. Computes hardav
//! via the parameterization and caches it. The hardness is *not* yet pushed
//! into the Blatter solver's internal 3D hardness array; that happens inside
//! linearize_at after Blatter::update() runs (see install_hardav_on_sigma_grid).
void IP_BlatterHardavForwardProblem::set_design(array::Scalar &new_zeta) {

  m_zeta = &new_zeta;

  // Convert zeta to hardav.
  m_design_param.convertToDesignVariable(*m_zeta, m_hardav);

  m_rebuild_J_state = true;
}

//! Installs the column-constant hardness into the Blatter sigma grid's 3D
//! hardness storage. Called after Blatter::update() has run the SNES with
//! the enthalpy-derived hardness — we then overwrite with our hardav and
//! mark the Jacobian for rebuild so subsequent adjoint solves use the right
//! viscosity.
//!
//! TODO: this is a placeholder. Blatter::init_ice_hardness (Blatter.cc:588)
//! is not virtual and the 3D hardness array it populates is a local
//! `DataAccess<double***>` inside that function rather than a persistent
//! member. To make hardav-as-design actually work for the Blatter solver,
//! one of two changes is required upstream:
//!
//!   (a) make `Blatter::init_ice_hardness` virtual, lift the
//!       DataAccess<double***> hardness storage to a member, and override
//!       in this class to use m_hardav instead of inputs.enthalpy.
//!
//!   (b) add a `set_ice_hardness_override(const array::Scalar *)` hook to
//!       Blatter that, when set, replaces the enthalpy-derived hardness
//!       with a column-constant value from the provided 2D field.
//!
//! Either is a small change to Blatter.{hh,cc} but is intentionally not done
//! here so this PR can land separately. Until then, linearize_at will throw.
void IP_BlatterHardavForwardProblem::install_hardav_on_sigma_grid() {
  throw RuntimeError(PISM_ERROR_LOCATION,
                     "IP_BlatterHardavForwardProblem::install_hardav_on_sigma_grid "
                     "is not yet implemented. Blatter::init_ice_hardness needs to "
                     "be made overridable (or a hardness-override hook added to "
                     "Blatter) before hardness inversion can be exercised. See "
                     "the comment in IP_BlatterHardavForwardProblem.cc for the "
                     "two minimal options.");
}

//! Sets `zeta`, runs the Blatter SNES with our hardness, and extracts the
//! surface velocity.
std::shared_ptr<TerminationReason> IP_BlatterHardavForwardProblem::linearize_at(
    array::Scalar &zeta) {

  this->set_design(zeta);

  // Build the standard Inputs and run a forward solve, identical to the tauc
  // forward problem. The hardness needs to be swapped *after* this for the
  // adjoint to be consistent — see install_hardav_on_sigma_grid().
  const auto &variables = m_grid->variables();

  Geometry geometry(m_grid);
  geometry.ice_thickness.copy_from(*variables.get_2d_scalar("land_ice_thickness"));
  geometry.bed_elevation.copy_from(*variables.get_2d_scalar("bedrock_altitude"));

  if (variables.is_available("sea_surface_height_above_reference_ellipsoid")) {
    geometry.sea_level_elevation.copy_from(
        *variables.get_2d_scalar("sea_surface_height_above_reference_ellipsoid"));
  } else {
    geometry.sea_level_elevation.copy_from(geometry.bed_elevation);
    geometry.sea_level_elevation.shift(-1000.0);
  }

  if (m_config->get_flag("geometry.part_grid.enabled") &&
      variables.is_available("ice_area_specific_volume")) {
    geometry.ice_area_specific_volume.copy_from(
        *variables.get_2d_scalar("ice_area_specific_volume"));
  } else {
    geometry.ice_area_specific_volume.set(0.0);
  }

  geometry.ensure_consistency(
      m_config->get_number("stress_balance.ice_free_thickness_standard"));

  stressbalance::Inputs inputs;
  inputs.geometry           = &geometry;
  inputs.basal_yield_stress = variables.get_2d_scalar("tauc");
  inputs.enthalpy           = variables.get_3d_scalar("enthalpy");
  inputs.basal_melt_rate    = NULL;
  inputs.age                = NULL;
  inputs.water_column_pressure = NULL;

  if (variables.is_available("vel_bc_mask")) {
    inputs.bc_mask = variables.get_2d_scalar("vel_bc_mask");
  }
  if (variables.is_available("vel_bc")) {
    inputs.bc_values = variables.get_2d_vector("vel_bc");
  }

  this->update(inputs, true);

  // Install our column-constant hardness on the sigma grid and mark the
  // Jacobian for rebuild. (See note above: this currently throws because
  // Blatter doesn't yet expose the right hook.)
  this->install_hardav_on_sigma_grid();

  this->extract_surface_velocity();

  SNESConvergedReason snes_reason;
  PetscErrorCode ierr = SNESGetConvergedReason(m_snes, &snes_reason);
  PISM_CHK(ierr, "SNESGetConvergedReason");

  if (snes_reason > 0) {
    return std::shared_ptr<TerminationReason>(
      new GenericTerminationReason(1, "Blatter solve converged"));
  } else {
    return std::shared_ptr<TerminationReason>(
      new GenericTerminationReason(-1, "Blatter solve failed to converge"));
  }
}

//! Extract the 2D surface velocity from the 3D Blatter solution. Identical
//! to the tauc version.
void IP_BlatterHardavForwardProblem::extract_surface_velocity() {
  PetscErrorCode ierr;
  Vector2d ***x = nullptr;

  ierr = DMDAVecGetArrayRead(m_da, m_x, &x);
  PISM_CHK(ierr, "DMDAVecGetArrayRead");

  int Mz = (int)m_u_sigma->levels().size();

  array::AccessScope list{m_surface_velocity.get()};

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();
    (*m_surface_velocity)(i, j) = x[j][i][Mz - 1]; // STORAGE_ORDER
  }

  ierr = DMDAVecRestoreArrayRead(m_da, m_x, &x);
  PISM_CHK(ierr, "DMDAVecRestoreArrayRead");
}

petsc::DM &IP_BlatterHardavForwardProblem::get_da() const {
  return *m_grid->get_dm(1, m_config->get_number("grid.max_stencil_width"));
}

// ============================================================================
// Design Jacobian — volume integral analog of IP_SSAHardavForwardProblem
// ============================================================================
//
// For hardness, the design Jacobian is supported on *every* sigma element,
// not just the basal face. The 3D analog of the SSA viscous Jacobian is:
//
//   J_design[k] = d/dB [ ∫ eta(B, |Du|^2) D(u) : D(psi_k) dV ]
//               = ∫ (d_eta/d_B)(B, |Du|^2) (D(u) : D(psi_k)) dV
//
// where D(u) is the 3D symmetric strain rate and eta is the effective
// viscosity from m_flow_law->effective_viscosity. The SSA version (lines
// 359-373 of IP_SSAHardavForwardProblem.cc) uses the same call with a 2D
// strain rate; we'd loop over sigma elements and quadrature points and use
// the 3D second invariant.
//
// Implementing this correctly requires reading Blatter's residual assembly
// in src/stressbalance/blatter/residual.cc to reuse the strain-rate
// machinery rather than reimplementing it. That work is deferred.

void IP_BlatterHardavForwardProblem::apply_jacobian_design(
    array::Vector &u, array::Scalar &dzeta, array::Vector &du) {
  (void)u;
  (void)dzeta;
  du.set(0.0);
  throw RuntimeError(PISM_ERROR_LOCATION,
                     "IP_BlatterHardavForwardProblem::apply_jacobian_design "
                     "is not yet implemented (volume Jacobian for hardness).");
}

void IP_BlatterHardavForwardProblem::apply_jacobian_design_transpose(
    array::Vector &u, array::Vector &du, array::Scalar &dzeta) {
  (void)u;
  (void)du;
  dzeta.set(0.0);
  throw RuntimeError(PISM_ERROR_LOCATION,
                     "IP_BlatterHardavForwardProblem::apply_jacobian_design_transpose "
                     "is not yet implemented (volume Jacobian transpose for hardness).");
}

void IP_BlatterHardavForwardProblem::apply_jacobian_design_3d(
    array::Scalar &dzeta, Vec result_3d) {
  (void)dzeta;
  (void)result_3d;
  throw RuntimeError(PISM_ERROR_LOCATION,
                     "IP_BlatterHardavForwardProblem::apply_jacobian_design_3d "
                     "is not yet implemented.");
}

void IP_BlatterHardavForwardProblem::apply_jacobian_design_transpose_3d(
    Vec lambda_3d, array::Scalar &dzeta) {
  (void)lambda_3d;
  (void)dzeta;
  throw RuntimeError(PISM_ERROR_LOCATION,
                     "IP_BlatterHardavForwardProblem::apply_jacobian_design_transpose_3d "
                     "is not yet implemented.");
}

void IP_BlatterHardavForwardProblem::apply_linearization(
    array::Scalar &dzeta, array::Vector &du) {
  (void)dzeta;
  (void)du;
  throw RuntimeError(PISM_ERROR_LOCATION,
                     "IP_BlatterHardavForwardProblem::apply_linearization "
                     "is not yet implemented (depends on apply_jacobian_design_3d).");
}

void IP_BlatterHardavForwardProblem::apply_linearization_transpose(
    array::Vector &du, array::Scalar &dzeta) {
  (void)du;
  (void)dzeta;
  throw RuntimeError(PISM_ERROR_LOCATION,
                     "IP_BlatterHardavForwardProblem::apply_linearization_transpose "
                     "is not yet implemented (depends on apply_jacobian_design_transpose_3d).");
}

} // end of namespace inverse
} // end of namespace pism
