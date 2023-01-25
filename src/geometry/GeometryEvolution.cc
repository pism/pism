/* Copyright (C) 2016--2022 PISM Authors
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

#include "GeometryEvolution.hh"

#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Staggered.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Mask.hh"

#include "pism/geometry/part_grid_threshold_thickness.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/Context.hh"

#include "flux_limiter.hh"

namespace pism {

using mask::floating_ice;
using mask::grounded_ice;
using mask::ice_free;
using mask::ice_free_land;
using mask::ice_free_ocean;
using mask::icy;

struct GeometryEvolution::Impl {
  Impl(IceGrid::ConstPtr g);

  const Profiling &profile;

  GeometryCalculator gc;

  double ice_density;

  //! True if the basal melt rate contributes to geometry evolution.
  bool use_bmr;

  //! True if the part-grid scheme is enabled.
  bool use_part_grid;

  //! Flux divergence (used to track thickness changes due to flow).
  array::Scalar flux_divergence;

  //! Conservation error due to enforcing non-negativity of ice thickness and enforcing
  //! thickness BC.
  array::Scalar conservation_error;

  //! Effective surface mass balance.
  array::Scalar effective_SMB;

  //! Effective basal mass balance.
  array::Scalar effective_BMB;

  //! Change in ice thickness due to flow during the last time step.
  array::Scalar thickness_change;

  //! Change in the ice area-specific volume due to flow during the last time step.
  array::Scalar ice_area_specific_volume_change;

  //! Flux through cell interfaces. Ghosted.
  array::Staggered1 flux_staggered;

  // Work space
  array::Vector1   input_velocity; // a ghosted copy; not modified
  array::Scalar1   bed_elevation; // a copy; not modified
  array::Scalar1   sea_level;   // a copy; not modified
  array::Scalar1   ice_thickness; // updated in place
  array::Scalar1   area_specific_volume; // updated in place
  array::Scalar1   surface_elevation; // updated to maintain consistency
  array::CellType1 cell_type;   // updated to maintain consistency
  array::Scalar1   residual;    // temporary storage
  array::Scalar1   thickness;   // temporary storage
};

GeometryEvolution::Impl::Impl(IceGrid::ConstPtr grid)
  : profile(grid->ctx()->profiling()),
    gc(*grid->ctx()->config()),
    flux_divergence(grid, "flux_divergence"),
    conservation_error(grid, "conservation_error"),
    effective_SMB(grid, "effective_SMB"),
    effective_BMB(grid, "effective_BMB"),
    thickness_change(grid, "thickness_change"),
    ice_area_specific_volume_change(grid, "ice_area_specific_volume_change"),
    flux_staggered(grid, "flux_staggered"),
    input_velocity(grid, "input_velocity"),
    bed_elevation(grid, "bed_elevation"),
    sea_level(grid, "sea_level"),
    ice_thickness(grid, "ice_thickness"),
    area_specific_volume(grid, "area_specific_volume"),
    surface_elevation(grid, "surface_elevation"),
    cell_type(grid, "cell_type"),
    residual(grid, "residual"),
    thickness(grid, "thickness") {

  Config::ConstPtr config = grid->ctx()->config();

  gc.set_icefree_thickness(config->get_number("geometry.ice_free_thickness_standard"));

  // constants
  {
    ice_density   = config->get_number("constants.ice.density");
    use_bmr       = config->get_flag("geometry.update.use_basal_melt_rate");
    use_part_grid = config->get_flag("geometry.part_grid.enabled");
  }

  // reported quantities
  {
    // This is the only reported field that is ghosted (we need ghosts to compute flux divergence).
    flux_staggered.set_attrs("diagnostic", "fluxes through cell interfaces (sides)"
                             " on the staggered grid (x-offset)",
                             "m2 s-1", "m2 year-1", "", 0);
    flux_staggered.set_attrs("diagnostic", "fluxes through cell interfaces (sides)"
                             " on the staggered grid (y-offset)",
                             "m2 s-1", "m2 year-1", "", 1);

    flux_divergence.set_attrs("diagnostic", "flux divergence", "m s-1", "m year-1", "", 0);

    conservation_error.set_attrs("diagnostic",
                                 "conservation error due to enforcing non-negativity of"
                                 " ice thickness (over the last time step)",
                                 "meters", "meters", "", 0);

    effective_SMB.set_attrs("internal", "effective surface mass balance over the last time step",
                            "meters", "meters", "", 0);

    effective_BMB.set_attrs("internal", "effective basal mass balance over the last time step",
                            "meters", "meters", "", 0);

    thickness_change.set_attrs("internal", "change in thickness due to flow",
                               "meters", "meters", "", 0);

    ice_area_specific_volume_change.set_attrs("interval",
                                              "change in area-specific volume due to flow",
                                              "meters3 / meters2", "meters3 / meters2", "", 0);
  }

  // internal storage
  {
    input_velocity.set_attrs("internal", "ghosted copy of the input velocity",
                             "meters / second", "meters / second", "", 0);

    bed_elevation.set_attrs("internal", "ghosted copy of the bed elevation",
                            "meters", "meters", "", 0);

    sea_level.set_attrs("internal", "ghosted copy of the sea level elevation",
                        "meters", "meters", "", 0);

    ice_thickness.set_attrs("internal", "working (ghosted) copy of the ice thickness",
                            "meters", "meters", "", 0);

    area_specific_volume.set_attrs("internal", "working (ghosted) copy of the area specific volume",
                                   "meters3 / meters2", "meters3 / meters2", "", 0);

    surface_elevation.set_attrs("internal", "working (ghosted) copy of the surface elevation",
                                "meters", "meters", "", 0);

    cell_type.set_attrs("internal", "working (ghosted) copy of the cell type mask",
                        "", "", "", 0);

    residual.set_attrs("internal", "residual area specific volume",
                       "meters3 / meters2", "meters3 / meters2", "", 0);

    thickness.set_attrs("internal", "thickness (temporary storage)",
                        "meters", "meters", "", 0);
  }
}

GeometryEvolution::GeometryEvolution(IceGrid::ConstPtr grid)
  : Component(grid) {
  m_impl = new Impl(grid);
}

GeometryEvolution::~GeometryEvolution() {
  delete m_impl;
}

void GeometryEvolution::init(const InputOptions &opts) {
  this->init_impl(opts);
}

void GeometryEvolution::init_impl(const InputOptions &opts) {
  (void) opts;
  // empty: the default implementation has no state
}

void GeometryEvolution::reset() {
  m_impl->conservation_error.set(0.0);
}

const array::Scalar& GeometryEvolution::flux_divergence() const {
  return m_impl->flux_divergence;
}

const array::Staggered1& GeometryEvolution::flux_staggered() const {
  return m_impl->flux_staggered;
}

const array::Scalar& GeometryEvolution::top_surface_mass_balance() const {
  return m_impl->effective_SMB;
}

const array::Scalar& GeometryEvolution::bottom_surface_mass_balance() const {
  return m_impl->effective_BMB;
}

const array::Scalar& GeometryEvolution::thickness_change_due_to_flow() const {
  return m_impl->thickness_change;
}

const array::Scalar& GeometryEvolution::area_specific_volume_change_due_to_flow() const {
  return m_impl->ice_area_specific_volume_change;
}

const array::Scalar& GeometryEvolution::conservation_error() const {
  return m_impl->conservation_error;
}

/*!
 * @param[in] geometry ice geometry
 * @param[in] dt time step, seconds
 * @param[in] advective_velocity advective (SSA) velocity
 * @param[in] diffusive_flux diffusive (SIA) flux
 * @param[in] velocity_bc_values advective velocity Dirichlet B.C. values
 * @param[in] thickness_bc_mask ice thickness Dirichlet B.C. mask
 *
 * Results are stored in internal fields accessible using getters.
 */
void GeometryEvolution::flow_step(const Geometry &geometry, double dt,
                                  const array::Vector    &advective_velocity,
                                  const array::Staggered &diffusive_flux,
                                  const array::Scalar  &thickness_bc_mask) {

  m_impl->profile.begin("ge.update_ghosted_copies");
  {
    // make ghosted copies of input fields
    m_impl->ice_thickness.copy_from(geometry.ice_thickness);
    m_impl->area_specific_volume.copy_from(geometry.ice_area_specific_volume);
    m_impl->sea_level.copy_from(geometry.sea_level_elevation);
    m_impl->bed_elevation.copy_from(geometry.bed_elevation);
    m_impl->input_velocity.copy_from(advective_velocity);

    // Compute cell_type and surface_elevation. Ghosts of results are updated.
    m_impl->gc.compute(m_impl->sea_level,          // in (uses ghosts)
                       m_impl->bed_elevation,      // in (uses ghosts)
                       m_impl->ice_thickness,      // in (uses ghosts)
                       m_impl->cell_type,          // out (ghosts are updated)
                       m_impl->surface_elevation); // out (ghosts are updated)
  }
  m_impl->profile.end("ge.update_ghosted_copies");

  // Derived classes can include modifications for regional runs.
  m_impl->profile.begin("ge.interface_fluxes");
  compute_interface_fluxes(m_impl->cell_type,          // in (uses ghosts)
                           m_impl->ice_thickness,      // in (uses ghosts)
                           m_impl->input_velocity,     // in (uses ghosts)
                           diffusive_flux,             // in
                           m_impl->flux_staggered);    // out
  m_impl->profile.end("ge.interface_fluxes");

  m_impl->flux_staggered.update_ghosts();

  {
    // allocate temporary storage (FIXME: at some point I should evaluate whether it's OK
    // to allocate this every time step)
    array::Staggered flux_limited(m_grid, "limited_ice_flux");

    make_nonnegative_preserving(dt,
                                m_impl->ice_thickness, // in (uses ghosts)
                                m_impl->flux_staggered, // in (uses ghosts)
                                flux_limited);

    m_impl->flux_staggered.copy_from(flux_limited);
  }

  m_impl->profile.begin("ge.flux_divergence");
  compute_flux_divergence(dt,                         // in
                          m_impl->flux_staggered,     // in (uses ghosts)
                          thickness_bc_mask,          // in
                          m_impl->conservation_error, // in/out
                          m_impl->flux_divergence);   // out
  m_impl->profile.end("ge.flux_divergence");

  // This is where part_grid is implemented.
  m_impl->profile.begin("ge.update_in_place");
  update_in_place(dt,                            // in
                  m_impl->bed_elevation,         // in
                  m_impl->sea_level,             // in
                  m_impl->flux_divergence,       // in
                  m_impl->ice_thickness,         // in/out
                  m_impl->area_specific_volume); // in/out
  m_impl->profile.end("ge.update_in_place");

  // Compute ice thickness and area specific volume changes.
  m_impl->profile.begin("ge.compute_changes");
  {
    m_impl->ice_thickness.add(-1.0, geometry.ice_thickness,
                              m_impl->thickness_change);
    m_impl->area_specific_volume.add(-1.0, geometry.ice_area_specific_volume,
                                     m_impl->ice_area_specific_volume_change);
  }
  m_impl->profile.end("ge.compute_changes");

  // Computes the numerical conservation error and corrects ice_thickness_change and
  // ice_area_specific_volume_change. We can do this here because
  // compute_surface_and_basal_mass_balance() preserves non-negativity.
  //
  // Note that here we use the "old" ice geometry.
  //
  // This computation is purely local.
  m_impl->profile.begin("ge.ensure_nonnegativity");
  ensure_nonnegativity(geometry.ice_thickness,                  // in
                       geometry.ice_area_specific_volume,       // in
                       m_impl->thickness_change,                // in/out
                       m_impl->ice_area_specific_volume_change, // in/out
                       m_impl->conservation_error);             // in/out
  m_impl->profile.end("ge.ensure_nonnegativity");

  // Now the caller can compute
  //
  // H_new    = H_old + thickness_change
  // Href_new = Href_old + ice_area_specific_volume_change.

  // calving is a separate issue
}

void GeometryEvolution::source_term_step(const Geometry &geometry, double dt,
                                         const array::Scalar  &thickness_bc_mask,
                                         const array::Scalar    &surface_mass_balance_rate,
                                         const array::Scalar    &basal_melt_rate) {

  m_impl->profile.begin("ge.source_terms");
  compute_surface_and_basal_mass_balance(dt,                        // in
                                         thickness_bc_mask,         // in
                                         geometry.ice_thickness,    // in
                                         geometry.cell_type,        // in
                                         surface_mass_balance_rate, // in
                                         basal_melt_rate,           // in
                                         m_impl->effective_SMB,     // out
                                         m_impl->effective_BMB);    // out
  m_impl->profile.end("ge.source_terms");

}

/*!
 * Apply changes due to flow to ice geometry and ice area specific volume.
 */
void GeometryEvolution::apply_flux_divergence(Geometry &geometry) const {
  geometry.ice_thickness.add(1.0, m_impl->thickness_change);
  geometry.ice_area_specific_volume.add(1.0, m_impl->ice_area_specific_volume_change);
}

/*!
 * Update geometry by applying changes due to surface and basal mass fluxes.
 *
 * Note: This method performs these changes in the same order as the code ensuring
 * non-negativity. This is important.
 */
void GeometryEvolution::apply_mass_fluxes(Geometry &geometry) const {

  const array::Scalar
    &dH_SMB  = top_surface_mass_balance(),
    &dH_BMB  = bottom_surface_mass_balance();
  array::Scalar &H = geometry.ice_thickness;

  array::AccessScope list{&H, &dH_SMB, &dH_BMB};
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // To preserve non-negativity of thickness we need to apply changes in this exact order.
      // (Recall that floating-point arithmetic is not associative.)
      const double H_new = (H(i, j) + dH_SMB(i, j)) + dH_BMB(i, j);

#if (Pism_DEBUG==1)
      if (H_new < 0.0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "H = %f (negative) at i=%d, j=%d",
                                      H_new, i, j);
      }
#endif

      H(i, j) = H_new;
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}


/*!
 * Prevent advective ice flow from floating ice to ice-free land, as well as in the
 * ice-free areas.
 *
 * Note: positive `input` corresponds to the flux from `current` to `neighbor`.
 */
static double limit_advective_flux(int current, int neighbor, double input) {

  // No flow from ice-free ocean:
  if ((ice_free_ocean(current) and input > 0.0) or
      (ice_free_ocean(neighbor) and input < 0.0)) {
    return 0.0;
  }

  // No flow from ice-free land:
  if ((ice_free_land(current) and input > 0.0) or
      (ice_free_land(neighbor) and input < 0.0)) {
    return 0.0;
  }

  // Case 1: Flow between grounded_ice and grounded_ice.
  if (grounded_ice(current) and grounded_ice(neighbor)) {
    return input;
  }

  // Cases 2 and 3: Flow between grounded_ice and floating_ice.
  if ((grounded_ice(current) and floating_ice(neighbor)) or
      (floating_ice(current) and grounded_ice(neighbor))) {
    return input;
  }

  // Cases 4 and 5: Flow between grounded_ice and ice_free_land.
  if ((grounded_ice(current) and ice_free_land(neighbor)) or
      (ice_free_land(current) and grounded_ice(neighbor))) {
    return input;
  }

  // Cases 6 and 7: Flow between grounded_ice and ice_free_ocean.
  if ((grounded_ice(current) and ice_free_ocean(neighbor)) or
      (ice_free_ocean(current) and grounded_ice(neighbor))) {
    return input;
  }

  // Case 8: Flow between floating_ice and floating_ice.
  if (floating_ice(current) and floating_ice(neighbor)) {
    return input;
  }

  // Cases 9 and 10: Flow between floating_ice and ice_free_land.
  if ((floating_ice(current) and ice_free_land(neighbor)) or
      (ice_free_land(current) and floating_ice(neighbor))) {
    // Disable all flow. This ensures that an ice shelf does not climb up a cliff.
    return 0.0;
  }

  // Cases 11 and 12: Flow between floating_ice and ice_free_ocean.
  if ((floating_ice(current) and ice_free_ocean(neighbor)) or
      (ice_free_ocean(current) and floating_ice(neighbor))) {
    return input;
  }

  // Case 13: Flow between ice_free_land and ice_free_land.
  if (ice_free_land(current) and ice_free_land(neighbor)) {
    return 0.0;
  }

  // Cases 14 and 15: Flow between ice_free_land and ice_free_ocean.
  if ((ice_free_land(current) and ice_free_ocean(neighbor)) or
      (ice_free_ocean(current) and ice_free_land(neighbor))) {
    return 0.0;
  }

  // Case 16: Flow between ice_free_ocean and ice_free_ocean.
  if (ice_free_ocean(current) and ice_free_ocean(neighbor)) {
    return 0.0;
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "cannot handle the case current=%d, neighbor=%d",
                                current, neighbor);
}

/*!
 * Prevent SIA-driven flow in ice shelves and ice-free areas.
 *
 * Note: positive `flux` corresponds to the flux from `current` to `neighbor`.
 */
static double limit_diffusive_flux(int current, int neighbor, double flux) {

  // No flow from ice-free ocean:
  if ((ice_free_ocean(current) and flux > 0.0) or
      (ice_free_ocean(neighbor) and flux < 0.0)) {
    return 0.0;
  }

  // No flow from ice-free land:
  if ((ice_free_land(current) and flux > 0.0) or
      (ice_free_land(neighbor) and flux < 0.0)) {
    return 0.0;
  }

  // Case 1: Flow between grounded_ice and grounded_ice.
  if (grounded_ice(current) and grounded_ice(neighbor)) {
    return flux;
  }

  // Cases 2 and 3: Flow between grounded_ice and floating_ice.
  if ((grounded_ice(current) and floating_ice(neighbor)) or
      (floating_ice(current) and grounded_ice(neighbor))) {
    return flux;
  }

  // Cases 4 and 5: Flow between grounded_ice and ice_free_land.
  if ((grounded_ice(current) and ice_free_land(neighbor)) or
      (ice_free_land(current) and grounded_ice(neighbor))) {
    return flux;
  }

  // Cases 6 and 7: Flow between grounded_ice and ice_free_ocean.
  if ((grounded_ice(current) and ice_free_ocean(neighbor)) or
      (ice_free_ocean(current) and grounded_ice(neighbor))) {
    return flux;
  }

  // Case 8: Flow between floating_ice and floating_ice.
  if (floating_ice(current) and floating_ice(neighbor)) {
    // no diffusive flux in ice shelves
    return 0.0;
  }

  // Cases 9 and 10: Flow between floating_ice and ice_free_land.
  if ((floating_ice(current) and ice_free_land(neighbor)) or
      (ice_free_land(current) and floating_ice(neighbor))) {
    // Disable all flow. This ensures that an ice shelf does not climb up a cliff.
    return 0.0;
  }

  // Cases 11 and 12: Flow between floating_ice and ice_free_ocean.
  if ((floating_ice(current) and ice_free_ocean(neighbor)) or
      (ice_free_ocean(current) and floating_ice(neighbor))) {
    return 0.0;
  }

  // Case 13: Flow between ice_free_land and ice_free_land.
  if (ice_free_land(current) and ice_free_land(neighbor)) {
    return 0.0;
  }

  // Cases 14 and 15: Flow between ice_free_land and ice_free_ocean.
  if ((ice_free_land(current) and ice_free_ocean(neighbor)) or
      (ice_free_ocean(current) and ice_free_land(neighbor))) {
    return 0.0;
  }

  // Case 16: Flow between ice_free_ocean and ice_free_ocean.
  if (ice_free_ocean(current) and ice_free_ocean(neighbor)) {
    return 0.0;
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "cannot handle the case current=%d, neighbor=%d",
                                current, neighbor);
}

/*!
 * Combine advective velocity and the diffusive flux on the staggered grid with the ice thickness to
 * compute the total flux through cell interfaces.
 *
 * Uses first-order upwinding to compute the advective flux.
 *
 * Limits the diffusive flux to prevent SIA-driven flow in the ocean and ice-free areas.
 */
void GeometryEvolution::compute_interface_fluxes(const array::CellType1 &cell_type,
                                                 const array::Scalar        &ice_thickness,
                                                 const array::Vector        &velocity,
                                                 const array::Staggered     &diffusive_flux,
                                                 array::Staggered           &output) {

  array::AccessScope list{&cell_type, &velocity, &ice_thickness,
                               &diffusive_flux, &output};

  ParallelSection loop(m_grid->com);
  try {
    // compute advective fluxes and put them in output
    for (Points p(*m_grid); p; p.next()) {
      const int
        i  = p.i(),
        j  = p.j(),
        M  = cell_type(i, j);

      const double H = ice_thickness(i, j);
      const Vector2d V  = velocity(i, j);

      for (int n = 0; n < 2; ++n) {
        const int
          oi  = 1 - n,               // offset in the i direction
          oj  = n,                   // offset in the j direction
          i_n = i + oi,              // i index of a neighbor
          j_n = j + oj;              // j index of a neighbor

        const int M_n = cell_type(i_n, j_n);

        // advective velocity at the current interface
        double v = 0.0;
        {
          Vector2d V_n  = velocity(i_n, j_n);
          int
            W   = icy(M),
            W_n = icy(M_n);

          auto v_staggered = (W * V + W_n * V_n) / std::max(W + W_n, 1);
          v = n == 0 ? v_staggered.u : v_staggered.v;
        }

        // advective flux
        const double
          H_n         = ice_thickness(i_n, j_n),
          Q_advective = v * (v > 0.0 ? H : H_n); // first order upwinding

        output(i, j, n) = Q_advective;
      } // end of the loop over neighbors (n)
    }

    // limit the advective flux and add the diffusive flux to it to get the total
    for (Points p(*m_grid); p; p.next()) {
      const int
        i = p.i(),
        j = p.j(),
        M = cell_type(i, j);

      for (int n = 0; n < 2; ++n) {
        const int
          oi  = 1 - n,               // offset in the i direction
          oj  = n,                   // offset in the j direction
          i_n = i + oi,              // i index of a neighbor
          j_n = j + oj;              // j index of a neighbor

        const int M_n = cell_type(i_n, j_n);

        // diffusive flux
        const double
          Q_diffusive = limit_diffusive_flux(M, M_n, diffusive_flux(i, j, n)),
          Q_advective = limit_advective_flux(M, M_n, output(i, j, n));

        output(i, j, n) = Q_diffusive + Q_advective;
      } // end of the loop over n
    }

  } catch (...) {
    loop.failed();
  }
  loop.check();
}

/*!
 * Compute flux divergence using cell interface fluxes on the staggered grid.
 *
 * The flux divergence at *ice thickness* Dirichlet B.C. locations is set to zero.
 */
void GeometryEvolution::compute_flux_divergence(double dt,
                                                const array::Staggered1 &flux,
                                                const array::Scalar &thickness_bc_mask,
                                                array::Scalar &conservation_error,
                                                array::Scalar &output) {
  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  array::AccessScope list{&flux, &thickness_bc_mask, &conservation_error, &output};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      auto Q = flux.star(i, j);

      double divQ = (Q.e - Q.w) / dx + (Q.n - Q.s) / dy;

      if (thickness_bc_mask(i, j) > 0.5) {
        // the thickness change would have been equal to -divQ*dt. By keeping ice
        // thickness fixed we *add* divQ*dt meters of ice.
        conservation_error(i, j) += divQ * dt; // units: meters
        output(i, j)              = 0.0;
      } else {
        output(i, j) = divQ;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

/*!
 * Update ice thickness and area_specific_volume *in place*.
 *
 * It would be better to compute the change in ice thickness and area_specific_volume and then apply
 * them, but it would require re-writing all the part-grid code from scratch. So, I make copies of
 * ice thickness and area_specific_volume, use this old code, then compute differences to get changes.
 * Compute ice thickness changes due to the flow of the ice.
 *
 * @param[in] dt time step, seconds
 * @param[in] bed_elevation bed elevation, meters
 * @param[in] sea_level sea level elevation
 * @param[in] ice_thickness ice thickness
 * @param[in] area_specific_volume area-specific volume (m3/m2)
 * @param[in] flux_divergence flux divergence
 * @param[out] thickness_change ice thickness change due to flow
 * @param[out] area_specific_volume_change area specific volume change due to flow
 */
void GeometryEvolution::update_in_place(double dt,
                                        const array::Scalar &bed_topography,
                                        const array::Scalar &sea_level,
                                        const array::Scalar &flux_divergence,
                                        array::Scalar &ice_thickness,
                                        array::Scalar &area_specific_volume) {

  m_impl->gc.compute(sea_level, bed_topography, ice_thickness,
                     m_impl->cell_type, m_impl->surface_elevation);

  array::AccessScope list{&ice_thickness, &flux_divergence};

  if (m_impl->use_part_grid) {
    m_impl->residual.set(0.0);

    // Store ice thickness. We need this copy to make sure that modifying ice_thickness in the loop
    // below does not affect the computation of the threshold thickness. (Note that
    // part_grid_threshold_thickness uses neighboring values of the mask, ice thickness, and surface
    // elevation.)
    m_impl->thickness.copy_from(ice_thickness);

    list.add({&area_specific_volume, &m_impl->residual, &m_impl->thickness,
          &m_impl->surface_elevation, &bed_topography, &m_impl->cell_type});
  }

#if (Pism_DEBUG==1)
  const double Lz = m_grid->Lz();
#endif

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double divQ = flux_divergence(i, j);

      if (m_impl->use_part_grid) {
        if (m_impl->cell_type.ice_free_ocean(i, j) and m_impl->cell_type.next_to_ice(i, j)) {
          assert(divQ <= 0.0);
          // Add the flow contribution to this partially filled cell.
          area_specific_volume(i, j) += -divQ * dt;

          double threshold = part_grid_threshold_thickness(m_impl->cell_type.star_int(i, j),
                                                           m_impl->thickness.star(i, j),
                                                           m_impl->surface_elevation.star(i, j),
                                                           bed_topography(i, j));

          // if threshold is zero, turn all the area specific volume into ice thickness, with zero
          // residual
          if (threshold == 0.0) {
            threshold = area_specific_volume(i, j);
          }

          if (area_specific_volume(i, j) >= threshold) {
            ice_thickness(i, j)        += threshold;
            m_impl->residual(i, j)      = area_specific_volume(i, j) - threshold;
            area_specific_volume(i, j)  = 0.0;
          }

          // In this case the flux goes into the area_specific_volume variable and does not directly
          // contribute to ice thickness at this location.
          divQ = 0.0;
        }
      } // end of if (use_part_grid)

      ice_thickness(i, j) += - dt * divQ;

#if (Pism_DEBUG==1)
      if (ice_thickness(i, j) > Lz) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "ice thickness exceeds Lz at i=%d, j=%d (H=%f, Lz=%f)",
                                      i, j, ice_thickness(i, j), Lz);
      }
#endif
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  ice_thickness.update_ghosts();

  // Compute the mask corresponding to the new thickness.
  m_impl->gc.compute_mask(sea_level, bed_topography, ice_thickness, m_impl->cell_type);

  /*
    Redistribute residual ice mass from subgrid-scale parameterization.

    See [@ref Albrechtetal2011].
  */
  if (m_impl->use_part_grid) {
    const int max_n_iterations = m_config->get_number("geometry.part_grid.max_iterations");

    bool done = false;
    for (int i = 0; i < max_n_iterations and not done; ++i) {
      m_log->message(4, "redistribution iteration %d\n", i);

      // this call may set done to true
      residual_redistribution_iteration(bed_topography,
                                        sea_level,
                                        m_impl->surface_elevation,
                                        ice_thickness,
                                        m_impl->cell_type,
                                        area_specific_volume,
                                        m_impl->residual,
                                        done);
    }

    if (not done) {
      m_log->message(2,
                     "WARNING: not done redistributing mass after %d iterations, remaining residual: %f m^3.\n",
                     max_n_iterations, sum(m_impl->residual) * m_grid->cell_area());

      // Add residual to ice thickness, preserving total ice mass. (This is not great, but
      // better than losing mass.)
      ice_thickness.add(1.0, m_impl->residual);
      m_impl->residual.set(0.0);
    }
  }
}

//! @brief Perform one iteration of the residual mass redistribution.
/*!
  @param[in] bed_topography bed elevation
  @param[in] sea_level sea level elevation
  @param[in,out] ice_surface_elevation surface elevation; used as temp. storage
  @param[in,out] ice_thickness ice thickness; updated
  @param[in,out] cell_type cell type mask; used as temp. storage
  @param[in,out] area_specific_volume area specific volume; updated
  @param[in,out] residual ice volume that still needs to be distributed; updated
  @param[in,out] done result flag: true if this iteration should be the last one
 */
void GeometryEvolution::residual_redistribution_iteration(const array::Scalar  &bed_topography,
                                                          const array::Scalar  &sea_level,
                                                          array::Scalar1       &ice_surface_elevation,
                                                          array::Scalar        &ice_thickness,
                                                          array::CellType1 &cell_type,
                                                          array::Scalar        &area_specific_volume,
                                                          array::Scalar        &residual,
                                                          bool &done) {

  m_impl->gc.compute_mask(sea_level, bed_topography, ice_thickness, cell_type);

  const Direction directions[4] = {North, East, South, West};

  // First step: distribute residual mass
  {
    // will be destroyed at the end of the block
    array::AccessScope list{&cell_type, &ice_thickness, &area_specific_volume, &residual};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (residual(i, j) <= 0.0) {
        continue;
      }

      auto m = cell_type.star(i, j);

      int N = 0; // number of empty or partially filled neighbors
      for (unsigned int n = 0; n < 4; ++n) {
        const Direction direction = directions[n];
        if (ice_free_ocean(m[direction])) {
          N++;
        }
      }

      if (N > 0)  {
        // Remaining ice mass will be redistributed equally among all adjacent
        // ice-free-ocean cells (is there a more physical way?)
        residual(i, j) /= N;
      } else {
        // Conserve mass, but (possibly) create a "ridge" at the shelf
        // front
        ice_thickness(i, j) += residual(i, j);
        residual(i, j) = 0.0;
      }
    }

    residual.update_ghosts();

    // update area_specific_volume using adjusted residuals
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (cell_type.ice_free_ocean(i, j)) {
        area_specific_volume(i, j) += (residual(i + 1, j) +
                                       residual(i - 1, j) +
                                       residual(i, j + 1) +
                                       residual(i, j - 1));
      }

    }

    residual.set(0.0);
  }

  ice_thickness.update_ghosts();

  // Store ice thickness. We need this copy to make sure that modifying ice_thickness in the loop
  // below does not affect the computation of the threshold thickness. (Note that
  // part_grid_threshold_thickness uses neighboring values of the mask, ice thickness, and surface
  // elevation.)
  m_impl->thickness.copy_from(ice_thickness);

  // The loop above updated ice_thickness, so we need to re-calculate the mask and the
  // surface elevation:
  m_impl->gc.compute(sea_level, bed_topography, ice_thickness, cell_type, ice_surface_elevation);

  double remaining_residual = 0.0;

  // Second step: we need to redistribute residual ice volume if
  // neighbors which gained redistributed ice also become full.
  {
    // will be destroyed at the end of the block
    array::AccessScope list{&m_impl->thickness, &ice_thickness,
        &ice_surface_elevation, &bed_topography, &cell_type};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (area_specific_volume(i, j) <= 0.0) {
        continue;
      }

      double threshold = part_grid_threshold_thickness(cell_type.star_int(i, j),
                                                       m_impl->thickness.star(i, j),
                                                       ice_surface_elevation.star(i, j),
                                                       bed_topography(i, j));

      // if threshold is zero, turn all the area specific volume into ice thickness, with zero
      // residual
      if (threshold == 0.0) {
        threshold = area_specific_volume(i, j);
      }

      if (area_specific_volume(i, j) >= threshold) {
        ice_thickness(i, j)        += threshold;
        residual(i, j)              = area_specific_volume(i, j) - threshold;
        area_specific_volume(i, j)  = 0.0;

        remaining_residual += residual(i, j);
      }
    }
  }

  // check if redistribution should be run once more
  remaining_residual = GlobalSum(m_grid->com, remaining_residual);

  if (remaining_residual > 0.0) {
    done = false;
  } else {
    done = true;
  }

  ice_thickness.update_ghosts();
}

/*!
 * Correct `thickness_change` and `area_specific_volume_change` so that applying them will not
 * result in negative `ice_thickness` and `area_specific_volume`.
 *
 * Compute the the amount of ice that is added to preserve non-negativity of ice thickness.
 *
 * @param[in] ice_thickness ice thickness (m)
 * @param[in] area_specific_volume area-specific volume (m3/m2)
 * @param[in,out] thickness_change "proposed" thickness change (m)
 * @param[in,out] area_specific_volume_change "proposed" area-specific volume change (m3/m2)
 * @param[in,out] conservation_error computed conservation error (m)
 *
 * This computation is purely local.
 */
void GeometryEvolution::ensure_nonnegativity(const array::Scalar &ice_thickness,
                                             const array::Scalar &area_specific_volume,
                                             array::Scalar &thickness_change,
                                             array::Scalar &area_specific_volume_change,
                                             array::Scalar &conservation_error) {

  array::AccessScope list{&ice_thickness, &area_specific_volume, &thickness_change,
      &area_specific_volume_change, &conservation_error};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double
        H  = ice_thickness(i, j),
        dH = thickness_change(i, j);

      // applying thickness_change will lead to negative thickness
      if (H + dH < 0.0) {
        thickness_change(i, j)    = -H;
        conservation_error(i, j) += - (H + dH);
      }

      const double
        V  = area_specific_volume(i, j),
        dV = area_specific_volume_change(i, j);

      if (V + dV < 0.0) {
        area_specific_volume_change(i, j)  = -V;
        conservation_error(i, j)          += - (V + dV);
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

/*!
 * Given ice thickness `H` and the "proposed" change `dH`, compute the corrected change preserving
 * non-negativity.
 */
static inline double effective_change(double H, double dH) {
  if (H + dH <= 0) {
    return -H;
  } else {
    return dH;
  }
}

/*!
 * Compute effective surface and basal mass balance.
 *
 * @param[in] dt time step, seconds
 * @param[in] thickness_bc_mask mask specifying ice thickness Dirichlet B.C. locations
 * @param[in] ice_thickness ice thickness, m
 * @param[in] thickness_change thickness change due to flow, m
 * @param[in] cell_type cell type mask
 * @param[in] smb_rate top surface mass balance rate, kg m-2 s-1
 * @param[in] basal_melt_rate basal melt rate, m s-1
 * @param[out] effective_smb effective top surface mass balance, m
 * @param[out] effective_bmb effective basal mass balance, m
 *
 * This computation is purely local.
 */
void GeometryEvolution::compute_surface_and_basal_mass_balance(double dt,
                                                               const array::Scalar      &thickness_bc_mask,
                                                               const array::Scalar        &ice_thickness,
                                                               const array::CellType &cell_type,
                                                               const array::Scalar        &smb_flux,
                                                               const array::Scalar        &basal_melt_rate,
                                                               array::Scalar              &effective_SMB,
                                                               array::Scalar              &effective_BMB) {

  array::AccessScope list{&ice_thickness,
      &smb_flux, &basal_melt_rate, &cell_type, &thickness_bc_mask,
      &effective_SMB, &effective_BMB};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // Don't modify ice thickness at Dirichlet B.C. locations and in the ice-free ocean.
      if (thickness_bc_mask.as_int(i, j) == 1 or cell_type.ice_free_ocean(i, j)) {
        effective_SMB(i, j) = 0.0;
        effective_BMB(i, j) = 0.0;
        continue;
      }

      const double H = ice_thickness(i, j);

      // Thickness change due to the surface mass balance
      //
      // Note that here we convert surface mass balance from [kg m-2 s-1] to [m s-1].
      double dH_SMB = effective_change(H, dt * smb_flux(i, j) / m_impl->ice_density);

      // Thickness change due to the basal mass balance
      //
      // Note that basal_melt_rate is in [m s-1]. Here the negative sign converts the melt rate into
      // mass balance.
      double dH_BMB = effective_change(H + dH_SMB,
                                       dt * (m_impl->use_bmr ? -basal_melt_rate(i, j) : 0.0));

      effective_SMB(i, j) = dH_SMB;
      effective_BMB(i, j) = dH_BMB;
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

namespace diagnostics {

/*! @brief Report the divergence of the ice flux. */
class FluxDivergence : public Diag<GeometryEvolution>
{
public:
  FluxDivergence(const GeometryEvolution *m)
    : Diag<GeometryEvolution>(m) {
    m_vars = {model->flux_divergence().metadata()};
  }
protected:
  array::Array::Ptr compute_impl() const {
    array::Scalar::Ptr result(new array::Scalar(m_grid, "flux_divergence"));
    result->metadata(0) = m_vars[0];

    result->copy_from(model->flux_divergence());

    return result;
  }
};

/*! @brief Report mass flux on the staggered grid. */
class FluxStaggered : public Diag<GeometryEvolution>
{
public:
  FluxStaggered(const GeometryEvolution *m)
    : Diag<GeometryEvolution>(m) {
    m_vars = {model->flux_staggered().metadata(0),
      model->flux_staggered().metadata(1)};
  }
protected:
  array::Array::Ptr compute_impl() const {
    auto result = std::make_shared<array::Staggered>(m_grid, "flux_staggered");

    result->metadata(0) = m_vars[0];
    result->metadata(1) = m_vars[1];

    const array::Staggered &input = model->flux_staggered();
    array::Staggered &output = *result.get();

    // FIXME: implement array::Staggered::copy_from()

    array::AccessScope list{&input, &output};

    ParallelSection loop(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        output(i, j, 0) = input(i, j, 0);
        output(i, j, 1) = input(i, j, 1);
      }
    } catch (...) {
      loop.failed();
    }
    loop.check();

    return result;
  }
};

} // end of namespace diagnostics

DiagnosticList GeometryEvolution::diagnostics_impl() const {
  using namespace diagnostics;
  typedef Diagnostic::Ptr Ptr;

  std::map<std::string, Ptr> result;
  result = {
    {"flux_staggered",               Ptr(new FluxStaggered(this))},
    {"flux_divergence",              Ptr(new FluxDivergence(this))},
  };
  return result;
}

RegionalGeometryEvolution::RegionalGeometryEvolution(IceGrid::ConstPtr grid)
  : GeometryEvolution(grid),
    m_no_model_mask(m_grid, "no_model_mask") {

  m_no_model_mask.set_interpolation_type(NEAREST);

  m_no_model_mask.set_attrs("model_mask", "'no model' mask", "", "", "", 0);
}

void RegionalGeometryEvolution::set_no_model_mask_impl(const array::Scalar &mask) {
  m_no_model_mask.copy_from(mask);
}

/*!
 * Disable ice flow in "no model" areas.
 */
void RegionalGeometryEvolution::compute_interface_fluxes(const array::CellType1 &cell_type,
                                                         const array::Scalar        &ice_thickness,
                                                         const array::Vector        &velocity,
                                                         const array::Staggered     &diffusive_flux,
                                                         array::Staggered           &output) {

  GeometryEvolution::compute_interface_fluxes(cell_type, ice_thickness,
                                              velocity, diffusive_flux,
                                              output);

  array::AccessScope list{&m_no_model_mask, &output};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const int M = m_no_model_mask.as_int(i, j);

      for (int n : {0, 1}) {
        const int
          oi  = 1 - n,               // offset in the i direction
          oj  = n,                   // offset in the j direction
          i_n = i + oi,              // i index of a neighbor
          j_n = j + oj;              // j index of a neighbor

        const int M_n = m_no_model_mask.as_int(i_n, j_n);

        if (not (M == 0 and M_n == 0)) {
          output(i, j, n) = 0.0;
        }
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

/*!
 * Set surface and basal mass balance to zero in "no model" areas.
 */
void RegionalGeometryEvolution::compute_surface_and_basal_mass_balance(double dt,
                                                                       const array::Scalar      &thickness_bc_mask,
                                                                       const array::Scalar        &ice_thickness,
                                                                       const array::CellType &cell_type,
                                                                       const array::Scalar        &surface_mass_flux,
                                                                       const array::Scalar        &basal_melt_rate,
                                                                       array::Scalar              &effective_SMB,
                                                                       array::Scalar              &effective_BMB) {
  GeometryEvolution::compute_surface_and_basal_mass_balance(dt,
                                                            thickness_bc_mask,
                                                            ice_thickness,
                                                            cell_type,
                                                            surface_mass_flux,
                                                            basal_melt_rate,
                                                            effective_SMB,
                                                            effective_BMB);

  array::AccessScope list{&m_no_model_mask, &effective_SMB, &effective_BMB};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_no_model_mask(i, j) > 0.5) {
        effective_SMB(i, j) = 0.0;
        effective_BMB(i, j) = 0.0;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

void GeometryEvolution::set_no_model_mask(const array::Scalar &mask) {
  this->set_no_model_mask_impl(mask);
}

void GeometryEvolution::set_no_model_mask_impl(const array::Scalar &mask) {
  (void) mask;
  // the default implementation is a no-op
}

void grounding_line_flux(const array::CellType1 &cell_type,
                         const array::Staggered1 &flux,
                         double dt,
                         bool add_values,
                         array::Scalar &output) {

  using mask::grounded;

  auto grid = output.grid();

  const double
    dx = grid->dx(),
    dy = grid->dy();

  auto cell_area = grid->cell_area();

  auto ice_density = grid->ctx()->config()->get_number("constants.ice.density");

  array::AccessScope list{&cell_type, &flux, &output};

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double result = 0.0;

      if (cell_type.ocean(i ,j)) {
        auto M = cell_type.star(i, j);
        auto Q = flux.star(i, j);

        if (grounded(M.n)) {
          result += Q.n * dx;
        }

        if (grounded(M.e)) {
          result += Q.e * dy;
        }

        if (grounded(M.s)) {
          result -= Q.s * dx;
        }

        if (grounded(M.w)) {
          result -= Q.w * dy;
        }

        // convert from "m^3 / s" to "kg / m^2"
        result *= dt * (ice_density / cell_area);
      }

      if (add_values) {
        output(i, j) += result;
      } else {
        output(i, j) = result;
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

/*!
 * Compute the total grounding line flux over a time step, in kg.
 */
double total_grounding_line_flux(const array::CellType1 &cell_type,
                                 const array::Staggered1 &flux,
                                 double dt) {
  using mask::grounded;

  auto grid = cell_type.grid();

  const double
    dx = grid->dx(),
    dy = grid->dy();

  auto ice_density = grid->ctx()->config()->get_number("constants.ice.density");

  double total_flux = 0.0;

  array::AccessScope list{&cell_type, &flux};

  ParallelSection loop(grid->com);
  try {
    for (Points p(*grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double volume_flux = 0.0;

      if (cell_type.ocean(i ,j)) {
        auto M = cell_type.star(i, j);
        auto Q = flux.star(i, j); // m^2 / s

        if (grounded(M.n)) {
          volume_flux += Q.n * dx;
        }

        if (grounded(M.e)) {
          volume_flux += Q.e * dy;
        }

        if (grounded(M.s)) {
          volume_flux -= Q.s * dx;
        }

        if (grounded(M.w)) {
          volume_flux -= Q.w * dy;
        }
      }

      // convert from "m^3 / s" to "kg" and sum up
      total_flux += volume_flux * dt * ice_density;
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  return GlobalSum(grid->com, total_flux);
}

} // end of namespace pism
