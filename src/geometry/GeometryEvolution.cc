/* Copyright (C) 2016, 2017, 2018, 2019 PISM Authors
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

#include "pism/util/iceModelVec.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Mask.hh"

#include "pism/geometry/part_grid_threshold_thickness.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Profiling.hh"

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
  IceModelVec2S flux_divergence;

  //! Conservation error due to enforcing non-negativity of ice thickness.
  IceModelVec2S conservation_error;

  //! Effective surface mass balance.
  IceModelVec2S effective_SMB;

  //! Effective basal mass balance.
  IceModelVec2S effective_BMB;

  //! Change in ice thickness due to flow during the last time step.
  IceModelVec2S thickness_change;

  //! Change in the ice area-specific volume due to flow during the last time step.
  IceModelVec2S ice_area_specific_volume_change;

  //! Flux through cell interfaces. Ghosted.
  IceModelVec2Stag flux_staggered;

  // Work space
  IceModelVec2V        input_velocity;       // ghosted copy; not modified
  IceModelVec2S        bed_elevation;        // ghosted copy; not modified
  IceModelVec2S        sea_level;            // ghosted copy; not modified
  IceModelVec2S        ice_thickness;        // ghosted; updated in place
  IceModelVec2S        area_specific_volume; // ghosted; updated in place
  IceModelVec2S        surface_elevation;    // ghosted; updated to maintain consistency
  IceModelVec2CellType cell_type;            // ghosted; updated to maintain consistency
  IceModelVec2S        residual;             // ghosted; temporary storage
  IceModelVec2S        thickness;            // ghosted; temporary storage
  IceModelVec2Int      velocity_bc_mask;
};

GeometryEvolution::Impl::Impl(IceGrid::ConstPtr grid)
  : profile(grid->ctx()->profiling()),
    gc(*grid->ctx()->config()) {

  Config::ConstPtr config = grid->ctx()->config();

  gc.set_icefree_thickness(config->get_double("geometry.ice_free_thickness_standard"));

  // constants
  {
    ice_density   = config->get_double("constants.ice.density");
    use_bmr       = config->get_boolean("geometry.update.use_basal_melt_rate");
    use_part_grid = config->get_boolean("geometry.part_grid.enabled");
  }

  // reported quantities
  {
    // This is the only reported field that is ghosted (we need ghosts to compute flux divergence).
    flux_staggered.create(grid, "flux_staggered", WITH_GHOSTS);
    flux_staggered.set_attrs("diagnostic", "fluxes through cell interfaces (sides)"
                             " on the staggered grid",
                             "m2 s-1", "");
    flux_staggered.metadata().set_string("glaciological_units", "m2 year-1");

    flux_divergence.create(grid, "flux_divergence", WITHOUT_GHOSTS);
    flux_divergence.set_attrs("diagnostic", "flux divergence", "m s-1", "");
    flux_divergence.metadata().set_string("glaciological_units", "m year-1");

    conservation_error.create(grid, "conservation_error", WITHOUT_GHOSTS);
    conservation_error.set_attrs("diagnostic",
                                 "conservation error due to enforcing non-negativity of"
                                 " ice thickness (over the last time step)", "meters", "");

    effective_SMB.create(grid, "effective_SMB", WITHOUT_GHOSTS);
    effective_SMB.set_attrs("internal", "effective surface mass balance over the last time step",
                            "meters", "");

    effective_BMB.create(grid, "effective_BMB", WITHOUT_GHOSTS);
    effective_BMB.set_attrs("internal", "effective basal mass balance over the last time step",
                            "meters", "");

    thickness_change.create(grid, "thickness_change", WITHOUT_GHOSTS);
    thickness_change.set_attrs("internal", "change in thickness due to flow", "meters", "");

    ice_area_specific_volume_change.create(grid, "ice_area_specific_volume_change", WITHOUT_GHOSTS);
    ice_area_specific_volume_change.set_attrs("interval",
                                              "change in area-specific volume due to flow",
                                              "meters3 / meters2", "");
  }

  // internal storage
  {
    input_velocity.create(grid, "input_velocity", WITH_GHOSTS);
    input_velocity.set_attrs("internal", "ghosted copy of the input velocity",
                             "meters / second", "");

    bed_elevation.create(grid, "bed_elevation", WITH_GHOSTS);
    bed_elevation.set_attrs("internal", "ghosted copy of the bed elevation",
                            "meters", "");

    sea_level.create(grid, "sea_level", WITH_GHOSTS);
    sea_level.set_attrs("internal", "ghosted copy of the sea level elevation",
                        "meters", "");

    ice_thickness.create(grid, "ice_thickness", WITH_GHOSTS);
    ice_thickness.set_attrs("internal", "working (ghosted) copy of the ice thickness",
                            "meters", "");

    area_specific_volume.create(grid, "area_specific_volume", WITH_GHOSTS);
    area_specific_volume.set_attrs("internal", "working (ghosted) copy of the area specific volume",
                                   "meters3 / meters2", "");

    surface_elevation.create(grid, "surface_elevation", WITH_GHOSTS);
    surface_elevation.set_attrs("internal", "working (ghosted) copy of the surface elevation",
                                "meters", "");

    cell_type.create(grid, "cell_type", WITH_GHOSTS);
    cell_type.set_attrs("internal", "working (ghosted) copy of the cell type mask",
                        "", "");

    residual.create(grid, "residual", WITH_GHOSTS);
    residual.set_attrs("internal", "residual area specific volume",
                       "meters3 / meters2", "");

    thickness.create(grid, "thickness", WITH_GHOSTS);
    thickness.set_attrs("internal", "thickness (temporary storage)",
                        "meters", "");

    velocity_bc_mask.create(grid, "velocity_bc_mask", WITH_GHOSTS);
    velocity_bc_mask.set_attrs("internal", "ghosted copy of the velocity B.C. mask"
                               " (1 at velocity B.C. location, 0 elsewhere)",
                               "", "");
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

const IceModelVec2S& GeometryEvolution::flux_divergence() const {
  return m_impl->flux_divergence;
}

const IceModelVec2Stag& GeometryEvolution::flux_staggered() const {
  return m_impl->flux_staggered;
}

const IceModelVec2S& GeometryEvolution::top_surface_mass_balance() const {
  return m_impl->effective_SMB;
}

const IceModelVec2S& GeometryEvolution::bottom_surface_mass_balance() const {
  return m_impl->effective_BMB;
}

const IceModelVec2S& GeometryEvolution::thickness_change_due_to_flow() const {
  return m_impl->thickness_change;
}

const IceModelVec2S& GeometryEvolution::area_specific_volume_change_due_to_flow() const {
  return m_impl->ice_area_specific_volume_change;
}

const IceModelVec2S& GeometryEvolution::conservation_error() const {
  return m_impl->conservation_error;
}

/*!
 * @param[in] geometry ice geometry
 * @param[in] dt time step, seconds
 * @param[in] advective_velocity advective (SSA) velocity
 * @param[in] diffusive_flux diffusive (SIA) flux
 * @param[in] velocity_bc_mask advective velocity Dirichlet B.C. mask
 * @param[in] velocity_bc_values advective velocity Dirichlet B.C. values
 * @param[in] thickness_bc_mask ice thickness Dirichlet B.C. mask
 *
 * Results are stored in internal fields accessible using getters.
 */
void GeometryEvolution::flow_step(const Geometry &geometry, double dt,
                                  const IceModelVec2V    &advective_velocity,
                                  const IceModelVec2Stag &diffusive_flux,
                                  const IceModelVec2Int  &velocity_bc_mask,
                                  const IceModelVec2Int  &thickness_bc_mask) {

  m_impl->profile.begin("ge.update_ghosted_copies");
  {
    // make ghosted copies of input fields
    m_impl->ice_thickness.copy_from(geometry.ice_thickness);
    m_impl->area_specific_volume.copy_from(geometry.ice_area_specific_volume);
    m_impl->sea_level.copy_from(geometry.sea_level_elevation);
    m_impl->bed_elevation.copy_from(geometry.bed_elevation);
    m_impl->input_velocity.copy_from(advective_velocity);
    m_impl->velocity_bc_mask.copy_from(velocity_bc_mask);

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
                           m_impl->velocity_bc_mask,   // in (uses ghosts)
                           diffusive_flux,             // in
                           m_impl->flux_staggered);    // out
  m_impl->profile.end("ge.interface_fluxes");

  m_impl->flux_staggered.update_ghosts();

  m_impl->profile.begin("ge.flux_divergence");
  compute_flux_divergence(m_impl->flux_staggered,   // in (uses ghosts)
                          thickness_bc_mask,        // in
                          m_impl->flux_divergence); // out
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
                       m_impl->conservation_error);             // out
  m_impl->profile.end("ge.ensure_nonnegativity");

  // Now the caller can compute
  //
  // H_new    = H_old + thickness_change
  // Href_new = Href_old + ice_area_specific_volume_change.

  // calving is a separate issue
}

void GeometryEvolution::source_term_step(const Geometry &geometry, double dt,
                                         const IceModelVec2Int  &thickness_bc_mask,
                                         const IceModelVec2S    &surface_mass_balance_rate,
                                         const IceModelVec2S    &basal_melt_rate) {

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

  const IceModelVec2S
    &dH_SMB  = top_surface_mass_balance(),
    &dH_BMB  = bottom_surface_mass_balance();
  IceModelVec2S &H = geometry.ice_thickness;

  IceModelVec::AccessList list{&H, &dH_SMB, &dH_BMB};
  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // To preserve non-negativity of thickness we need to apply changes in this exact order.
      // (Recall that floating-point arithmetic is not associative.)
      const double H_new = (H(i, j) + dH_SMB(i, j)) + dH_BMB(i, j);

#if (PISM_DEBUG==1)
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
 * Prevent advective ice flow from floating ice to ice-free land, as well as in the ice-free areas.
 */
static double limit_advective_velocity(int current, int neighbor, double velocity) {

  // Case 1: Flow between grounded_ice and grounded_ice.
  if (grounded_ice(current) and grounded_ice(neighbor)) {
    return velocity;
  }

  // Cases 2 and 3: Flow between grounded_ice and floating_ice.
  if ((grounded_ice(current) and floating_ice(neighbor)) or
      (floating_ice(current) and grounded_ice(neighbor))) {
    return velocity;
  }

  // Cases 4 and 5: Flow between grounded_ice and ice_free_land.
  if ((grounded_ice(current) and ice_free_land(neighbor)) or
      (ice_free_land(current) and grounded_ice(neighbor))) {
    return velocity;
  }

  // Cases 6 and 7: Flow between grounded_ice and ice_free_ocean.
  if ((grounded_ice(current) and ice_free_ocean(neighbor)) or
      (ice_free_ocean(current) and grounded_ice(neighbor))) {
    return velocity;
  }

  // Case 8: Flow between floating_ice and floating_ice.
  if (floating_ice(current) and floating_ice(neighbor)) {
    return velocity;
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
    return velocity;
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
 */
static double limit_diffusive_flux(int current, int neighbor, double flux) {

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
void GeometryEvolution::compute_interface_fluxes(const IceModelVec2CellType &cell_type,
                                                 const IceModelVec2S        &ice_thickness,
                                                 const IceModelVec2V        &velocity,
                                                 const IceModelVec2Int      &velocity_bc_mask,
                                                 const IceModelVec2Stag     &diffusive_flux,
                                                 IceModelVec2Stag           &output) {

  IceModelVec::AccessList list{&cell_type, &velocity, &velocity_bc_mask, &ice_thickness,
      &diffusive_flux, &output};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int
        i  = p.i(),
        j  = p.j(),
        M  = cell_type(i, j),
        BC = velocity_bc_mask.as_int(i, j);

      const double H = ice_thickness(i, j);
      const Vector2 V  = velocity(i, j);

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
          const Vector2 V_n  = velocity(i_n, j_n);

          // Regular case
          {
            if (icy(M) and icy(M_n)) {
              // Case 1: both sides of the interface are icy
              v = (n == 0 ? 0.5 * (V.u + V_n.u) : 0.5 * (V.v + V_n.v));

            } else if (icy(M) and ice_free(M_n)) {
              // Case 2: icy cell next to an ice-free cell
              v = (n == 0 ? V.u : V.v);

            } else if (ice_free(M) and icy(M_n)) {
              // Case 3: ice-free cell next to icy cell
              v = (n == 0 ? V_n.u : V_n.v);

            } else if (ice_free(M) and ice_free(M_n)) {
              // Case 4: both sides of the interface are ice-free
              v = 0.0;

            }
          }

          // The Dirichlet B.C. case:
          {
            const int BC_n = velocity_bc_mask.as_int(i_n, j_n);

            if (BC == 1 and BC_n == 1) {
              // Case 1: both sides of the interface are B.C. locations: average from
              // the regular grid onto the staggered grid.
              v = (n == 0 ? 0.5 * (V.u + V_n.u) : 0.5 * (V.v + V_n.v));

            } else if (BC == 1 and BC_n == 0) {
              // Case 2: at a Dirichlet B.C. location next to a regular location
              v = (n == 0 ? V.u : V.v);

            } else if (BC == 0 and BC_n == 1) {

              // Case 3: at a regular location next to a Dirichlet B.C. location
              v = (n == 0 ? V_n.u : V_n.v);

            } else {
              // Case 4: elsewhere.
              // No Dirichlet B.C. adjustment here.
            }

          } // end of the Dirichlet B.C. case

          // finally, limit advective velocities
          v = limit_advective_velocity(M, M_n, v);
        }

        // advective flux
        const double
          H_n         = ice_thickness(i_n, j_n),
          Q_advective = v * (v > 0.0 ? H : H_n); // first order upwinding

        // diffusive flux
        const double
          Q_diffusive = limit_diffusive_flux(M, M_n, diffusive_flux(i, j, n));

        output(i, j, n) = Q_diffusive + Q_advective;
      } // end of the loop over neighbors (n)
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
void GeometryEvolution::compute_flux_divergence(const IceModelVec2Stag &flux,
                                                const IceModelVec2Int &thickness_bc_mask,
                                                IceModelVec2S &output) {
  const double
    dx = m_grid->dx(),
    dy = m_grid->dy();

  IceModelVec::AccessList list{&flux, &thickness_bc_mask, &output};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (thickness_bc_mask(i, j) > 0.5) {
        output(i, j) = 0.0;
      } else {
        StarStencil<double> Q = flux.star(i, j);

        output(i, j) = (Q.e - Q.w) / dx + (Q.n - Q.s) / dy;
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
                                        const IceModelVec2S &bed_topography,
                                        const IceModelVec2S &sea_level,
                                        const IceModelVec2S &flux_divergence,
                                        IceModelVec2S &ice_thickness,
                                        IceModelVec2S &area_specific_volume) {

  m_impl->gc.compute(sea_level, bed_topography, ice_thickness,
                     m_impl->cell_type, m_impl->surface_elevation);

  IceModelVec::AccessList list{&ice_thickness, &flux_divergence};

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

#if (PISM_DEBUG==1)
  const double Lz = m_grid->Lz();
#endif

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double divQ = flux_divergence(i, j);

      if (m_impl->use_part_grid) {
        if (m_impl->cell_type.ice_free_ocean(i, j) and m_impl->cell_type.next_to_ice(i, j)) {
          // Add the flow contribution to this partially filled cell.
          area_specific_volume(i, j) += -divQ * dt;

          double threshold = part_grid_threshold_thickness(m_impl->cell_type.int_star(i, j),
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

#if (PISM_DEBUG==1)
      if (ice_thickness(i, j) > Lz) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "ice thickness exceeds Lz at i=%d, j=%d (H=%f, Lz=%f)",
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
    const int max_n_iterations = m_config->get_double("geometry.part_grid.max_iterations");

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
                     max_n_iterations, m_impl->residual.sum() * m_grid->cell_area());

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
void GeometryEvolution::residual_redistribution_iteration(const IceModelVec2S  &bed_topography,
                                                          const IceModelVec2S  &sea_level,
                                                          IceModelVec2S        &ice_surface_elevation,
                                                          IceModelVec2S        &ice_thickness,
                                                          IceModelVec2CellType &cell_type,
                                                          IceModelVec2S        &area_specific_volume,
                                                          IceModelVec2S        &residual,
                                                          bool &done) {

  m_impl->gc.compute_mask(sea_level, bed_topography, ice_thickness, cell_type);

  const Direction directions[4] = {North, East, South, West};

  // First step: distribute residual mass
  {
    // will be destroyed at the end of the block
    IceModelVec::AccessList list{&cell_type, &ice_thickness, &area_specific_volume, &residual};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (residual(i, j) <= 0.0) {
        continue;
      }

      StarStencil<int> m = cell_type.int_star(i, j);

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
    IceModelVec::AccessList list{&m_impl->thickness, &ice_thickness,
        &ice_surface_elevation, &bed_topography, &cell_type};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (area_specific_volume(i, j) <= 0.0) {
        continue;
      }

      double threshold = part_grid_threshold_thickness(cell_type.int_star(i, j),
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
 * Compute the `conservation_error`, i.e. the amount of ice that is added to preserve
 * non-negativity.
 *
 * @param[in] ice_thickness ice thickness (m)
 * @param[in] area_specific_volume area-specific volume (m3/m2)
 * @param[in,out] thickness_change "proposed" thickness change (m)
 * @param[in,out] area_specific_volume_change "proposed" area-specific volume change (m3/m2)
 * @param[out] conservation_error computed conservation error (m)
 *
 * This computation is purely local.
 */
void GeometryEvolution::ensure_nonnegativity(const IceModelVec2S &ice_thickness,
                                             const IceModelVec2S &area_specific_volume,
                                             IceModelVec2S &thickness_change,
                                             IceModelVec2S &area_specific_volume_change,
                                             IceModelVec2S &conservation_error) {

  IceModelVec::AccessList list{&ice_thickness, &area_specific_volume, &thickness_change,
      &area_specific_volume_change, &conservation_error};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      conservation_error(i, j) = 0.0;

      const double
        H  = ice_thickness(i, j),
        dH = thickness_change(i, j);

      // applying thickness_change will lead to negative thickness
      if (H + dH < 0.0) {
        thickness_change(i, j)    = H;
        conservation_error(i, j) += - (H + dH);
      }

      const double
        V  = area_specific_volume(i, j),
        dV = area_specific_volume_change(i, j);

      if (V + dV < 0.0) {
        area_specific_volume_change(i, j)  = V;
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
                                                               const IceModelVec2Int      &thickness_bc_mask,
                                                               const IceModelVec2S        &ice_thickness,
                                                               const IceModelVec2CellType &cell_type,
                                                               const IceModelVec2S        &smb_flux,
                                                               const IceModelVec2S        &basal_melt_rate,
                                                               IceModelVec2S              &effective_SMB,
                                                               IceModelVec2S              &effective_BMB) {

  IceModelVec::AccessList list{&ice_thickness,
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
  IceModelVec::Ptr compute_impl() const {
    IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "flux_divergence", WITHOUT_GHOSTS));
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
    m_vars = {model->flux_staggered().metadata()};
  }
protected:
  IceModelVec::Ptr compute_impl() const {
    IceModelVec2Stag::Ptr result(new IceModelVec2Stag(m_grid, "flux_staggered", WITHOUT_GHOSTS));
    result->metadata(0) = m_vars[0];

    const IceModelVec2Stag &input = model->flux_staggered();
    IceModelVec2Stag &output = *result.get();

    // FIXME: implement IceModelVec2Stag::copy_from()

    IceModelVec::AccessList list{&input, &output};

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
  : GeometryEvolution(grid) {

  m_no_model_mask.create(m_grid, "no_model_mask", WITH_GHOSTS);
  m_no_model_mask.set_attrs("model_mask", "'no model' mask", "", "");
}

void RegionalGeometryEvolution::set_no_model_mask_impl(const IceModelVec2Int &mask) {
  m_no_model_mask.copy_from(mask);
}

/*!
 * Disable ice flow in "no model" areas.
 */
void RegionalGeometryEvolution::compute_interface_fluxes(const IceModelVec2CellType &cell_type,
                                                         const IceModelVec2S        &ice_thickness,
                                                         const IceModelVec2V        &velocity,
                                                         const IceModelVec2Int      &velocity_bc_mask,
                                                         const IceModelVec2Stag     &diffusive_flux,
                                                         IceModelVec2Stag           &output) {

  GeometryEvolution::compute_interface_fluxes(cell_type, ice_thickness,
                                              velocity, velocity_bc_mask, diffusive_flux,
                                              output);

  IceModelVec::AccessList list{&m_no_model_mask, &output};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const int M = m_no_model_mask.as_int(i, j);

      for (unsigned int n = 0; n < 2; ++n) {
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
                                                                       const IceModelVec2Int      &thickness_bc_mask,
                                                                       const IceModelVec2S        &ice_thickness,
                                                                       const IceModelVec2CellType &cell_type,
                                                                       const IceModelVec2S        &surface_mass_flux,
                                                                       const IceModelVec2S        &basal_melt_rate,
                                                                       IceModelVec2S              &effective_SMB,
                                                                       IceModelVec2S              &effective_BMB) {
  GeometryEvolution::compute_surface_and_basal_mass_balance(dt,
                                                            thickness_bc_mask,
                                                            ice_thickness,
                                                            cell_type,
                                                            surface_mass_flux,
                                                            basal_melt_rate,
                                                            effective_SMB,
                                                            effective_BMB);

  IceModelVec::AccessList list{&m_no_model_mask, &effective_SMB, &effective_BMB};

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

void GeometryEvolution::set_no_model_mask(const IceModelVec2Int &mask) {
  this->set_no_model_mask_impl(mask);
}

void GeometryEvolution::set_no_model_mask_impl(const IceModelVec2Int &mask) {
  (void) mask;
  // the default implementation is a no-op
}


} // end of namespace pism
