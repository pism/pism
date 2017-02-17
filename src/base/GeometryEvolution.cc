/* Copyright (C) 2016, 2017 PISM Authors
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

#include "base/util/iceModelVec.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"

namespace pism {

struct GeometryEvolution::Impl {
  Impl(IceGrid::ConstPtr g);

  IceGrid::ConstPtr grid;

  //! Flux divergence (used to track thickness changes due to flow).
  IceModelVec2S flux_divergence;
  //! Conservation error due to enforcing non-negativity of ice thickness.
  IceModelVec2S conservation_error;
  //! Effective surface mass balance.
  IceModelVec2S effective_SMB;
  //! Effective basal mass balance.
  IceModelVec2S effective_BMB;
  //! Change in ice thickness during the last time step.
  IceModelVec2S thickness_change;
  //! Change in the ice area-specific volume during the last time step.
  IceModelVec2S ice_area_specific_volume_change;

  // Work space
  IceModelVec2Stag velocity, flux;

  double ice_thickness_threshold;
};

GeometryEvolution::Impl::Impl(IceGrid::ConstPtr g)
  : grid(g) {

  flux.create(grid, "flux", WITH_GHOSTS);
  flux.set_attrs("internal", "fluxes through cell interfaces (sides)", "m2 s-1", "");
  flux.metadata().set_string("glaciological_units", "m2 year-1");
  flux.write_in_glaciological_units = true;

  flux_divergence.create(grid, "flux_divergence", WITHOUT_GHOSTS);
  flux_divergence.set_attrs("internal", "flux divergence", "m s-1", "");

  conservation_error.create(grid, "conservation_error", WITHOUT_GHOSTS);
  conservation_error.set_attrs("internal",
                               "conservation error due to enforcing non-negativity of ice thickness",
                               "m", "");

  effective_SMB.create(grid, "effective_SMB", WITHOUT_GHOSTS);
  effective_SMB.set_attrs("internal", "effective surface mass balance rate", "m s-1", "");
  effective_SMB.metadata().set_string("glaciological_units", "m year-1");
  effective_SMB.write_in_glaciological_units = true;

  effective_BMB.create(grid, "effective_BMB", WITHOUT_GHOSTS);
  effective_BMB.set_attrs("internal", "effective basal mass balance rate", "m s-1", "");
  effective_BMB.metadata().set_string("glaciological_units", "m year-1");
  effective_BMB.write_in_glaciological_units = true;

  thickness_change.create(grid, "thickness_change", WITHOUT_GHOSTS);
  thickness_change.set_attrs("internal", "change in thickness during the last time step", "m", "");

  ice_area_specific_volume_change.create(grid, "ice_area_specific_volume_change", WITHOUT_GHOSTS);
  ice_area_specific_volume_change.set_attrs("interval",
                                            "change in area-specific volume during the last time step",
                                            "m3 / m2", "");
}

GeometryEvolution::GeometryEvolution(IceGrid::ConstPtr grid) {
  m_impl->grid = grid;
}

const IceModelVec2S& GeometryEvolution::flux_divergence() const {
  return m_impl->flux_divergence;
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

const IceModelVec2S& GeometryEvolution::conservation_error() const {
  return m_impl->conservation_error;
}


/*!
 * @param[in,out] geometry ice geometry
 * @param[in] dt time step, seconds
 * @param[in] advective_velocity advective (SSA) velocity
 * @param[in] diffusive_flux diffusive (SIA) flux
 * @param[in] velocity_bc_mask advective velocity Dirichlet B.C. mask
 * @param[in] velocity_bc_values advective velocity Dirichlet B.C. values
 * @param[in] thickness_bc_mask ice thickness Dirichlet B.C. mask
 * @param[in] surface_mass_balance_rate top surface mass balance rate (m / second)
 * @param[in] basal_mass_balance_rate basal (bottom surface) mass balance rate
 */
void GeometryEvolution::step(Geometry &geometry, double dt,
                             const IceModelVec2V    &advective_velocity,
                             const IceModelVec2Stag &diffusive_flux,
                             const IceModelVec2Int  &velocity_bc_mask,
                             const IceModelVec2V    &velocity_bc_values,
                             const IceModelVec2Int  &thickness_bc_mask,
                             const IceModelVec2S    &surface_mass_balance_rate,
                             const IceModelVec2S    &basal_mass_balance_rate) {

  // uses B.C. data and includes modifications for regional runs
  compute_interface_velocity(geometry.cell_type(),
                             advective_velocity,
                             velocity_bc_mask,
                             velocity_bc_values,
                             m_impl->velocity);

  // upwinding shows up here; includes modifications for regional runs
  compute_interface_fluxes(geometry.cell_type(),
                           m_impl->velocity,
                           geometry.ice_thickness(),
                           diffusive_flux,
                           m_impl->flux);

  // simple finite differences
  compute_flux_divergence(m_impl->flux,
                          thickness_bc_mask,
                          m_impl->flux_divergence);

  // this is where part_grid should be implemented; uses thickness_bc_mask
  compute_thickness_change_due_to_flow(geometry.cell_type(),
                                       geometry.ice_thickness(),
                                       geometry.ice_surface_elevation(),
                                       geometry.bed_elevation(),
                                       m_impl->thickness_change,
                                       m_impl->ice_area_specific_volume_change);

  // Computes the numerical conservation error and corrects ice thickness. We can do this here
  // because compute_surface_and_basal_mass_balance() preserves non-negativity.
  ensure_nonnegativity(geometry.ice_thickness(),
                       geometry.ice_area_specific_volume(),
                       m_impl->thickness_change,
                       m_impl->ice_area_specific_volume_change,
                       m_impl->conservation_error);

  // ice extent may change due to the flow, so we need to update the mask in preparation for
  // applying the source terms (surface and basal mass balance)
  geometry.ensure_consistency(m_impl->ice_thickness_threshold);

  // uses thickness_bc_mask
  compute_surface_and_basal_mass_balance(dt,
                                         thickness_bc_mask,
                                         geometry.ice_thickness(),
                                         m_impl->thickness_change,
                                         geometry.cell_type(),
                                         surface_mass_balance_rate,
                                         basal_mass_balance_rate,
                                         m_impl->effective_SMB,
                                         m_impl->effective_BMB);

  // FIXME: Now H_new = H_old + thickness_change + effective_SMB + effective_BMB.

  geometry.ensure_consistency(m_impl->ice_thickness_threshold);

  // calving is a separate issue
}

static double limit_advective_velocity(int mask_current, int mask_neighbor, double velocity) {

  using mask::grounded_ice;
  using mask::floating_ice;
  using mask::ice_free_land;
  using mask::ice_free_ocean;

  // Case 1: Flow between grounded_ice and grounded_ice.
  if (grounded_ice(mask_current) && grounded_ice(mask_neighbor)) {
    return velocity;
  }

  // Cases 2 and 3: Flow between grounded_ice and floating_ice.
  if ((grounded_ice(mask_current) && floating_ice(mask_neighbor)) ||
      (floating_ice(mask_current) && grounded_ice(mask_neighbor))) {
    return velocity;
  }

  // Cases 4 and 5: Flow between grounded_ice and ice_free_land.
  if ((grounded_ice(mask_current) && ice_free_land(mask_neighbor)) ||
      (ice_free_land(mask_current) && grounded_ice(mask_neighbor))) {
    return velocity;
  }

  // Cases 6 and 7: Flow between grounded_ice and ice_free_ocean.
  if ((grounded_ice(mask_current) && ice_free_ocean(mask_neighbor)) ||
      (ice_free_ocean(mask_current) && grounded_ice(mask_neighbor))) {
    return velocity;
  }

  // Case 8: Flow between floating_ice and floating_ice.
  if (floating_ice(mask_current) && floating_ice(mask_neighbor)) {
    return velocity;
  }

  // Cases 9 and 10: Flow between floating_ice and ice_free_land.
  if ((floating_ice(mask_current) && ice_free_land(mask_neighbor)) ||
      (ice_free_land(mask_current) && floating_ice(mask_neighbor))) {
    // Disable all flow. This ensures that an ice shelf does not climb up a cliff.
    return 0.0;
  }

  // Cases 11 and 12: Flow between floating_ice and ice_free_ocean.
  if ((floating_ice(mask_current) && ice_free_ocean(mask_neighbor)) ||
      (ice_free_ocean(mask_current) && floating_ice(mask_neighbor))) {
    return velocity;
  }

  // Case 13: Flow between ice_free_land and ice_free_land.
  if (ice_free_land(mask_current) && ice_free_land(mask_neighbor)) {
    return 0.0;
  }

  // Cases 14 and 15: Flow between ice_free_land and ice_free_ocean.
  if ((ice_free_land(mask_current) && ice_free_ocean(mask_neighbor)) ||
      (ice_free_ocean(mask_current) && ice_free_land(mask_neighbor))) {
    return 0.0;
  }

  // Case 16: Flow between ice_free_ocean and ice_free_ocean.
  if (ice_free_ocean(mask_current) && ice_free_ocean(mask_neighbor)) {
    return 0.0;
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "cannot handle the mask_current=%d, mask_neighbor=%d case",
                                mask_current, mask_neighbor);
}

void GeometryEvolution::compute_interface_velocity(const IceModelVec2CellType &cell_type,
                                                   const IceModelVec2V &velocity,
                                                   const IceModelVec2Int &bc_mask,
                                                   const IceModelVec2V &bc_values,
                                                   IceModelVec2Stag &output) {

  using mask::icy;
  using mask::ice_free;

  IceModelVec::AccessList list{&cell_type, &velocity, &bc_mask, &bc_values, &output};

  ParallelSection loop(m_impl->grid->com);
  try {
    for (Points p(*m_impl->grid); p; p.next()) {
      const int
        i  = p.i(),
        j  = p.j(),
        M  = cell_type(i, j),
        BC = bc_mask.as_int(i, j);

      const Vector2
        V    = velocity(i, j),
        V_bc = bc_values(i, j);

      for (int n = 0; n < 2; ++n) {
        const int
          oi  = 1 - n,               // offset in the i direction
          oj  = n,                   // offset in the j direction
          i_n = i + oi,              // i index of a neighbor
          j_n = j + oj;              // j index of a neighbor

        const int M_n = cell_type(i_n, j_n); // mask of a neighbor

        // Regular case
        {
          const Vector2 V_n  = velocity(i_n, j_n);

          if (icy(M) and icy(M_n)) {
            // Case 1: both sides of the interface are icy
            output(i, j, n) = (n == 0 ? 0.5 * (V.u + V_n.u) : 0.5 * (V.v + V_n.v));

          } else if (icy(M) and ice_free(M_n)) {
            // Case 2: icy cell next to an ice-free cell
            output(i, j, n) = (n == 0 ? V.u : V.v);

          } else if (ice_free(M) and icy(M_n)) {
            // Case 3: ice-free cell next to icy cell
            output(i, j, n) = (n == 0 ? V_n.u : V_n.v);

          } else if (ice_free(M) and ice_free(M_n)) {
            // Case 4: both sides of the interface are ice-free
            output(i, j, n) = 0.0;

          }
        }

        // The Dirichlet B.C. case:
        {
          const int     BC_n   = bc_mask.as_int(i_n, j_n);
          const Vector2 V_bc_n = bc_values(i_n, j_n);

          if (BC == 1 and BC_n == 1) {
            // Case 1: both sides of the interface are B.C. locations: average from
            // the regular grid onto the staggered grid.
            output(i, j, n) = (n == 0 ? 0.5 * (V_bc.u + V_bc_n.u) : 0.5 * (V_bc.v + V_bc_n.v));

          } else if (BC == 1 and BC_n == 0) {
            // Case 2: at a Dirichlet B.C. location next to a regular location
            output(i, j, n) = (n == 0 ? V_bc.u : V_bc.v);

          } else if (BC == 0 and BC_n == 1) {

            // Case 3: at a regular location next to a Dirichlet B.C. location
            output(i, j, n) = (n == 0 ? V_bc_n.u : V_bc_n.v);

          } else {
            // Case 4: elsewhere.
            // No Dirichlet B.C. adjustment here.
          }

        } // end of the Dirichlet B.C. case

        // finally, limit advective velocities
        output(i, j, n) = limit_advective_velocity(M, M_n, output(i, j, n));

      } // staggered grid offset (n) loop

    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

static double limit_diffusive_flux(int mask_current, int mask_neighbor, double flux) {

  using mask::grounded_ice;
  using mask::floating_ice;
  using mask::ice_free_land;
  using mask::ice_free_ocean;

  // Case 1: Flow between grounded_ice and grounded_ice.
  if (grounded_ice(mask_current) && grounded_ice(mask_neighbor)) {
    return flux;
  }

  // Cases 2 and 3: Flow between grounded_ice and floating_ice.
  if ((grounded_ice(mask_current) && floating_ice(mask_neighbor)) ||
      (floating_ice(mask_current) && grounded_ice(mask_neighbor))) {
    return flux;
  }

  // Cases 4 and 5: Flow between grounded_ice and ice_free_land.
  if ((grounded_ice(mask_current) && ice_free_land(mask_neighbor)) ||
      (ice_free_land(mask_current) && grounded_ice(mask_neighbor))) {
    return flux;
  }

  // Cases 6 and 7: Flow between grounded_ice and ice_free_ocean.
  if ((grounded_ice(mask_current) && ice_free_ocean(mask_neighbor)) ||
      (ice_free_ocean(mask_current) && grounded_ice(mask_neighbor))) {
    return flux;
  }

  // Case 8: Flow between floating_ice and floating_ice.
  if (floating_ice(mask_current) && floating_ice(mask_neighbor)) {
    // no diffusive flux in ice shelves
    return 0.0;
  }

  // Cases 9 and 10: Flow between floating_ice and ice_free_land.
  if ((floating_ice(mask_current) && ice_free_land(mask_neighbor)) ||
      (ice_free_land(mask_current) && floating_ice(mask_neighbor))) {
    // Disable all flow. This ensures that an ice shelf does not climb up a cliff.
    return 0.0;
  }

  // Cases 11 and 12: Flow between floating_ice and ice_free_ocean.
  if ((floating_ice(mask_current) && ice_free_ocean(mask_neighbor)) ||
      (ice_free_ocean(mask_current) && floating_ice(mask_neighbor))) {
    return 0.0;
  }

  // Case 13: Flow between ice_free_land and ice_free_land.
  if (ice_free_land(mask_current) && ice_free_land(mask_neighbor)) {
    return 0.0;
  }

  // Cases 14 and 15: Flow between ice_free_land and ice_free_ocean.
  if ((ice_free_land(mask_current) && ice_free_ocean(mask_neighbor)) ||
      (ice_free_ocean(mask_current) && ice_free_land(mask_neighbor))) {
    return 0.0;
  }

  // Case 16: Flow between ice_free_ocean and ice_free_ocean.
  if (ice_free_ocean(mask_current) && ice_free_ocean(mask_neighbor)) {
    return 0.0;
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "cannot handle the mask_current=%d, mask_neighbor=%d case",
                                mask_current, mask_neighbor);
}

/*!
 * Combine advective velocity and the diffusive flux on the staggered grid with the ice thickness to
 * compute the total flux through cell interfaces.
 */
void GeometryEvolution::compute_interface_fluxes(const IceModelVec2CellType &cell_type,
                                                 const IceModelVec2Stag &velocity,
                                                 const IceModelVec2S &ice_thickness,
                                                 const IceModelVec2Stag &diffusive_flux,
                                                 IceModelVec2Stag &output) {

  ParallelSection loop(m_impl->grid->com);
  try {
    for (Points p(*m_impl->grid); p; p.next()) {
      const int
        i  = p.i(),
        j  = p.j(),
        M  = cell_type(i, j);

      for (int n = 0; n < 2; ++n) {
        const int
          oi  = 1 - n,               // offset in the i direction
          oj  = n,                   // offset in the j direction
          i_n = i + oi,              // i index of a neighbor
          j_n = j + oj;              // j index of a neighbor

        const int M_n = cell_type(i_n, j_n);

        // diffusive
        const double
          Q_diffusive = limit_diffusive_flux(M, M_n, diffusive_flux(i, j, n));

        // advective
        const double
          v           = velocity(i, j, n),
          H           = ice_thickness(i, j),
          H_n         = ice_thickness(i_n, j_n),
          Q_advective = v * (v > 0.0 ? H : H_n); // first order upwinding

        output(i, j, n) = Q_diffusive + Q_advective;
      } // end of the loop over neighbors (n)
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

void GeometryEvolution::compute_flux_divergence(const IceModelVec2Stag &flux,
                                                const IceModelVec2Int &thickness_bc_mask,
                                                IceModelVec2S &output) {
  const double
    dx = m_impl->grid->dx(),
    dy = m_impl->grid->dy();

  ParallelSection loop(m_impl->grid->com);
  try {
    for (Points p(*m_impl->grid); p; p.next()) {
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

void GeometryEvolution::compute_thickness_change_due_to_flow(const IceModelVec2CellType &cell_type,
                                                             const IceModelVec2S &ice_thickness,
                                                             const IceModelVec2S &ice_surface_elevation,
                                                             const IceModelVec2S &bed_elevation,
                                                             IceModelVec2S &thickness_change,
                                                             IceModelVec2S &area_specific_volume_change) {
  IceModelVec::AccessList list{&cell_type, &ice_thickness, &ice_surface_elevation,
      &bed_elevation, &thickness_change, &area_specific_volume_change};

  ParallelSection loop(m_impl->grid->com);
  try {
    for (Points p(*m_impl->grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // FIXME
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();


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
 */
void GeometryEvolution::ensure_nonnegativity(const IceModelVec2S &ice_thickness,
                                             const IceModelVec2S &area_specific_volume,
                                             IceModelVec2S &thickness_change,
                                             IceModelVec2S &area_specific_volume_change,
                                             IceModelVec2S &conservation_error) {
  ParallelSection loop(m_impl->grid->com);
  try {
    for (Points p(*m_impl->grid); p; p.next()) {
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
  return H + dH < 0.0 ? -H : dH;
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
 */
void GeometryEvolution::compute_surface_and_basal_mass_balance(double dt,
                                                               const IceModelVec2Int &thickness_bc_mask,
                                                               const IceModelVec2S &ice_thickness,
                                                               const IceModelVec2S &thickness_change,
                                                               const IceModelVec2CellType &cell_type,
                                                               const IceModelVec2S &smb_rate,
                                                               const IceModelVec2S &basal_melt_rate,
                                                               IceModelVec2S &effective_SMB,
                                                               IceModelVec2S &effective_BMB) {

  Config::ConstPtr config = m_impl->grid->ctx()->config();
  const double ice_density = config->get_double("constants.ice.density");

  const bool use_bmr = config->get_boolean("geometry.update.use_basal_melt_rate");

  IceModelVec::AccessList list{&ice_thickness, &thickness_change,
      &smb_rate, &basal_melt_rate,
      &cell_type, &thickness_bc_mask, &effective_SMB, &effective_BMB};

  ParallelSection loop(m_impl->grid->com);
  try {
    for (Points p(*m_impl->grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // Don't modify ice thickness at Dirichlet B.C. locations and in the ice-free ocean.
      if (thickness_bc_mask.as_int(i, j) == 1 or
          cell_type.ice_free_ocean(i, j)) {
        effective_SMB(i, j) = 0.0;
        effective_BMB(i, j) = 0.0;
        continue;
      }

      const double
        H       = ice_thickness(i, j),
        dH_flow = thickness_change(i, j);

      // Thickness change due to the surface mass balance
      //
      // Note that here we convert surface mass balance from [kg m-2 s-1] to [m s-1].
      const double dH_SMB = effective_change(H + dH_flow,
                                             dt * smb_rate(i, j) / ice_density);
      effective_SMB(i, j) = dH_SMB;

      // Thickness change due to the basal mass balance
      //
      // Note that basal_melt_rate is in [m s-1]. Here the negative sign converts the melt rate into
      // mass balance.
      effective_BMB(i, j) = effective_change(H + dH_flow + dH_SMB,
                                             dt * (use_bmr ? -basal_melt_rate(i, j) : 0.0));
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

} // end of namespace pism
