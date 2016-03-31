// Copyright (C) 2004-2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cmath>
#include <cstring>
#include <petscdmda.h>
#include <cassert>
#include <algorithm>

#include "iceModel.hh"
#include "base/calving/PISMIcebergRemover.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_const.hh"
#include "coupler/PISMOcean.hh"
#include "coupler/PISMSurface.hh"
#include "earth/PISMBedDef.hh"
#include "base/util/pism_utilities.hh"

#include "base/grounded_cell_fraction.hh"

namespace pism {

using namespace mask;

//! \file iMgeometry.cc Methods of IceModel which update and maintain consistency of ice sheet geometry.


//! Update the surface elevation and the flow-type mask when the geometry has changed.
/*!
  This procedure should be called whenever necessary to maintain consistency of geometry.

  For instance, it should be called when either ice thickness or bed elevation change.
  In particular we always want \f$h = H + b\f$ to apply at grounded points, and, on the
  other hand, we want the mask to reflect that the ice is floating if the flotation
  criterion applies at a point.

  Also calls the code which removes icebergs, to avoid stress balance
  solver problems associated to not-attached-to-grounded ice.
*/
void IceModel::updateSurfaceElevationAndMask() {

  assert(m_ocean != NULL);
  const double sea_level = m_ocean->sea_level_elevation();

  const IceModelVec2S &bed_topography = m_beddef->bed_elevation();

  GeometryCalculator gc(*m_config);

  if (m_config->get_boolean("kill_icebergs") && iceberg_remover != NULL) {
    // the iceberg remover has to use the same mask as the stress balance code, hence the
    // stress-balance-related threshold here
    gc.set_icefree_thickness(m_config->get_double("mask_icefree_thickness_stress_balance_standard"));

    gc.compute_mask(sea_level, bed_topography, m_ice_thickness, m_cell_type);

    iceberg_remover->update(m_cell_type, m_ice_thickness);
    // the call above modifies ice thickness and updates the mask accordingly, but we re-compute the
    // mask (we need to use the different threshold)
  }

  gc.set_icefree_thickness(m_config->get_double("mask_icefree_thickness_standard"));
  gc.compute_mask(sea_level, bed_topography, m_ice_thickness, m_cell_type);
  gc.compute_surface(sea_level, bed_topography, m_ice_thickness, m_ice_surface_elevation);

  check_minimum_ice_thickness();
}

//! \brief Adjust ice flow through interfaces of the cell i,j.
/*!
 *
 * From the point of view of the code solving the mass continuity equation
 * there are 4 kinds of cells: "grounded ice", "floating ice", "ice-free land",
 * and "ice-free ocean". This makes 16 kinds of interfaces between grid cells.

 * In some of the cases PISM computes a SIA flux (or the SSA velocity) across
 * an interface that ice should not cross. (A silly example: there should be no
 * flow between two adjacent ice-free cells, even if a stress balance solver
 * computes velocity components for the whole computational domain. Please see
 * comments below for more examples.)

 * Note that 6 cases are paired. This is crucial: to be consistent (and
 * conserve mass) the interface between A and B has to get the same treatment
 * as the one between B and A.
 *
 * @param[in] mask Mask used to check for icy/ice-free and floatation
 * @param[in,out] SSA_velocity SSA velocity to be adjusted
 * @param[in,out] SIA_flux SIA flux to be adjusted
 */
void IceModel::adjust_flow(StarStencil<int> mask,
                           StarStencil<double> &SSA_velocity,
                           StarStencil<double> &SIA_flux) {

  // Prepare to loop over neighbors:
  // directions
  Direction dirs[] = {North, East, South, West};
  // mask values
  int mask_current = mask.ij;

  // Loop over neighbors:
  for (int n = 0; n < 4; ++n) {
    Direction direction = dirs[n];
    int mask_neighbor = mask[direction];

    // Only one of the cases below applies at any given time/location, so
    // "continuing" the for-loop allows us to avoid checking conditions we know
    // will fail. This also means that the code executed in *all* cases should
    // go here and not after if-conditions.

    // Case 1: Flow between grounded_ice and grounded_ice.
    if (grounded_ice(mask_current) && grounded_ice(mask_neighbor)) {
      // no adjustment; kept for completeness
      continue;
    }


    // Cases 2 and 3: Flow between grounded_ice and floating_ice.
    if ((grounded_ice(mask_current) && floating_ice(mask_neighbor)) ||
        (floating_ice(mask_current) && grounded_ice(mask_neighbor))) {
      // no adjustment
      continue;
    }


    // Cases 4 and 5: Flow between grounded_ice and ice_free_land.
    if ((grounded_ice(mask_current) && ice_free_land(mask_neighbor)) ||
        (ice_free_land(mask_current) && grounded_ice(mask_neighbor))) {
      // no adjustment
      continue;
    }


    // Cases 6 and 7: Flow between grounded_ice and ice_free_ocean.
    if ((grounded_ice(mask_current) && ice_free_ocean(mask_neighbor)) ||
        (ice_free_ocean(mask_current) && grounded_ice(mask_neighbor))) {
      // no adjustment
      continue;
    }

    // Case 8: Flow between floating_ice and floating_ice.
    if (floating_ice(mask_current) && floating_ice(mask_neighbor)) {
      // disable SIA flow
      SIA_flux[direction] = 0.0;
      continue;
    }


    // Cases 9 and 10: Flow between floating_ice and ice_free_land.
    if ((floating_ice(mask_current) && ice_free_land(mask_neighbor)) ||
        (ice_free_land(mask_current) && floating_ice(mask_neighbor))) {
      // disable all flow

      // This ensures that an ice shelf does not climb up a cliff.

      SIA_flux[direction] = 0.0;
      SSA_velocity[direction] = 0.0;
      continue;
    }


    // Cases 11 and 12: Flow between floating_ice and ice_free_ocean.
    if ((floating_ice(mask_current) && ice_free_ocean(mask_neighbor)) ||
        (ice_free_ocean(mask_current) && floating_ice(mask_neighbor))) {

      SIA_flux[direction] = 0.0;

      // The SSA flow may be later adjusted by the code implementing the
      // partially-filled cell parameterization.
      continue;
    }

    // Case 13: Flow between ice_free_land and ice_free_land.
    if (ice_free_land(mask_current) && ice_free_land(mask_neighbor)) {

      SIA_flux[direction] = 0.0;
      SSA_velocity[direction] = 0.0;
      continue;
    }


    // Cases 14 and 15: Flow between ice_free_land and ice_free_ocean.
    if ((ice_free_land(mask_current) && ice_free_ocean(mask_neighbor)) ||
        (ice_free_ocean(mask_current) && ice_free_land(mask_neighbor))) {

      SIA_flux[direction] = 0.0;
      SSA_velocity[direction] = 0.0;
      continue;
    }

    // Case 16: Flow between ice_free_ocean and ice_free_ocean.
    if (ice_free_ocean(mask_current) && ice_free_ocean(mask_neighbor)) {

      SIA_flux[direction] = 0.0;
      SSA_velocity[direction] = 0.0;
      continue;
    }
  } // end of the loop over neighbors (n)

}

//! \brief Compute fluxes through interfaces of a cell i,j.
/*!
 * This method implements two steps:
 *
 * 1) Compute SSA velocities through interfaces of a cell by averaging values
 * from regular grid neighbors, making sure that velocities from ice-free areas
 * are not used.
 *
 * Note that the input parameter `in_SSA_velocity` contains both components of
 * the velocity field in the neighborhood of i,j, while `out_SSA_velocity`
 * contains \b scalars: projections of velocity vectors onto normals to
 * corresponding cell interfaces.
 *
 * The SIA flux `in_SIA_flux` is computed on the staggered grid by SIAFD, so
 * the loop below just copies it to `out_SIA_flux`.
 *
 * 2) Adjust the flow using the mask by calling adjust_flow().
 *
 * @param[in] dirichlet_bc true if Dirichlet B.C. are set.
 * @param[in] i i-index of the current cell
 * @param[in] j j-index of the current cell
 * @param[in] in_SSA_velocity SSA velocity on the regular grid in the neighborhood of i,j
 * @param[in] in_SIA_flux SIA flux on the staggered grid (at interfaces of the cell i,j)
 * @param[out] out_SSA_velocity SSA velocities through interfaces of the cell i,j
 * @param[out] out_SIA_flux SIA flux through interfaces of the cell i,j
 */
void IceModel::cell_interface_fluxes(bool dirichlet_bc,
                                     int i, int j,
                                     StarStencil<Vector2> in_SSA_velocity,
                                     StarStencil<double> in_SIA_flux,
                                     StarStencil<double> &out_SSA_velocity,
                                     StarStencil<double> &out_SIA_flux) {

  StarStencil<int> mask = m_cell_type.int_star(i,j);
  Direction dirs[4] = {North, East, South, West};

  StarStencil<int> bc_mask;
  StarStencil<Vector2> bc_velocity;
  if (dirichlet_bc) {
    bc_mask = m_ssa_dirichlet_bc_mask.int_star(i,j);
    bc_velocity = m_ssa_dirichlet_bc_values.star(i,j);
  }

  out_SSA_velocity.ij = 0.0;
  out_SIA_flux.ij = 0.0;

  for (int n = 0; n < 4; ++n) {
    Direction direction = dirs[n];
    int mask_current = mask.ij,
      mask_neighbor = mask[direction];

    // The in_SIA_flux is already on the staggered grid, so we can just
    // copy it to out_SIA_flux:
    out_SIA_flux[direction] = in_SIA_flux[direction];

    // Compute the out_SSA_velocity (SSA):
    if (icy(mask_current) && icy(mask_neighbor)) {

      // Case 1: both sides of the interface are icy
      if (direction == East || direction == West) {
        out_SSA_velocity[direction] = 0.5 * (in_SSA_velocity.ij.u + in_SSA_velocity[direction].u);
      } else {
        out_SSA_velocity[direction] = 0.5 * (in_SSA_velocity.ij.v + in_SSA_velocity[direction].v);
      }

    } else if (icy(mask_current) && ice_free(mask_neighbor)) {

      // Case 2: icy cell next to an ice-free cell
      if (direction == East || direction == West) {
        out_SSA_velocity[direction] = in_SSA_velocity.ij.u;
      } else {
        out_SSA_velocity[direction] = in_SSA_velocity.ij.v;
      }

    } else if (ice_free(mask_current) && icy(mask_neighbor)) {

      // Case 3: ice-free cell next to icy cell
      if (direction == East || direction == West) {
        out_SSA_velocity[direction] = in_SSA_velocity[direction].u;
      } else {
        out_SSA_velocity[direction] = in_SSA_velocity[direction].v;
      }

    } else if (ice_free(mask_current) && ice_free(mask_neighbor)) {

      // Case 4: both sides of the interface are ice-free
      out_SSA_velocity[direction] = 0.0;

    }

    // The Dirichlet B.C. case:
    if (dirichlet_bc) {

      if (bc_mask.ij == 1 && bc_mask[direction] == 1) {

        // Case 1: both sides of the interface are B.C. locations: average from
        // the regular grid onto the staggered grid.
        if (direction == East || direction == West) {
          out_SSA_velocity[direction] = 0.5 * (bc_velocity.ij.u + bc_velocity[direction].u);
        } else {
          out_SSA_velocity[direction] = 0.5 * (bc_velocity.ij.v + bc_velocity[direction].v);
        }

      } else if (bc_mask.ij == 1 && bc_mask[direction] == 0) {

        // Case 2: at a Dirichlet B.C. location
        if (direction == East || direction == West) {
          out_SSA_velocity[direction] = bc_velocity.ij.u;
        } else {                    // North or South
          out_SSA_velocity[direction] = bc_velocity.ij.v;
        }

      } else if (bc_mask.ij == 0 && bc_mask[direction] == 1) {

        // Case 3: next to a Dirichlet B.C. location
        if (direction == East || direction == West) {
          out_SSA_velocity[direction] = bc_velocity[direction].u;
        } else {                  // North or South
          out_SSA_velocity[direction] = bc_velocity[direction].v;
        }

      } else {
        // Case 4: elsewhere.
        // No Dirichlet B.C. adjustment here.
      }

    } // end of "if (dirichlet_bc)"

  } // end of the loop over neighbors

  adjust_flow(mask, out_SSA_velocity, out_SIA_flux);

}


//! Update the thickness from the diffusive flux and sliding velocity, and the surface and basal mass balance rates.
/*!
  The partial differential equation describing the conservation of mass in the
  map-plane (parallel to the geoid) is
  \f[ \frac{\partial H}{\partial t} = M - S - \nabla\cdot \mathbf{q} \f]
  where
  \f[ \mathbf{q} = \bar{\mathbf{U}} H. \f]
  In these equations \f$H\f$ is the ice thickness,
  \f$M\f$ is the surface mass balance (accumulation or ablation), \f$S\f$ is the
  basal mass balance (e.g. basal melt or freeze-on), and \f$\bar{\mathbf{U}}\f$ is
  the vertically-averaged horizontal velocity of the ice.  This procedure uses
  this conservation of mass equation to update the ice thickness.

  The SurfaceModel IceModel::surface provides \f$M\f$.  The
  OceanModel IceModel::ocean provides \f$S\f$ at locations below
  floating ice (ice shelves).

  Even if we regard the map-plane flux as defined by the formula
  \f$\mathbf{q} = \bar{\mathbf{U}} H\f$, the flow of ice is at least somewhat
  diffusive in almost all cases.  In the non-sliding SIA model it
  is exactly true that \f$\mathbf{q} = - D \nabla h\f$.  In the current method the
  flux is split into the part from the diffusive non-sliding SIA model
  and a part which is a less-diffusive, presumably membrane-stress-dominated
  2D advective velocity, which generally describes sliding:
  \f[ \mathbf{q} = - D \nabla h + \mathbf{U}_b H.\f]
  Here \f$D\f$ is the (positive, scalar) effective diffusivity of the non-sliding
  SIA and \f$\mathbf{U}_b\f$ is the less-diffusive sliding velocity.
  We interpret \f$\mathbf{U}_b\f$ as a basal sliding velocity in the hybrid.

  The main ice-dynamical inputs to this method come from the outputs from
  StressBalance *stress_balance:
  \code
  const IceModelVec2Stag &Qdiff = stress_balance->diffusive_flux();
  const IceModelVec2V &vel_advective = stress_balance->advective_velocity();
  \endcode
  The diffusive flux \f$-D\nabla h\f$ is thus stored in a `IceModelVec2Stag`
  while the less-diffusive velocity \f$\mathbf{U}_b\f$ is stored in a
  `IceModelVec2V`.

  The methods used here are first-order and explicit in time.  The derivatives in
  \f$\nabla \cdot (D \nabla h)\f$ are computed by centered finite difference
  methods.  The diffusive flux `Qdiff` is already stored on the staggered grid
  and it is differenced in a centered way here.  The time-stepping for this part
  of the explicit scheme is controlled by equation (25) in [\ref BBL], so that
  \f$\Delta t \sim \Delta x^2 / \max D\f$; see also [\ref MortonMayers].

  The divergence of the flux from velocity \f$\mathbf{U}_b\f$ is computed by
  the upwinding technique [equation (25) in \ref Winkelmannetal2011] which
  is the donor cell upwind (i.e. Gudunov) method [\ref LeVeque].
  The CFL condition for this advection scheme is checked; see
  computeMax2DSlidingSpeed() and determineTimeStep().  This method implements the
  direct-superposition (PIK) hybrid which adds the SSA velocity to the SIA velocity
  [equation (15) in \ref Winkelmannetal2011].  The hybrid described by equations
  (21) and (22) in \ref BBL is no longer used.

We also compute total ice fluxes in kg s-1 at 3 interfaces:

  \li the ice-atmosphere interface: gets surface mass balance rate from
      SurfaceModel *surface,
  \li the ice-ocean interface at the bottom of ice shelves: gets ocean-imposed
      basal melt rate from OceanModel *ocean, and
  \li the ice-bedrock interface: gets basal melt rate from IceModelVec2S basal_melt_rate.

A unit-conversion occurs for all three quantities, from ice-equivalent m s-1
to kg s-1.  The sign convention about these fluxes is that positve flux means
ice is being \e added to the ice fluid volume at that interface.

These quantities should be understood as *instantaneous at the beginning of
the time-step.*  Multiplying by dt will \b not necessarily give the
corresponding change from the beginning to the end of the time-step.

FIXME:  The calving rate can be computed by post-processing:
dimassdt = surface_ice_flux + basal_ice_flux + sub_shelf_ice_flux + discharge_flux_mass_rate + nonneg_rule_flux

Removed commented-out code using the coverage ratio to compute the surface mass
balance contribution (to reduce clutter). Please see the commit 26330a7 and
earlier. (CK)

*/
void IceModel::massContExplicitStep() {

  double
    // totals over the processor's domain:
    proc_H_to_Href_flux           = 0,
    proc_Href_to_H_flux           = 0,
    proc_grounded_basal_ice_flux  = 0,
    proc_nonneg_rule_flux         = 0,
    proc_sub_shelf_ice_flux       = 0,
    proc_sum_divQ_SIA             = 0,
    proc_sum_divQ_SSA             = 0,
    proc_surface_ice_flux         = 0,
    // totals over all processors:
    total_H_to_Href_flux          = 0,
    total_Href_to_H_flux          = 0,
    total_grounded_basal_ice_flux = 0,
    total_nonneg_rule_flux        = 0,
    total_sub_shelf_ice_flux      = 0,
    total_sum_divQ_SIA            = 0,
    total_sum_divQ_SSA            = 0,
    total_surface_ice_flux        = 0;

  const bool
    include_bmr_in_continuity = m_config->get_boolean("include_bmr_in_continuity");

  const double
    dx                       = m_grid->dx(),
    dy                       = m_grid->dy(),
    ice_density              = m_config->get_double("ice_density"),
    meter_per_s_to_kg_per_m2 = m_dt * ice_density;

  assert(m_surface != NULL);
  m_surface->ice_surface_mass_flux(m_climatic_mass_balance);

  IceModelVec2S &H_new = vWork2d[0];
  H_new.copy_from(m_ice_thickness);

  IceModelVec2S &H_residual = vWork2d[1];

  const IceModelVec2Stag &Qdiff = m_stress_balance->diffusive_flux();

  const IceModelVec2V &vel_advective = m_stress_balance->advective_velocity();

  const IceModelVec2S &bed_topography = m_beddef->bed_elevation();

  IceModelVec::AccessList list;
  list.add(m_cell_area);
  list.add(m_ice_thickness);
  list.add(m_ice_surface_elevation);
  list.add(bed_topography);
  list.add(m_basal_melt_rate);
  list.add(Qdiff);
  list.add(vel_advective);
  list.add(m_climatic_mass_balance);
  list.add(m_cell_type);
  list.add(H_new);

  // related to PIK part_grid mechanism; see Albrecht et al 2011
  const bool do_part_grid = m_config->get_boolean("part_grid"),
    do_redist = m_config->get_boolean("part_redist"),
    reduce_frontal_thickness = m_config->get_boolean("part_grid_reduce_frontal_thickness");
  if (do_part_grid) {
    list.add(vHref);
    if (do_redist) {
      list.add(H_residual);
      // FIXME: next line causes mass loss if max_loopcount in redistResiduals()
      //        was not sufficient to zero-out H_residual already
      H_residual.set(0.0);
    }
  }
  const bool dirichlet_bc = m_config->get_boolean("ssa_dirichlet_bc");
  if (dirichlet_bc) {
    list.add(m_ssa_dirichlet_bc_mask);
    list.add(m_ssa_dirichlet_bc_values);
  }

  list.add(m_climatic_mass_balance_cumulative);
  list.add(m_nonneg_flux_2D_cumulative);
  list.add(m_grounded_basal_flux_2D_cumulative);
  list.add(m_floating_basal_flux_2D_cumulative);
  list.add(m_flux_divergence);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // These constants are used to convert ice equivalent
      // thicknesses and thickening rates to kg, for accounting of
      // fluxes during the current time-step.
      const double
        meter_to_kg       = m_cell_area(i,j) * ice_density,
        meter_per_s_to_kg = meter_to_kg * m_dt;

      // Divergence terms:
      double
        divQ_SIA = 0.0,         // units: [m s-1]
        divQ_SSA = 0.0;         // units: [m s-1]

      // Source terms:
      // Note: here we convert surface mass balance from [kg m-2 s-1] to [m s-1]:
      double
        surface_mass_balance = m_climatic_mass_balance(i, j) / ice_density, // units: [m s-1]
        basal_melt_rate      = 0.0, // units: [m s-1]
        H_to_Href_flux       = 0.0, // units: [m]
        Href_to_H_flux       = 0.0, // units: [m]
        nonneg_rule_flux     = 0.0; // units: [m]

      if (include_bmr_in_continuity) {
        basal_melt_rate = m_basal_melt_rate(i, j);
      }

      StarStencil<double> Q, v;
      cell_interface_fluxes(dirichlet_bc, i, j,
                            vel_advective.star(i, j), Qdiff.star(i, j),
                            v, Q);

      // Compute divergence terms:
      {
        // Staggered grid Div(Q) for diffusive non-sliding SIA deformation part:
        // divQ_SIA = - D grad h
        divQ_SIA = (Q.e - Q.w) / dx + (Q.n - Q.s) / dy;

        // Plug flow part (i.e. basal sliding; from SSA): upwind by staggered grid
        // PIK method;  this is   \nabla \cdot [(u, v) H]
        divQ_SSA += (v.e * (v.e > 0 ? m_ice_thickness(i, j) : m_ice_thickness(i + 1, j))
                     - v.w * (v.w > 0 ? m_ice_thickness(i - 1, j) : m_ice_thickness(i, j))) / dx;
        divQ_SSA += (v.n * (v.n > 0 ? m_ice_thickness(i, j) : m_ice_thickness(i, j + 1))
                     - v.s * (v.s > 0 ? m_ice_thickness(i, j - 1) : m_ice_thickness(i, j))) / dy;
      }

      // Set source terms

      if (m_cell_type.ice_free_ocean(i, j)) {
        // Decide whether to apply Albrecht et al 2011 subgrid-scale
        // parameterization
        if (do_part_grid && m_cell_type.next_to_ice(i, j)) {

          // Add the flow contribution to this partially filled cell.
          H_to_Href_flux  = -(divQ_SSA + divQ_SIA) * m_dt;
          vHref(i, j)    += H_to_Href_flux;

          if (vHref(i, j) < 0) {
            m_log->message(2,
                       "PISM WARNING: negative Href at (%d,%d)\n",
                       i, j);

            // Note: this adds mass!
            nonneg_rule_flux += vHref(i, j);
            vHref(i, j) = 0;
          }

          double H_threshold = part_grid_threshold_thickness(m_cell_type.int_star(i, j),
                                                             m_ice_thickness.star(i, j),
                                                             m_ice_surface_elevation.star(i, j),
                                                             bed_topography(i,j),
                                                             dx,
                                                             reduce_frontal_thickness);
          double coverage_ratio = 1.0;
          if (H_threshold > 0.0) {
            coverage_ratio = vHref(i, j) / H_threshold;
          }

          if (coverage_ratio >= 1.0) {
            // A partially filled grid cell is now considered to be full.
            if (do_redist) {
              H_residual(i, j) = vHref(i, j) - H_threshold; // residual ice thickness
            }

            vHref(i, j)    = 0.0;
            Href_to_H_flux = H_threshold;

            // A cell that became "full" experiences both SMB and basal melt.
          } else {
            // A not-full partially-filled cell experiences neither.
            surface_mass_balance = 0.0;
            basal_melt_rate      = 0.0;
          }

          // In this case the SSA flux goes into the Href variable and does not
          // directly contribute to ice thickness at this location.
          proc_sum_divQ_SIA += - divQ_SIA * meter_per_s_to_kg;
          proc_sum_divQ_SSA += - divQ_SSA * meter_per_s_to_kg;
          divQ_SIA           = 0.0;
          divQ_SSA           = 0.0;
        } else { // end of "if (part_grid...)

          // Standard ice-free ocean case:
          surface_mass_balance = 0.0;
          basal_melt_rate      = 0.0;
        }
      } // end of "if (ice_free_ocean)"

      // Dirichlet BC case (should go last to override previous settings):
      if (dirichlet_bc && m_ssa_dirichlet_bc_mask.as_int(i,j) == 1) {
        surface_mass_balance = 0.0;
        basal_melt_rate      = 0.0;
        Href_to_H_flux       = 0.0;
        divQ_SIA             = 0.0;
        divQ_SSA             = 0.0;
      }

      m_flux_divergence(i, j) = divQ_SIA + divQ_SSA;

      // mass transport
      H_new(i, j) += - m_dt * (divQ_SIA + divQ_SSA) + Href_to_H_flux;

      if (H_new(i, j) < 0.0) {
        nonneg_rule_flux += -H_new(i, j);

        // convert from [m] to [kg m-2]:
        m_nonneg_flux_2D_cumulative(i, j) += nonneg_rule_flux * ice_density; // units: [kg m-2]

        // this has to go *after* accounting above!
        H_new(i, j) = 0.0;
      }

      // surface mass balance
      if (H_new(i, j) + m_dt * surface_mass_balance < 0.0) {
        // applying the surface mass balance results in negative thickness
        //
        // modify the surface mass balance so that the resulting thickness is zero
        surface_mass_balance = - H_new(i, j) / m_dt;
        H_new(i, j)          = 0.0;
      } else {
        H_new(i, j) += m_dt * surface_mass_balance;
      }

      // basal mass balance
      if (H_new(i, j) - m_dt * basal_melt_rate < 0.0) {
        // applying the basal melt rate results in negative thickness
        //
        // modify the basal melt rate so that the resulting thickness is zero
        basal_melt_rate = H_new(i, j) / m_dt;
        H_new(i, j)     = 0.0;
      } else {
        H_new(i, j) += - m_dt * basal_melt_rate;
      }

      // surface_mass_balance has the units of [m s-1]; convert to [kg m-2]
      m_climatic_mass_balance_cumulative(i, j) += surface_mass_balance * meter_per_s_to_kg_per_m2;

      // basal_melt_rate has the units of [m s-1]; convert to [kg m-2]
      m_grounded_basal_flux_2D_cumulative(i, j) += -basal_melt_rate * meter_per_s_to_kg_per_m2;

      // basal_melt_rate has the units of [m s-1]; convert to [kg m-2]
      m_floating_basal_flux_2D_cumulative(i, j) += -basal_melt_rate * meter_per_s_to_kg_per_m2;

      // time-series accounting:
      {
        // all these are in units of [kg]
        if (m_cell_type.grounded(i,j)) {
          proc_grounded_basal_ice_flux += - basal_melt_rate * meter_per_s_to_kg;
        } else {
          proc_sub_shelf_ice_flux      += - basal_melt_rate * meter_per_s_to_kg;
        }

        proc_surface_ice_flux        +=   surface_mass_balance * meter_per_s_to_kg;
        proc_sum_divQ_SIA            += - divQ_SIA             * meter_per_s_to_kg;
        proc_sum_divQ_SSA            += - divQ_SSA             * meter_per_s_to_kg;
        proc_nonneg_rule_flux        +=   nonneg_rule_flux     * meter_to_kg;
        proc_H_to_Href_flux          += - H_to_Href_flux       * meter_to_kg;
        proc_Href_to_H_flux          +=   Href_to_H_flux       * meter_to_kg;
      }

    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  // flux accounting
  {
    // combine data to perform one reduction call instead of 8:
    double tmp_local[8] = {proc_grounded_basal_ice_flux,
                           proc_nonneg_rule_flux,
                           proc_sub_shelf_ice_flux,
                           proc_surface_ice_flux,
                           proc_sum_divQ_SIA,
                           proc_sum_divQ_SSA,
                           proc_Href_to_H_flux,
                           proc_H_to_Href_flux};
    double tmp_global[8];
    GlobalSum(m_grid->com, tmp_local, tmp_global, 8);

    total_grounded_basal_ice_flux = tmp_global[0];
    total_nonneg_rule_flux        = tmp_global[1];
    total_sub_shelf_ice_flux      = tmp_global[2];
    total_surface_ice_flux        = tmp_global[3];
    total_sum_divQ_SIA            = tmp_global[4];
    total_sum_divQ_SSA            = tmp_global[5];
    total_Href_to_H_flux          = tmp_global[6];
    total_H_to_Href_flux          = tmp_global[7];

    grounded_basal_ice_flux_cumulative += total_grounded_basal_ice_flux;
    sub_shelf_ice_flux_cumulative      += total_sub_shelf_ice_flux;
    surface_ice_flux_cumulative        += total_surface_ice_flux;
    sum_divQ_SIA_cumulative            += total_sum_divQ_SIA;
    sum_divQ_SSA_cumulative            += total_sum_divQ_SSA;
    nonneg_rule_flux_cumulative        += total_nonneg_rule_flux;
    Href_to_H_flux_cumulative          += total_Href_to_H_flux;
    H_to_Href_flux_cumulative          += total_H_to_Href_flux;
  }

  // finally copy H_new into ice_thickness and communicate ghosted values
  H_new.update_ghosts(m_ice_thickness);

  // distribute residual ice mass if desired
  if (do_redist) {
    residual_redistribution(H_residual);
  }

  // Check if the ice thickness exceeded the height of the computational box and stop if it did.
  check_maximum_ice_thickness();
}

/**
   @brief Updates the fractional "flotation mask".
 */
void IceModel::update_grounded_cell_fraction() {

  const double
    ice_density   = m_config->get_double("ice_density"),
    ocean_density = m_config->get_double("sea_water_density");

  assert(m_ocean != NULL);
  const double sea_level = m_ocean->sea_level_elevation();

  assert(m_beddef != NULL);
  const IceModelVec2S &bed_topography = m_beddef->bed_elevation();

  compute_grounded_cell_fraction(ice_density, ocean_density, sea_level,
                          m_ice_thickness, bed_topography, m_cell_type,
                          m_gl_mask, NULL, NULL);
}

} // end of namespace pism
