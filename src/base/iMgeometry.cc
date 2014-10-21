// Copyright (C) 2004-2014 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <assert.h>
#include <algorithm>

#include "iceModel.hh"
#include "Mask.hh"
#include "PISMOcean.hh"
#include "PISMSurface.hh"
#include "PISMStressBalance.hh"
#include "PISMIcebergRemover.hh"
#include "error_handling.hh"

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
PetscErrorCode IceModel::updateSurfaceElevationAndMask() {
  PetscErrorCode ierr;

  ierr = update_mask(bed_topography, ice_thickness, vMask); CHKERRQ(ierr);
  ierr = update_surface_elevation(bed_topography, ice_thickness, ice_surface_elevation); CHKERRQ(ierr);

  if (config.get_flag("kill_icebergs") && iceberg_remover != NULL) {
    ierr = iceberg_remover->update(vMask, ice_thickness); CHKERRQ(ierr);
    // the call above modifies ice thickness and updates the mask
    // accordingly
    ierr = update_surface_elevation(bed_topography, ice_thickness, ice_surface_elevation); CHKERRQ(ierr);
  }

  return 0;
}

/**
 * Update ice cover mask using the floatation criterion, sea level
 * elevation, ice thickness, and bed topography.
 *
 * @param[in]  bed bedrock surface elevation
 * @param[in]  thickness ice thicnness
 * @param[out] result cell type mask
 *
 * @return 0 on success.
 */
PetscErrorCode IceModel::update_mask(IceModelVec2S &bed,
                                     IceModelVec2S &thickness,
                                     IceModelVec2Int &result) {
  PetscErrorCode ierr;
  double sea_level;

  assert(ocean != NULL);
  ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);

  GeometryCalculator gc(sea_level, config);

  IceModelVec::AccessList list(result);
  list.add(bed);
  list.add(thickness);

  unsigned int GHOSTS = result.get_stencil_width();
  assert(bed.get_stencil_width() >= result.get_stencil_width());
  assert(thickness.get_stencil_width() >= result.get_stencil_width());

  for (PointsWithGhosts p(grid, GHOSTS); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i, j) = gc.mask(bed(i, j), thickness(i, j));
  }

  return 0;
}

/**
 * Update ice surface elevation using the floatation criterion, sea
 * level elevation, ice thickness, and bed topography.
 *
 * Uses ghosts of `bed_topography`, `ice_thickness`. Updates ghosts of
 * the `result`.
 *
 * @param[in] bed bedrock surface elevation
 * @param[in] thickness ice thickness
 * @param[out] result computed surface elevation
 *
 * @return 0 on success.
 */
PetscErrorCode IceModel::update_surface_elevation(IceModelVec2S &bed,
                                                  IceModelVec2S &thickness,
                                                  IceModelVec2S &result) {
  PetscErrorCode ierr;
  double sea_level;

  assert(ocean != NULL);
  ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);

  GeometryCalculator gc(sea_level, config);

  IceModelVec::AccessList list(result);
  list.add(bed);
  list.add(thickness);

  unsigned int GHOSTS = result.get_stencil_width();

  assert(bed.get_stencil_width() >= result.get_stencil_width());
  assert(thickness.get_stencil_width() >= result.get_stencil_width());

  for (PointsWithGhosts p(grid, GHOSTS); p; p.next()) {
    const int i = p.i(), j = p.j();

    // take this opportunity to check that thickness(i, j) >= 0
    if (thickness(i, j) < 0) {
      throw RuntimeError::formatted("Thickness negative at point i=%d, j=%d", i, j);
    }
    result(i, j) = gc.surface(bed(i, j), thickness(i, j));
  }

  return 0;
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
void IceModel::adjust_flow(planeStar<int> mask,
                           planeStar<double> &SSA_velocity,
                           planeStar<double> &SIA_flux) {

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
 * Note that the input parameter `input_velocity` contains both components of
 * the velocity field in the neighborhood of i,j, while `output_velocity`
 * contains \b scalars: projections of velocity vectors onto normals to
 * corresponding cell interfaces.
 *
 * The SIA flux `input_flux` is computed on the staggered grid by SIAFD, so
 * the loop below just copies it to `output_flux`.
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
                                     planeStar<Vector2> in_SSA_velocity,
                                     planeStar<double> in_SIA_flux,
                                     planeStar<double> &out_SSA_velocity,
                                     planeStar<double> &out_SIA_flux) {

  planeStar<int> mask = vMask.int_star(i,j);
  Direction dirs[4] = {North, East, South, West};

  planeStar<int> bc_mask;
  planeStar<Vector2> bc_velocity;
  if (dirichlet_bc) {
    bc_mask = vBCMask.int_star(i,j);
    bc_velocity = vBCvel.star(i,j);
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
      if (direction == East || direction == West)
        out_SSA_velocity[direction] = 0.5 * (in_SSA_velocity.ij.u + in_SSA_velocity[direction].u);
      else
        out_SSA_velocity[direction] = 0.5 * (in_SSA_velocity.ij.v + in_SSA_velocity[direction].v);

    } else if (icy(mask_current) && ice_free(mask_neighbor)) {

      // Case 2: icy cell next to an ice-free cell
      if (direction == East || direction == West)
        out_SSA_velocity[direction] = in_SSA_velocity.ij.u;
      else
        out_SSA_velocity[direction] = in_SSA_velocity.ij.v;

    } else if (ice_free(mask_current) && icy(mask_neighbor)) {

      // Case 3: ice-free cell next to icy cell
      if (direction == East || direction == West)
        out_SSA_velocity[direction] = in_SSA_velocity[direction].u;
      else
        out_SSA_velocity[direction] = in_SSA_velocity[direction].v;

    } else if (ice_free(mask_current) && ice_free(mask_neighbor)) {

      // Case 4: both sides of the interface are ice-free
      out_SSA_velocity[direction] = 0.0;

    }

    // The Dirichlet B.C. case:
    if (dirichlet_bc) {

      if (bc_mask.ij == 1 && bc_mask[direction] == 1) {

        // Case 1: both sides of the interface are B.C. locations: average from
        // the regular grid onto the staggered grid.
        if (direction == East || direction == West)
          out_SSA_velocity[direction] = 0.5 * (bc_velocity.ij.u + bc_velocity[direction].u);
        else
          out_SSA_velocity[direction] = 0.5 * (bc_velocity.ij.v + bc_velocity[direction].v);

      } else if (bc_mask.ij == 1 && bc_mask[direction] == 0) {

        // Case 2: at a Dirichlet B.C. location
        if (direction == East || direction == West)
          out_SSA_velocity[direction] = bc_velocity.ij.u;
        else                    // North or South
          out_SSA_velocity[direction] = bc_velocity.ij.v;

      } else if (bc_mask.ij == 0 && bc_mask[direction] == 1) {

        // Case 3: next to a Dirichlet B.C. location
        if (direction == East || direction == West)
          out_SSA_velocity[direction] = bc_velocity[direction].u;
        else                  // North or South
          out_SSA_velocity[direction] = bc_velocity[direction].v;

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

  The main ice-dynamical inputs to this method are identified in these lines,
  which get outputs from StressBalance *stress_balance:
  \code
  IceModelVec2Stag *Qdiff;
  stress_balance->get_diffusive_flux(Qdiff);
  IceModelVec2V *vel_advective;
  stress_balance->get_advective_2d_velocity(vel_advective);
  \endcode
  The diffusive flux \f$-D\nabla h\f$ is thus stored in `IceModelVec2Stag`
  `*Qdiff` while the less-diffusive velocity \f$\mathbf{U}_b\f$ is stored in
  `IceModelVec2V` `*vel_advective`.

  The methods used here are first-order and explicit in time.  The derivatives in
  \f$\nabla \cdot (D \nabla h)\f$ are computed by centered finite difference
  methods.  The diffusive flux `Qdiff` is already stored on the staggered grid
  and it is differenced in a centered way here.  The time-stepping for this part
  of the explicit scheme is controlled by equation (25) in [\ref BBL], so that
  \f$\Delta t \sim \Delta x^2 / \max D\f$; see also [\ref MortonMayers].

  The divergence of the flux from velocity \f$\mathbf{U}_b\f$ is computed by
  the upwinding technique [equation (25) in \ref Winkelmannetal2011] which
  is the donor cell upwind method [\ref LeVeque].
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
PetscErrorCode IceModel::massContExplicitStep() {
  PetscErrorCode ierr;

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

  const double dx = grid.dx, dy = grid.dy;
  bool
    include_bmr_in_continuity = config.get_flag("include_bmr_in_continuity"),
    compute_cumulative_climatic_mass_balance = climatic_mass_balance_cumulative.was_created(),
    compute_cumulative_nonneg_flux = nonneg_flux_2D_cumulative.was_created(),
    compute_cumulative_grounded_basal_flux = grounded_basal_flux_2D_cumulative.was_created(),
    compute_cumulative_floating_basal_flux = floating_basal_flux_2D_cumulative.was_created(),
    compute_flux_divergence = flux_divergence.was_created();

  double ice_density = config.get("ice_density"),
    meter_per_s_to_kg_per_m2 = dt * ice_density;

  assert(surface != NULL);
  ierr = surface->ice_surface_mass_flux(climatic_mass_balance); CHKERRQ(ierr);

  IceModelVec2S &vHnew = vWork2d[0];
  ierr = ice_thickness.copy_to(vHnew); CHKERRQ(ierr);

  IceModelVec2S &H_residual = vWork2d[1];

  IceModelVec2Stag *Qdiff;
  ierr = stress_balance->get_diffusive_flux(Qdiff); CHKERRQ(ierr);

  IceModelVec2V *vel_advective;
  ierr = stress_balance->get_2D_advective_velocity(vel_advective); CHKERRQ(ierr);

  IceModelVec::AccessList list;
  list.add(cell_area);
  list.add(ice_thickness);
  list.add(ice_surface_elevation);
  list.add(bed_topography);
  list.add(basal_melt_rate);
  list.add(*Qdiff);
  list.add(*vel_advective);
  list.add(climatic_mass_balance);
  list.add(vMask);
  list.add(vHnew);

  // related to PIK part_grid mechanism; see Albrecht et al 2011
  const bool do_part_grid = config.get_flag("part_grid"),
    do_redist = config.get_flag("part_redist"),
    reduce_frontal_thickness = config.get_flag("part_grid_reduce_frontal_thickness");
  if (do_part_grid) {
    list.add(vHref);
    if (do_redist) {
      list.add(H_residual);
      // FIXME: next line causes mass loss if max_loopcount in redistResiduals()
      //        was not sufficient to zero-out H_residual already
      ierr = H_residual.set(0.0); CHKERRQ(ierr);
    }
  }
  const bool dirichlet_bc = config.get_flag("ssa_dirichlet_bc");
  if (dirichlet_bc) {
    list.add(vBCMask);
    list.add(vBCvel);
  }

  if (compute_cumulative_climatic_mass_balance) {
    list.add(climatic_mass_balance_cumulative);
  }

  if (compute_cumulative_nonneg_flux) {
    list.add(nonneg_flux_2D_cumulative);
  }

  if (compute_cumulative_grounded_basal_flux) {
    list.add(grounded_basal_flux_2D_cumulative);
  }

  if (compute_cumulative_floating_basal_flux) {
    list.add(floating_basal_flux_2D_cumulative);
  }

  if (compute_flux_divergence) {
    list.add(flux_divergence);
  }

  MaskQuery mask(vMask);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // These constants are used to convert ice equivalent
    // thicknesses and thickening rates to kg, for accounting of
    // fluxes during the current time-step.
    const double
      meter_to_kg       = cell_area(i,j) * ice_density,
      meter_per_s_to_kg = meter_to_kg * dt;

    // Divergence terms:
    double
      divQ_SIA = 0.0,         // units: [m s-1]
      divQ_SSA = 0.0;         // units: [m s-1]

    // Source terms:
    // Note: here we convert surface mass balance from [kg m-2 s-1] to [m s-1]:
    double
      surface_mass_balance = climatic_mass_balance(i, j) / ice_density, // units: [m s-1]
      my_basal_melt_rate   = 0.0, // units: [m s-1]
      H_to_Href_flux       = 0.0, // units: [m]
      Href_to_H_flux       = 0.0, // units: [m]
      nonneg_rule_flux     = 0.0; // units: [m]

    if (include_bmr_in_continuity) {
      my_basal_melt_rate = basal_melt_rate(i, j);
    }

    planeStar<double> Q, v;
    cell_interface_fluxes(dirichlet_bc, i, j,
                          vel_advective->star(i, j), Qdiff->star(i, j),
                          v, Q);

    // Compute divergence terms:
    {
      // Staggered grid Div(Q) for diffusive non-sliding SIA deformation part:
      // divQ_SIA = - D grad h
      divQ_SIA = (Q.e - Q.w) / dx + (Q.n - Q.s) / dy;

      // Plug flow part (i.e. basal sliding; from SSA): upwind by staggered grid
      // PIK method;  this is   \nabla \cdot [(u, v) H]
      divQ_SSA += (v.e * (v.e > 0 ? ice_thickness(i, j) : ice_thickness(i + 1, j))
                   - v.w * (v.w > 0 ? ice_thickness(i - 1, j) : ice_thickness(i, j))) / dx;
      divQ_SSA += (v.n * (v.n > 0 ? ice_thickness(i, j) : ice_thickness(i, j + 1))
                   - v.s * (v.s > 0 ? ice_thickness(i, j - 1) : ice_thickness(i, j))) / dy;
    }

    // Set source terms

    if (mask.ice_free_ocean(i, j)) {
      // Decide whether to apply Albrecht et al 2011 subgrid-scale
      // parameterization
      if (do_part_grid && mask.next_to_ice(i, j)) {

        // Add the flow contribution to this partially filled cell.
        H_to_Href_flux = -(divQ_SSA + divQ_SIA) * dt;
        vHref(i, j) += H_to_Href_flux;
        if (vHref(i, j) < 0) {
          ierr = verbPrintf(3, grid.com,
                            "PISM WARNING: negative Href at (%d,%d)\n",
                            i, j); CHKERRQ(ierr);

          // Note: this adds mass!
          nonneg_rule_flux += vHref(i, j);
          vHref(i, j) = 0;
        }

        double H_threshold = get_threshold_thickness(vMask.int_star(i, j),
                                                     ice_thickness.star(i, j),
                                                     ice_surface_elevation.star(i, j),
                                                     bed_topography(i,j),
                                                     reduce_frontal_thickness);
        double coverage_ratio = 1.0;
        if (H_threshold > 0.0)
          coverage_ratio = vHref(i, j) / H_threshold;

        if (coverage_ratio >= 1.0) {
          // A partially filled grid cell is now considered to be full.
          if (do_redist)
            H_residual(i, j) = vHref(i, j) - H_threshold; // residual ice thickness

          vHref(i, j) = 0.0;

          Href_to_H_flux = H_threshold;

          // A cell that became "full" experiences both SMB and basal melt.
        } else {
          // An empty of partially-filled cell experiences neither.
          surface_mass_balance = 0.0;
          my_basal_melt_rate   = 0.0;
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
        my_basal_melt_rate   = 0.0;
      }
    } // end of "if (ice_free_ocean)"

      // Dirichlet BC case (should go last to override previous settings):
    if (dirichlet_bc && vBCMask.as_int(i,j) == 1) {
      surface_mass_balance = 0.0;
      my_basal_melt_rate   = 0.0;
      Href_to_H_flux       = 0.0;
      divQ_SIA             = 0.0;
      divQ_SSA             = 0.0;
    }

    if (compute_flux_divergence == true) {
      flux_divergence(i, j) = divQ_SIA + divQ_SSA;
    }

    vHnew(i, j) += (dt * (surface_mass_balance // accumulation/ablation
                          - my_basal_melt_rate // basal melt rate (grounded or floating)
                          - (divQ_SIA + divQ_SSA)) // flux divergence
                    + Href_to_H_flux); // corresponds to a cell becoming "full"

    if (vHnew(i, j) < 0.0) {
      nonneg_rule_flux += -vHnew(i, j);

      if (compute_cumulative_nonneg_flux)
        // convert from [m] to [kg m-2]:
        nonneg_flux_2D_cumulative(i, j) += nonneg_rule_flux * ice_density; // units: [kg m-2]

      // this has to go *after* accounting above!
      vHnew(i, j) = 0.0;
    }

    // Track cumulative surface mass balance. Note that this keeps track of
    // cumulative climatic_mass_balance at all the grid cells (including ice-free cells).
    if (compute_cumulative_climatic_mass_balance) {
      // surface_mass_balance has the units of [m s-1]; convert to [kg m-2]
      climatic_mass_balance_cumulative(i, j) += surface_mass_balance * meter_per_s_to_kg_per_m2;
    }

    if (compute_cumulative_grounded_basal_flux && mask.grounded(i, j)) {
      // my_basal_melt_rate has the units of [m s-1]; convert to [kg m-2]
      grounded_basal_flux_2D_cumulative(i, j) += -my_basal_melt_rate * meter_per_s_to_kg_per_m2;
    }

    if (compute_cumulative_floating_basal_flux && mask.ocean(i, j)) {
      // my_basal_melt_rate has the units of [m s-1]; convert to [kg m-2]
      floating_basal_flux_2D_cumulative(i, j) += -my_basal_melt_rate * meter_per_s_to_kg_per_m2;
    }

    // time-series accounting:
    {
      // all these are in units of [kg]
      if (mask.grounded(i,j)) {
        proc_grounded_basal_ice_flux += - my_basal_melt_rate * meter_per_s_to_kg;
      } else {
        proc_sub_shelf_ice_flux      += - my_basal_melt_rate * meter_per_s_to_kg;
      }

      proc_surface_ice_flux        +=   surface_mass_balance * meter_per_s_to_kg;
      proc_sum_divQ_SIA            += - divQ_SIA             * meter_per_s_to_kg;
      proc_sum_divQ_SSA            += - divQ_SSA             * meter_per_s_to_kg;
      proc_nonneg_rule_flux        +=   nonneg_rule_flux     * meter_to_kg;
      proc_H_to_Href_flux          += - H_to_Href_flux       * meter_to_kg;
      proc_Href_to_H_flux          +=   Href_to_H_flux       * meter_to_kg;
    }

  }

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
    ierr = GlobalSum(grid.com, tmp_local, tmp_global, 8); CHKERRQ(ierr);

    proc_grounded_basal_ice_flux = tmp_global[0];
    proc_nonneg_rule_flux        = tmp_global[1];
    proc_sub_shelf_ice_flux      = tmp_global[2];
    proc_surface_ice_flux        = tmp_global[3];
    proc_sum_divQ_SIA            = tmp_global[4];
    proc_sum_divQ_SSA            = tmp_global[5];
    proc_Href_to_H_flux          = tmp_global[6];
    proc_H_to_Href_flux          = tmp_global[7];

    grounded_basal_ice_flux_cumulative += total_grounded_basal_ice_flux;
    sub_shelf_ice_flux_cumulative      += total_sub_shelf_ice_flux;
    surface_ice_flux_cumulative        += total_surface_ice_flux;
    sum_divQ_SIA_cumulative            += total_sum_divQ_SIA;
    sum_divQ_SSA_cumulative            += total_sum_divQ_SSA;
    nonneg_rule_flux_cumulative        += total_nonneg_rule_flux;
    Href_to_H_flux_cumulative          += total_Href_to_H_flux;
    H_to_Href_flux_cumulative          += total_H_to_Href_flux;
  }

  // finally copy vHnew into ice_thickness and communicate ghosted values
  ierr = vHnew.update_ghosts(ice_thickness); CHKERRQ(ierr);

  // distribute residual ice mass if desired
  if (do_redist) {
    ierr = residual_redistribution(H_residual); CHKERRQ(ierr);
  }

  // Check if the ice thickness exceeded the height of the computational box
  // and extend the grid if necessary:
  ierr = check_maximum_thickness(); CHKERRQ(ierr);

  return 0;
}

/**
   @brief Updates the fractional "floatation mask".

   This mask ranges from 0 to 1 and is equal to the fraction of the
   cell (by area) that is grounded.

   Currently it is used to adjust the basal drag near the grounding
   line in the SSAFD stress balance model and the basal melt rate
   computation in the energy balance code (IceModel::temperatureStep()
   and IceModel::enthalpyAndDrainageStep()).

   We use the 1D (flow line) parameterization of the sub-grid
   grounding line position due to [@ref Gladstoneetal2010], (section
   3.1.1) and generalize it to the case of arbitrary sea level
   elevation. Then this sub-grid grounding line position is used to
   compute the grounded area fraction for each cell.

   We use the "LI" 1D grounding line position parameterization in "x"
   and "y" directions *independently*. The grounding line positions in
   the "x" and "y" directions (@f$ \lambda_x @f$ and @f$ \lambda_y
   @f$) are then interpreted as *width* and *height* of a rectangular
   sub-set of the cell. The grounded area is computed as the product
   of these two: @f$ A_g = \lambda_x \times \lambda_y. @f$

   Consider a cell at `(i,j)` and assume that the ice is grounded
   there and floating at `(i+1,j)`.

   Assume that the ice thickness and bedrock elevation change linearly
   from `(i,j)` to `(i+1,j)` and drop the `j` index for clarity:

   @f{align*}{
   H(\lambda) &= H_{i}(1 - \lambda) + H_{i+1}\lambda,\\
   b(\lambda) &= b_{i}(1 - \lambda) + b_{i+1}\lambda.\\
   @f}

   Here @f$ \lambda @f$ is the dimensionless variable parameterizing
   the sub-grid grounding line position, ranging from 0 to 1.

   Now, substituting @f$ b(\lambda) @f$ and @f$ H(\lambda) @f$ into
   the floatation criterion

   @f{align*}{
   \mu\cdot H(\lambda) &= z_{\text{sea level}} - b(\lambda),\\
   \mu &= \rho_{\text{ice}} / \rho_{\text{sea water}}.
   @f}

   and solving for @f$ \lambda @f$, we get

   @f{align*}{
   \lambda_{g} &= \frac{\alpha}{\alpha - \beta}\\
   & \text{where} \\
   \alpha &= \mu\cdot H_{i} + b_{i} - z_{\text{sea level}}\\
   \beta &= \mu\cdot H_{i+1} + b_{i+1} - z_{\text{sea level}}.
   @f}

   Note that [@ref Gladstoneetal2010] describe a parameterization of
   the grounding line position within a cell defined as the interval
   from the grid point `(i)` to the grid point `(i+1)`, with the ice
   thickness @f$ H @f$ and the bed elevation @f$ b @f$ defined *at
   grid points* (boundaries of a 1D cell).

   Here we compute a grounded fraction of the **cell centered at the
   grid point**. Grid-point-centered cells and cells with grid points
   at their corners are shifted by 0.5 of a cell width relative to
   each other.

   This explains if-clauses like this one:
   ~~~ c++
   if (lambda_g < 0.5)
     gl_mask_gr_x += (lambda_g - 0.5);
   ~~~

   they convert the sub-grid grounding line position into the form
   used by PISM.

   FIXME: sometimes alpha<0 (slightly below flotation) even though the
   mask says grounded
 */
PetscErrorCode IceModel::update_floatation_mask() {
  PetscErrorCode ierr;

  MaskQuery mask(vMask);

  double
    ice_density   = config.get("ice_density"),
    ocean_density = config.get("sea_water_density"),
    mu            = ice_density / ocean_density,
    sea_level     = 0.0;

  assert(ocean != NULL);
  ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);

  ierr = gl_mask.set(0.0);   CHKERRQ(ierr);
  ierr = gl_mask_x.set(0.0); CHKERRQ(ierr);
  ierr = gl_mask_y.set(0.0); CHKERRQ(ierr);

  IceModelVec::AccessList list;
  list.add(ice_thickness);
  list.add(bed_topography);
  list.add(vMask);
  list.add(gl_mask);
  list.add(gl_mask_x);
  list.add(gl_mask_y);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double
      alpha        = 0.0,
      beta         = 0.0,
      lambda_g     = 0.0,
      gl_mask_gr_x = 1.0,
      gl_mask_gr_y = 1.0,
      gl_mask_fl_x = 1.0,
      gl_mask_fl_y = 1.0;

    // grounded part
    if (mask.grounded(i, j)) {
      alpha = mu*ice_thickness(i, j) + bed_topography(i, j) - sea_level;

      if (mask.ocean(i + 1, j)) {

        beta = mu*ice_thickness(i + 1, j) + bed_topography(i + 1, j) - sea_level;

        assert(alpha - beta != 0.0);
        lambda_g = alpha / (alpha - beta);
        lambda_g = std::max(0.0, std::min(lambda_g, 1.0));

        if (lambda_g < 0.5)
          gl_mask_gr_x += (lambda_g - 0.5);

      } else if (mask.ocean(i - 1, j)) {

        beta = mu*ice_thickness(i - 1, j) + bed_topography(i - 1, j) - sea_level;

        assert(alpha - beta != 0.0);
        lambda_g = alpha / (alpha - beta);
        lambda_g = std::max(0.0, std::min(lambda_g, 1.0));

        if (lambda_g < 0.5)
          gl_mask_gr_x += (lambda_g - 0.5);

      } else if (mask.ocean(i, j + 1)) {

        beta = mu*ice_thickness(i, j + 1) + bed_topography(i, j + 1) - sea_level;

        assert(alpha - beta != 0.0);
        lambda_g = alpha / (alpha - beta);
        lambda_g = std::max(0.0, std::min(lambda_g, 1.0));

        if (lambda_g < 0.5)
          gl_mask_gr_y += (lambda_g - 0.5);

      } else if (mask.ocean(i, j - 1)) {

        beta = mu*ice_thickness(i, j - 1) + bed_topography(i, j - 1) - sea_level;

        assert(alpha - beta != 0.0);
        lambda_g = alpha / (alpha - beta);
        lambda_g = std::max(0.0, std::min(lambda_g, 1.0));

        if (lambda_g < 0.5)
          gl_mask_gr_y += (lambda_g - 0.5);

      }

      gl_mask_x(i,j) = gl_mask_gr_x;
      gl_mask_y(i,j) = gl_mask_gr_y;
      gl_mask(i,j)   = gl_mask_gr_x * gl_mask_gr_y;
    }

    // floating part
    if (mask.floating_ice(i, j)) {
      beta = mu*ice_thickness(i, j) + bed_topography(i, j) - sea_level;

      if (mask.grounded(i - 1, j)) {

        alpha = mu*ice_thickness(i - 1, j) + bed_topography(i - 1, j) - sea_level;

        assert(alpha - beta != 0.0);
        lambda_g = alpha / (alpha - beta);
        lambda_g = std::max(0.0, std::min(lambda_g, 1.0));

        if (lambda_g >= 0.5)
          gl_mask_fl_x -= (lambda_g - 0.5);

      } else if (mask.grounded(i + 1, j)) {

        alpha = mu*ice_thickness(i + 1, j) + bed_topography(i + 1, j) - sea_level;

        assert(alpha - beta != 0.0);
        lambda_g = alpha / (alpha - beta);
        lambda_g = std::max(0.0, std::min(lambda_g, 1.0));

        if (lambda_g >= 0.5)
          gl_mask_fl_x -= (lambda_g - 0.5);

      } else if (mask.grounded(i, j - 1)) {

        alpha = mu*ice_thickness(i, j - 1) + bed_topography(i, j - 1) - sea_level;

        assert(alpha - beta != 0.0);
        lambda_g = alpha / (alpha - beta);
        lambda_g = std::max(0.0, std::min(lambda_g, 1.0));

        if (lambda_g >= 0.5)
          gl_mask_fl_y -= (lambda_g - 0.5);

      } else if (mask.grounded(i, j + 1)) {

        alpha = mu*ice_thickness(i, j + 1) + bed_topography(i, j + 1) - sea_level;

        assert(alpha - beta != 0.0);
        lambda_g = alpha / (alpha - beta);
        lambda_g = std::max(0.0, std::min(lambda_g, 1.0));

        if (lambda_g >= 0.5)
          gl_mask_fl_y -= (lambda_g - 0.5);

      }

      gl_mask_x(i,j) = 1.0 - gl_mask_fl_x;
      gl_mask_y(i,j) = 1.0 - gl_mask_fl_y;
      gl_mask(i,j)   = 1.0 - gl_mask_fl_x * gl_mask_fl_y;
    }
  }

  return 0;
}

} // end of namespace pism
