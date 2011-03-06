// Copyright (C) 2004-2011 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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
#include <petscda.h>
#include "iceModel.hh"

//! Update the surface elevation and the flow-type mask when the geometry has changed.
/*!
This procedure should be called whenever necessary to maintain consistency of geometry.

For instance, it should be called when either ice thickness or bed elevation change. 
In particular we always want \f$h = H + b\f$ to apply at grounded points, and, on the
other hand, we want the mask to reflect that the ice is floating if the flotation 
criterion applies at a point.

There is one difficult case.  When a point was floating and becomes grounded we generally
do not know whether to mark it as \c MASK_SHEET so that the SIA applies or \c MASK_DRAGGING
so that the SSA applies.  For now there is a vote-by-neighbors scheme (among the grounded 
neighbors).  When the \c MASK_DRAGGING points have plastic till bases this is not an issue.
 */
PetscErrorCode IceModel::updateSurfaceElevationAndMask() {
  PetscErrorCode ierr;

  ierr = update_mask(); CHKERRQ(ierr);
  ierr = update_surface_elevation(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::update_mask() {
  PetscErrorCode ierr;

  if (ocean == PETSC_NULL) {  SETERRQ(1,"PISM ERROR: ocean == PETSC_NULL");  }
  PetscReal currentSeaLevel;
  ierr = ocean->sea_level_elevation(grid.year, dt / secpera, currentSeaLevel); CHKERRQ(ierr);

  bool use_ssa_when_grounded = config.get_flag("use_ssa_when_grounded"),
    use_ssa_velocity = config.get_flag("use_ssa_velocity");

  double ocean_rho = config.get("sea_water_density");

  ierr =    vH.begin_access();    CHKERRQ(ierr);
  ierr =  vbed.begin_access();  CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);

  PetscInt GHOSTS = 2;
  for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {

      const PetscScalar hgrounded = vbed(i,j) + vH(i,j), // FIXME task #7297
	hfloating = currentSeaLevel + (1.0 - ice->rho/ocean_rho) * vH(i,j);

      // points marked as "ocean at time zero" are not updated
      if (vMask.value(i,j) == MASK_OCEAN_AT_TIME_0)
	continue;

      if (vMask.is_floating(i,j)) { // floating

	if (hgrounded > hfloating+1.0) { // flotation criterion says it is grounded
	  if (use_ssa_velocity) {
	    if (use_ssa_when_grounded) {
	      vMask(i,j) = MASK_DRAGGING_SHEET;
	    } else {
	      vMask(i,j) = MASK_SHEET;
	    }
	  } else {
	    // we do not have any ice handled by SSA, so it must be SHEET
	    vMask(i,j) = MASK_SHEET;
	  }
	}

      } else {   // grounded

	// apply the flotation criterion:
	if (hgrounded > hfloating-1.0) { // flotation criterion says it is grounded

	  // we are using SSA-as-a-sliding-law, so grounded points become DRAGGING
	  if (use_ssa_velocity && use_ssa_when_grounded)
	    vMask(i,j) = MASK_DRAGGING_SHEET;

	} else {
	  vMask(i,j) = MASK_FLOATING;
	}

      }
    } // inner for loop (j)
  } // outer for loop (i)

  ierr =         vH.end_access(); CHKERRQ(ierr);
  ierr =       vbed.end_access(); CHKERRQ(ierr);
  ierr =      vMask.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::update_surface_elevation() {
  PetscErrorCode ierr;

  if (ocean == PETSC_NULL) {  SETERRQ(1,"PISM ERROR: ocean == PETSC_NULL");  }
  PetscReal currentSeaLevel;
  ierr = ocean->sea_level_elevation(grid.year, dt / secpera, currentSeaLevel); CHKERRQ(ierr);

  bool is_dry_simulation = config.get_flag("is_dry_simulation");

  double ocean_rho = config.get("sea_water_density");

  ierr =    vh.begin_access();    CHKERRQ(ierr);
  ierr =    vH.begin_access();    CHKERRQ(ierr);
  ierr =  vbed.begin_access();  CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);

  PetscInt GHOSTS = 2;
  for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
      // take this opportunity to check that vH(i,j) >= 0
      if (vH(i,j) < 0) {
        SETERRQ2(1,"Thickness negative at point i=%d, j=%d",i,j);
      }

      const PetscScalar hgrounded = vbed(i,j) + vH(i,j), // FIXME task #7297
	hfloating = currentSeaLevel + (1.0 - ice->rho/ocean_rho) * vH(i,j);

      if (is_dry_simulation) {
        // Don't update mask; potentially one would want to do SSA dragging ice
        //   shelf in dry case and/or ignore mean sea level elevation.
        vh(i,j) = hgrounded;
	continue;		// go to the next grid point
      }

      if (vMask.value(i,j) == MASK_OCEAN_AT_TIME_0) {
        // mask takes priority over bed in this case (note sea level may change).
        // Example Greenland case: if mask say Ellesmere is OCEAN0,
        //   then never want ice on Ellesmere.
        // If mask says OCEAN0 then don't change the mask and also don't change
        // the thickness; massContExplicitStep() is in charge of that.
        // Almost always the next line is equivalent to vh(i,j) = 0.
        vh(i,j) = hfloating;  // ignore bed and treat it like deep ocean
	continue;	      // go to the next grid point
      }

      if (vMask.is_floating(i,j)) {
	vh(i,j) = hfloating; // actually floating so update h
      } else { 
	vh(i,j) = hgrounded; // actually grounded so update h
      }
    }
        
  }

  ierr =         vh.end_access(); CHKERRQ(ierr);
  ierr =         vH.end_access(); CHKERRQ(ierr);
  ierr =       vbed.end_access(); CHKERRQ(ierr);
  ierr =      vMask.end_access(); CHKERRQ(ierr);

  return 0;
}

//! Update the thickness from the diffusive flux, also additional horizontal velocity, and the surface and basal mass balance rates.
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
conservation of mass to update the ice thickness.

The PISMSurfaceModel pointed to by IceModel::surface provides \f$M\f$.  The
PISMOceanModel pointed to by IceModel::ocean provides \f$S\f$ at locations below
floating ice (ice shelves).

Even if we regard the map-plane flux as defined by the formula
\f$\mathbf{q} = \bar{\mathbf{U}} H\f$, the flow of ice is at least somewhat
diffusive in almost all cases, and in simplified models (%e.g.~the SIA model) it
is exactly true that \f$\mathbf{q} = - D \nabla h\f$.  In the current method the
flux is split into the part from the diffusive non-sliding SIA model
and a part which is a less-diffusive, possibly membrane-stress-dominated
2D advective velocity.  That is, the flux is split
this way, which is common in the literature:
  \f[ \mathbf{q} = - D \nabla h + \mathbf{U}_b H.\f]
Here \f$D\f$ is the (positive, scalar) effective diffusivity of the non-sliding
SIA and \f$\mathbf{U}_b\f$ is the less-diffusive velocity from the membrane
stress balance.  We may interpret \f$\mathbf{U}_b\f$ as a basal sliding velocity
in the hybrid or in classical SIA sliding schemes (though the latter are not
recommended).  

The main ice-dynamical inputs to this method are identified in these lines,
which get outputs from PISMStressBalance *stress_balance:
\code
  IceModelVec2Stag *Qdiff;
  stress_balance->get_diffusive_flux(Qdiff);
  IceModelVec2V *vel_advective;
  stress_balance->get_advective_2d_velocity(vel_advective);
\endcode
The diffusive flux \f$-D\nabla h\f$ is thus stored in \c IceModelVec2Stag
\c *Qdiff while the less-diffusive velocity \f$\mathbf{U}_b\f$ is stored in
\c IceModelVec2V \c *vel_advective.

The methods used here are first-order and explicit in time.  The derivatives in 
\f$\nabla \cdot (D \nabla h)\f$ are computed by centered finite difference
methods.  The diffusive flux \c Qdiff is already stored on the staggered grid.
The time-stepping for this part of the explicit scheme is controlled by equation
(25) in [\ref BBL], so that \f$\Delta t \sim \Delta x^2 / \max D\f$; see
also [\ref MortonMayers].

The divergence of the flux from the less-diffusive velocity is computed by
the PIK upwinding technique [\ref Albrechtetal2011, \ref Winkelmannetal2010TCD].
The CFL condition for this advection scheme is checked; see 
computeMax2DSlidingSpeed() and determineTimeStep().

Checks are made which can generate zero thickness according to minimal calving
relations, specifically the mechanisms turned-on by options \c -ocean_kill and
\c -float_kill.

The rate of thickness change \f$\partial H/\partial t\f$ is computed and saved,
as is the rate of volume loss or gain.
 */
PetscErrorCode IceModel::massContExplicitStep() {
  PetscErrorCode ierr;
  PetscScalar my_nonneg_rule_flux = 0, my_ocean_kill_flux = 0, my_float_kill_flux = 0;

  const PetscScalar   dx = grid.dx, dy = grid.dy;
  bool do_ocean_kill = config.get_flag("ocean_kill"),
    floating_ice_killed = config.get_flag("floating_ice_killed"),
    include_bmr_in_continuity = config.get_flag("include_bmr_in_continuity");

  if (surface != NULL) {
    ierr = surface->ice_surface_mass_flux(grid.year, dt / secpera, acab); CHKERRQ(ierr);
  } else { SETERRQ(1,"PISM ERROR: surface == NULL"); }

  if (ocean != NULL) {
    ierr = ocean->shelf_base_mass_flux(grid.year, dt / secpera, shelfbmassflux); CHKERRQ(ierr);
  } else { SETERRQ(2,"PISM ERROR: ocean == NULL"); }

  IceModelVec2S vHnew = vWork2d[0];
  ierr = vH.copy_to(vHnew); CHKERRQ(ierr);

  IceModelVec2Stag *Qdiff;
  ierr = stress_balance->get_diffusive_flux(Qdiff); CHKERRQ(ierr);

  IceModelVec2V *vel_advective;
  ierr = stress_balance->get_advective_2d_velocity(vel_advective); CHKERRQ(ierr);
  IceModelVec2V vel = *vel_advective; // just an alias

  PetscScalar **bmr_gnded;
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vbmr.get_array(bmr_gnded); CHKERRQ(ierr);
  ierr = Qdiff->begin_access(); CHKERRQ(ierr);
  ierr = vel.begin_access(); CHKERRQ(ierr);
  ierr = acab.begin_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access();  CHKERRQ(ierr);
  ierr = vHnew.begin_access(); CHKERRQ(ierr);
  
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      PetscScalar divQ = 0.0;

      if (!vMask.is_floating(i,j)) {
        // staggered grid Div(Q) for diffusive non-sliding SIA deformation part:
        //    Qdiff = - D grad h
        divQ =   ((*Qdiff)(i,j,0) - (*Qdiff)(i-1,j,0)) / dx
               + ((*Qdiff)(i,j,1) - (*Qdiff)(i,j-1,1)) / dy;
      }

      // membrane stress (and/or basal sliding) part: upwind by staggered grid
      // PIK method;  this is   \nabla \cdot [(u,v) H]
      const PetscScalar // compute (i,j)-centered "face" velocity components by average
          velE = 0.5 * (vel(i,j).u + vel(i+1,j).u),
          velW = 0.5 * (vel(i-1,j).u + vel(i,j).u),
          velN = 0.5 * (vel(i,j).v + vel(i,j+1).v),
          velS = 0.5 * (vel(i,j-1).v + vel(i,j).v);
      divQ += (  velE * (velE > 0 ? vH(i,j) : vH(i+1,j))
               - velW * (velW > 0 ? vH(i-1,j) : vH(i,j)) ) / dx;
      divQ += (  velN * (velN > 0 ? vH(i,j) : vH(i,j+1))
               - velS * (velS > 0 ? vH(i,j-1) : vH(i,j)) ) / dy;

      vHnew(i,j) += (acab(i,j) - divQ) * dt; // include M

      if (include_bmr_in_continuity) { // include S
        if (vMask.is_floating(i,j)) {
	  vHnew(i,j) -= shelfbmassflux(i,j) * dt;
        } else {
          vHnew(i,j) -= bmr_gnded[i][j] * dt;
        }
      }

      // apply free boundary rule: negative thickness becomes zero
      if (vHnew(i,j) < 0) {
        my_nonneg_rule_flux += (-vHnew(i,j));
        vHnew(i,j) = 0.0;
      }

      // the following conditionals, both -ocean_kill and -float_kill, are also applied in 
      //   IceModel::computeMax2DSlidingSpeed() when determining CFL
      
      // force zero thickness at points which were originally ocean (if "-ocean_kill");
      //   this is calving at original calving front location
      if ( do_ocean_kill && (vMask.value(i,j) == MASK_OCEAN_AT_TIME_0) ) {
        my_ocean_kill_flux -= vHnew(i,j);
        vHnew(i,j) = 0.0;
      }

      // force zero thickness at points which are floating (if "-float_kill");
      //   this is calving at grounding line
      if ( floating_ice_killed && vMask.is_floating(i,j) ) {
        my_float_kill_flux -= vHnew(i,j);
        vHnew(i,j) = 0.0;
      }

    } // end of the inner for loop
  } // end of the outer for loop

  ierr = vbmr.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = Qdiff->end_access(); CHKERRQ(ierr);
  ierr = vel.end_access(); CHKERRQ(ierr);
  ierr = acab.end_access(); CHKERRQ(ierr);
  ierr = shelfbmassflux.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vHnew.end_access(); CHKERRQ(ierr);

  {
    ierr = PetscGlobalSum(&my_nonneg_rule_flux, &nonneg_rule_flux, grid.com); CHKERRQ(ierr);
    ierr = PetscGlobalSum(&my_ocean_kill_flux,  &ocean_kill_flux,  grid.com); CHKERRQ(ierr);
    ierr = PetscGlobalSum(&my_float_kill_flux,  &float_kill_flux,  grid.com); CHKERRQ(ierr);

    // FIXME: use corrected cell areas (when available)
    PetscScalar ice_density = config.get("ice_density"),
      factor = ice_density * (dx * dy) / dt;
    nonneg_rule_flux *= factor;
    ocean_kill_flux  *= factor;
    float_kill_flux  *= factor;
  }

  // compute dH/dt (thickening rate) for viewing and for saving at end; only diagnostic
  ierr = vHnew.add(-1.0, vH, vdHdt); CHKERRQ(ierr); // vdHdt = vHnew - vH
  ierr = vdHdt.scale(1.0/dt); CHKERRQ(ierr);	    // vdHdt = vdHdt / dt

  // d(volume)/dt
  {
    PetscScalar dvol=0.0;
  
    ierr = vdHdt.begin_access(); CHKERRQ(ierr);
    ierr = cell_area.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        dvol += vdHdt(i,j) * cell_area(i,j);
      }
    }  
    ierr = cell_area.end_access(); CHKERRQ(ierr);
    ierr = vdHdt.end_access(); CHKERRQ(ierr);

    ierr = PetscGlobalSum(&dvol, &dvoldt, grid.com); CHKERRQ(ierr);
  }

  // average value of dH/dt; 
  PetscScalar ice_area;
  ierr = compute_ice_area(ice_area); CHKERRQ(ierr);

  ierr = vdHdt.sum(gdHdtav); CHKERRQ(ierr);
  gdHdtav = gdHdtav / ice_area; // m/s
  
  // finally copy vHnew into vH and communicate ghosted values
  ierr = vHnew.beginGhostComm(vH); CHKERRQ(ierr);
  ierr = vHnew.endGhostComm(vH); CHKERRQ(ierr);

  // Check if the ice thickness exceeded the height of the computational box
  // and extend the grid if necessary:
  ierr = check_maximum_thickness(); CHKERRQ(ierr);

  return 0;
}

