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

Also calls the (PIK) routines which remove icebergs, to avoid stress balance
solver problems associated to not-attached-to-grounded ice.
 */
PetscErrorCode IceModel::updateSurfaceElevationAndMask() {
  PetscErrorCode ierr;

  ierr = update_mask(); CHKERRQ(ierr);
  ierr = update_surface_elevation(); CHKERRQ(ierr);

  if (config.get_flag("kill_icebergs") == true) {

    ierr = findIceBergCandidates(); CHKERRQ(ierr);
	ierr = identifyNotAnIceBerg(); CHKERRQ(ierr);
	ierr = killIceBergs(); CHKERRQ(ierr);	
	if (config.get_flag("do_eigen_calving") || config.get_flag("do_thickness_calving") ) {
    ierr = killEasyIceBergs(); CHKERRQ(ierr);
	}
  }


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

      bool is_floating = hfloating > hgrounded + 1.0,
        // note: the following implies that ice-free cells with bed evelation
        // exactly at sea level are considered grounded
        is_grounded = ! is_floating,
        ice_free = vH(i,j) < 0.01;

      // points marked as "ocean at time zero" are not updated
      if (vMask.value(i,j) == MASK_OCEAN_AT_TIME_0)
        continue;

      if (is_floating) {
        if (ice_free) {
          vMask(i,j) = MASK_ICE_FREE_OCEAN;
          // added just for clarity, needs to be tested for interference in run
        } else {					
          vMask(i,j) = MASK_FLOATING; // to enable for floating front propagation
        }
      }

      if (is_grounded) {

        if (ice_free) {
          vMask(i,j) = MASK_ICE_FREE_BEDROCK;
        } else if (use_ssa_velocity && use_ssa_when_grounded)
          vMask(i,j) = MASK_GROUNDED;
        else
          vMask(i,j) = MASK_GROUNDED; // for historical reasons

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
methods.  The diffusive flux \c Qdiff is already stored on the staggered grid
and it is differenced in a centered way here.  The time-stepping for this part
of the explicit scheme is controlled by equation (25) in [\ref BBL], so that
\f$\Delta t \sim \Delta x^2 / \max D\f$; see also [\ref MortonMayers].

The divergence of the flux from velocity \f$\mathbf{U}_b\f$ is computed by
the PIK upwinding technique [equation (25) in \ref Winkelmannetal2010TCD; see
also \ref Albrechtetal2011]. The CFL condition for this advection scheme is checked; see 
computeMax2DSlidingSpeed() and determineTimeStep().  This method implements the
direct-superposition (PIK) hybrid which adds the SSA velocity, as a basal
sliding velocity, to the SIA velocity; [see equation (15) in \ref Winkelmannetal2010TCD].
The hybrid described by equations (21) and (22) in \ref BBL is no longer used
for this purpose.

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


// the update of H is altered along the ice front boundary in terms of partially filled grid cells as in PISM-PIK
PetscErrorCode IceModel::massContExplicitStepPartGrids() {
  PetscErrorCode ierr;

  ierr = verbPrintf(4,grid.com, "massContExplicitStepPartGrids() is called\n"); CHKERRQ(ierr);

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
  ierr = vHref.begin_access(); CHKERRQ(ierr); 

  const bool do_redist = config.get_flag("part_redist");
   if (do_redist) {
		ierr = vHresidual.begin_access(); CHKERRQ(ierr);
		ierr = vHresidual.set(0.0); CHKERRQ(ierr); //mass loss if max_loopcount for redistribution was not sufficient
	}

  
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      // mask values at the current cell and its four immediate neighbors
      // (east, west, north, south)
      int Mo = vMask.value(i,j),
        Me = vMask.value(i+1,j),
        Mw = vMask.value(i-1,j),
        Mn = vMask.value(i,j+1),
        Ms = vMask.value(i,j-1);

      // advective velocities at cell interfaces (sides)
      PetscScalar		
        velE = 0.0,
        velW = 0.0,
        velN = 0.0,
        velS = 0.0;

      //just to make sure, that no velocities on the ice free ocean are used
      //case1: [i][j] in the middle of ice or bedrock: default scheme
      if (Mo <= MASK_FLOATING && 
          Me <= MASK_FLOATING &&
          Mw <= MASK_FLOATING &&
          Mn <= MASK_FLOATING &&
          Ms <= MASK_FLOATING) {

        // compute (i,j)-centered "face" velocity components by average
        velE = 0.5 * (vel(i,j).u + vel(i+1,j).u);
        velW = 0.5 * (vel(i-1,j).u + vel(i,j).u);
        velN = 0.5 * (vel(i,j).v + vel(i,j+1).v);
        velS = 0.5 * (vel(i,j-1).v + vel(i,j).v);


        //case2: [i][j] on floating or grounded ice, but next to a ice-free ocean grid cell
      } else if (Mo <= MASK_FLOATING && 
                 (Me > MASK_FLOATING ||
                  Mw > MASK_FLOATING ||
                  Mn > MASK_FLOATING ||
                  Ms > MASK_FLOATING)) {

        // velocities on ice-free ocean may not be valid for averaging staggered velocity
        velE = (Me > MASK_FLOATING ? vel(i,j).u : 0.5 * (vel(i,j).u + vel(i+1,j).u));
        velW = (Mw > MASK_FLOATING ? vel(i,j).u : 0.5 * (vel(i-1,j).u + vel(i,j).u));
        velN = (Mn > MASK_FLOATING ? vel(i,j).v : 0.5 * (vel(i,j).v + vel(i,j+1).v));
        velS = (Ms > MASK_FLOATING ? vel(i,j).v : 0.5 * (vel(i,j-1).v + vel(i,j).v));

        //case3: [i][j] on ice-free ocean (or partially filled), but next to ice grid cell
      } else if (Mo > MASK_FLOATING && 
                 (Me <= MASK_FLOATING ||
                  Mw <= MASK_FLOATING ||
                  Mn <= MASK_FLOATING ||
                  Ms <= MASK_FLOATING)) {

        // velocities on ice-free ocean may not be valid for averaging staggered velocity
        velE = (Me <= MASK_FLOATING ? vel(i+1,j).u : 0.0);
        velW = (Mw <= MASK_FLOATING ? vel(i-1,j).u : 0.0);
        velN = (Mn <= MASK_FLOATING ? vel(i,j+1).v : 0.0);
        velS = (Ms <= MASK_FLOATING ? vel(i,j-1).v : 0.0);

        //case4: [i][j] on ice-free ocean, and no ice neighbors, and else
      } else {		
        velE = 0.0;
        velW = 0.0;
        velN = 0.0;
        velS = 0.0;
      }
	
      // here divQ is calculated
      PetscScalar divQ = 0.0;
      // staggered grid Div(Q) for diffusive non-sliding SIA deformation part:
      //    Qdiff = - D grad h
      if (vMask.is_grounded(i,j)) {
        divQ = ((*Qdiff)(i,j,0) - (*Qdiff)(i-1,j,0)) / dx
          + ((*Qdiff)(i,j,1) - (*Qdiff)(i,j-1,1)) / dy;
      }

      // membrane stress (and/or basal sliding) part: upwind by staggered grid
      // PIK method;  this is   \nabla \cdot [(u,v) H]
      divQ += (  velE * (velE > 0 ? vH(i,j) : vH(i+1,j))
                 - velW * (velW > 0 ? vH(i-1,j) : vH(i,j)) ) / dx;
      divQ += (  velN * (velN > 0 ? vH(i,j) : vH(i,j+1))
                 - velS * (velS > 0 ? vH(i,j-1) : vH(i,j)) ) / dy;


        // ice-free (or partially-filled) cells adjacent to "full" floating cells
      if ((Mo > MASK_FLOATING) &&
          (Me == MASK_FLOATING ||
           Mw == MASK_FLOATING ||
           Mn == MASK_FLOATING ||
           Ms == MASK_FLOATING)) { //does in this form not account for grounded tributaries: thin ice shelves may evolve from grounded tongue

          PetscScalar Hav = 0.0;
		///*	
		  PetscInt countIceNeighbors=0; // counting existing ice neighbors
		
		  // mean ice thickness over adjacent floating ice shelf neighbors
          if (Me == MASK_FLOATING) { Hav+=vH(i+1,j); countIceNeighbors+=1;} 
          if (Mw == MASK_FLOATING) { Hav+=vH(i-1,j); countIceNeighbors+=1;}
          if (Mn == MASK_FLOATING) { Hav+=vH(i,j+1); countIceNeighbors+=1;}
          if (Ms == MASK_FLOATING) { Hav+=vH(i,j-1); countIceNeighbors+=1;}

 		  if (countIceNeighbors>0){ 
	      Hav=Hav/countIceNeighbors;
   	   	  if (config.get_flag("part_redist") == true) {	
   		   	const PetscReal  mslope = 2.4511e-18*grid.dx/(300*600/secpera);
   		  //    for declining front C/Q0 according to analytical flowline profile in vandeveen with v0=300m/yr and H0=600m	    
   		    Hav-=0.8*mslope*pow(Hav,5); //reduces the guess at the front
   		  }
 		  } else {
   		  ierr = verbPrintf(4, grid.com,"!!! PISM_WARNING: no ice shelf neighbors at %d,%d\n",i,j); CHKERRQ(ierr);}
		//*/
		/*
		  // alternative: flux-weighted average over floating ice-shelf neighbors
		  if (Me == MASK_FLOATING && vel(i+1,j).u < 0.0) { Hav-= vel(i+1,j).u * vH(i+1,j);} 
          if (Mw == MASK_FLOATING && vel(i-1,j).u > 0.0) { Hav+= vel(i-1,j).u * vH(i-1,j);}
          if (Mn == MASK_FLOATING && vel(i,j+1).v < 0.0) { Hav-= vel(i,j+1).v * vH(i,j+1);}
          if (Ms == MASK_FLOATING && vel(i,j-1).v > 0.0) { Hav+= vel(i,j-1).v * vH(i-1,j);}
		
		  // velocity magnitude
		  //PetscScalar velRoot = sqrt((vel(i+1,j).u < 0 ? PetscSqr(vel(i+1,j).u) : 0.0) + 
		  //						 (vel(i-1,j).u > 0 ? PetscSqr(vel(i-1,j).u) : 0.0) +
	      //						 (vel(i,j+1).v < 0 ? PetscSqr(vel(i,j+1).v) : 0.0) +
		  //						 (vel(i,j-1).v > 0 ? PetscSqr(vel(i,j-1).v) : 0.0));
		  // or just the sum:							
		  PetscScalar velRoot =     ((vel(i+1,j).u < 0 ? vel(i+1,j).u : 0.0) + 
									( vel(i-1,j).u > 0 ? vel(i-1,j).u : 0.0) +
									( vel(i,j+1).v < 0 ? vel(i,j+1).v : 0.0) +
									( vel(i,j-1).v > 0 ? vel(i,j-1).v : 0.0));
			
		 PetscScalar minHav=50.0;	
		 if (velRoot>0 && Hav>=minHav) { 
            Hav=Hav/velRoot;
          } else { //can occur when flux comes from grounded cell, but next to floating cell
			Hav=minHav;
			//ierr = verbPrintf(3, grid.com,"!!! PISM_WARNING: no flux into partially filled grid cell at %d,%d\n",i,j); CHKERRQ(ierr);
		  }
		*/
          vHref(i,j)  -= divQ * dt;	
          vHnew(i,j) = 0.0; // redundant

          // applying M and S also to partially filled cells	
          // to calculate the S-M contribution with respect to the final coverage ratio, let's assume	
		  // $ vHref_new = vHref_old + (acab-shelfbmassflux) * dt * coverageRatio $, which equals
		  // $ vHref_new = vHref_old + (acab-shelfbmassflux) * dt * vHref_new / Hav $.
		  // Hence we get $ vHref_new = vHref_old / (1.0 - (acab-shelfbmassflux) * dt * Hav))
		
		  PetscScalar denominator =  1.0 - dt / Hav * (acab(i,j) - (include_bmr_in_continuity ? shelfbmassflux(i,j) : 0.0));
          vHref(i,j) = vHref(i,j) / denominator;
			

          PetscScalar coverageRatio = vHref(i,j)/Hav;
          if (coverageRatio>1.0) { // partially filled grid cell is considered to be full
			if (do_redist) {
				vHresidual(i,j)=vHref(i,j)-Hav; //residual ice thickness
				//ierr = verbPrintf(3,grid.com,"!!! Hresidual=%.2f for Href=%.2f, Hav=%.2f and acab-melt-factor=%.3f at %d,%d \n",vHresidual(i,j),vHref(i,j),Hav,1/denominator,i,j); CHKERRQ(ierr);
				//ierr = verbPrintf(3,grid.com,"!!! Hresidual=%.5f for Href=%.5f, Hav=%.5f, velRoot=%.5f and acab-melt-factor=%.5f at %d,%d \n",vHresidual(i,j),vHref(i,j),Hav,velRoot*secpera,1/denominator,i,j); CHKERRQ(ierr);
			}
            vHnew(i,j) = Hav; //gets a "real" ice thickness
            vHref(i,j)=0.0;
          }	
		
      } else if ((Mo <= MASK_FLOATING) || // grounded/floating default case
          ((Mo > MASK_FLOATING) &&
           (Me < MASK_FLOATING ||
            Mw < MASK_FLOATING ||
            Mn < MASK_FLOATING ||
            Ms < MASK_FLOATING))) { // ice-free ocean, adjacent to a grounded cell: FIXME
        
        vHnew(i,j) += (acab(i,j) - divQ) * dt; // include M
		
        if (include_bmr_in_continuity) { // include S
          if (vMask.is_floating(i,j)) {
            vHnew(i,j) -= shelfbmassflux(i,j) * dt;
          } else {
            vHnew(i,j) -= bmr_gnded[i][j] * dt;
          }
        }
      }

	else { // basically ice-free, not adjacent to a "full" cell, no matter what kind of
		vHnew(i,j)=0.0;
	}

      // FIXME: take care of ice-free cells adjacent to grounded ice


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
  ierr = vHref.end_access(); CHKERRQ(ierr);
  if (config.get_flag("part_redist") == true) {
	ierr = vHresidual.end_access(); CHKERRQ(ierr);
  }

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
  } //FIXME: Reporting not yet adjusted to usage of partially filled grid cell scheme


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

  // These are new routines adopted from PISM-PIK. The Place and order is not clear yet!
  // There is no reporting of single ice fluxes yet in comparison to total ice thickness change.

	if (config.get_flag("part_redist") == true) {
		// distribute residual ice mass, FIXME: Reporting!
		ierr = calculateRedistResiduals(); CHKERRQ(ierr); //while loop?
		PetscInt loopcount=0;
		while (repeatRedist==PETSC_TRUE && loopcount<3) {
			ierr = calculateRedistResiduals(); CHKERRQ(ierr);
			loopcount+=1;
			ierr = verbPrintf(4,grid.com, "distribution loopcount = %d\n",loopcount); CHKERRQ(ierr);
		}
	}



	// maybe calving should be applied before the redistribution part?
	if (config.get_flag("do_eigen_calving") == true) {
		ierr = stress_balance->get_principle_strain_rates(
		                         vPrinStrain1,vPrinStrain2); CHKERRQ(ierr);
		ierr = eigenCalving(); CHKERRQ(ierr);
	}

	if (config.get_flag("do_thickness_calving")==true) { 
		if (config.get_flag("part_grid")==true) { 
			ierr = calvingAtThickness(); CHKERRQ(ierr);
		} else {
			ierr = verbPrintf(2,grid.com, "PISM-WARNING: calving at certian terminal ice thickness without application of partially filled grid cell scheme may lead to non-moving ice shef front!\n"); CHKERRQ(ierr);
		}
	}

  return 0;
}

// This routine takes care of redistributing carry over ice mass when using -redist option
PetscErrorCode IceModel::calculateRedistResiduals() {
  	PetscErrorCode ierr;
	ierr = verbPrintf(4,grid.com, "calculateRedistResiduals() is called\n"); CHKERRQ(ierr);

  	IceModelVec2S vHnew = vWork2d[0];
  	ierr = vH.copy_to(vHnew); CHKERRQ(ierr);
  	ierr = vHnew.begin_access(); CHKERRQ(ierr);
  	ierr = vH.begin_access(); CHKERRQ(ierr);

  	ierr = vHref.begin_access(); CHKERRQ(ierr); 
	ierr = vbed.begin_access(); CHKERRQ(ierr);
	
  	IceModelVec2S vHresidualnew = vWork2d[1];
  	ierr = vHresidual.copy_to(vHresidualnew); CHKERRQ(ierr);
	ierr = vHresidual.begin_access(); CHKERRQ(ierr);
	ierr = vHresidualnew.begin_access(); CHKERRQ(ierr);

	
	if (ocean == PETSC_NULL) {  SETERRQ(1,"PISM ERROR: ocean == PETSC_NULL");  }
    PetscReal currentSeaLevel=0.0; //FIXME
    ierr = ocean->sea_level_elevation(grid.year, dt / secpera, currentSeaLevel); CHKERRQ(ierr);
	
	PetscScalar minHRedist = 0.0; // to avoid the propagation of thin ice shelf tongues
	
	for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    	for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
			// first step: distributing residual ice masses
			if (vHresidual(i,j)>0.0) {
	
	          PetscInt countEmptyNeighbors=0; // counting empty/partially filled neighbors
	          PetscTruth exEast=PETSC_FALSE,exWest=PETSC_FALSE,exNorth=PETSC_FALSE,exSouth=PETSC_FALSE;

			  // check for partially filled/empty grid cell neighbors (mask not updated yet, but vH is)
	          if (vH(i+1,j)==0.0 && vbed(i+1,j)<currentSeaLevel) {countEmptyNeighbors+=1; exEast=PETSC_TRUE;} 
	          if (vH(i-1,j)==0.0 && vbed(i-1,j)<currentSeaLevel) {countEmptyNeighbors+=1; exWest=PETSC_TRUE;}
	          if (vH(i,j+1)==0.0 && vbed(i,j+1)<currentSeaLevel) {countEmptyNeighbors+=1; exNorth=PETSC_TRUE;}
	          if (vH(i,j-1)==0.0 && vbed(i,j-1)<currentSeaLevel) {countEmptyNeighbors+=1; exSouth=PETSC_TRUE;}
	
			  if (countEmptyNeighbors>0 && vH(i,j)>minHRedist)  {
			    //remainder ice mass will be redistributed equally to all adjacent imfrac boxes (is there a more physical way?)
			    if (exEast) vHref(i+1,j)+=vHresidual(i,j)/countEmptyNeighbors;
			    if (exWest) vHref(i-1,j)+=vHresidual(i,j)/countEmptyNeighbors;
			    if (exNorth) vHref(i,j+1)+=vHresidual(i,j)/countEmptyNeighbors;
			    if (exSouth) vHref(i,j-1)+=vHresidual(i,j)/countEmptyNeighbors;
			
				//ierr = verbPrintf(3,grid.com,"!!! Hresidual has been redistributed to %d neighbors around %d,%d \n",countEmptyNeighbors,i,j ); CHKERRQ(ierr);
			    vHresidualnew(i,j)=0.0;
			  } else {
				vHnew(i,j)+=vHresidual(i,j); // mass conservation, but thick ice at one grid cell possible
				vHresidualnew(i,j)=0.0;
				ierr = verbPrintf(4,grid.com,"!!! PISM WARNING: Hresidual has %d partially filled neighbors, set ice thickness to vHnew=%.2e at %d,%d \n",countEmptyNeighbors,vHnew(i,j),i,j ); CHKERRQ(ierr);
			  }
		    } 
		}
	}
	ierr = vHnew.beginGhostComm(vH); CHKERRQ(ierr);
  	ierr = vHnew.endGhostComm(vH); CHKERRQ(ierr);
	
	
	double 	ocean_rho = config.get("sea_water_density");
	double	ice_rho = config.get("ice_density");
	PetscScalar	Hav;
	PetscScalar	Hcut=0.0;
	for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    	for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
	
			// second step: if neighbors, which gained redistributed ice, also become full, this needs to be redistributed in a repeated loop
			if (vHref(i,j)>0.0) {
				Hav=0.0;
	          	PetscInt countIceNeighbors=0; // counting current full floating ice neighbors (mask not yet updated), and calculate new average Hav

          	  	if (vH(i+1,j)>0.0 && (vbed(i+1,j) < (currentSeaLevel - ice_rho/ocean_rho*vH(i+1,j)))) { Hav+=vH(i+1,j); countIceNeighbors+=1;} 
			  	if (vH(i-1,j)>0.0 && (vbed(i-1,j) < (currentSeaLevel - ice_rho/ocean_rho*vH(i-1,j)))) { Hav+=vH(i-1,j); countIceNeighbors+=1;} 
			  	if (vH(i,j+1)>0.0 && (vbed(i,j+1) < (currentSeaLevel - ice_rho/ocean_rho*vH(i,j+1)))) { Hav+=vH(i,j+1); countIceNeighbors+=1;} 
			  	if (vH(i,j-1)>0.0 && (vbed(i,j-1) < (currentSeaLevel - ice_rho/ocean_rho*vH(i,j-1)))) { Hav+=vH(i,j-1); countIceNeighbors+=1;} 

				if (countIceNeighbors>0){
	            	Hav=Hav/countIceNeighbors;
	
	          		PetscScalar coverageRatio = vHref(i,j)/Hav;
	          		if (coverageRatio>1.0) { // partially filled grid cell is considered to be full
						vHresidualnew(i,j)=vHref(i,j)-Hav;
						Hcut+=vHresidualnew(i,j); // summed up to decide, if methods needs to be run once more 
						vHnew(i,j) += Hav;
			    		vHref(i,j) = 0.0;
						//ierr = verbPrintf(3,grid.com,"!!! Partially filled grid cell became full after redistribution at %d,%d \n",i,j ); CHKERRQ(ierr);
					} 
			  	} else { // no full floating ice neighbor 
					vHnew(i,j) += vHref(i,j); // mass conservation, but thick ice at one grid cell possible
					vHref(i,j) = 0.0;
					vHresidualnew(i,j)=0.0;
					ierr = verbPrintf(4, grid.com,"!!! PISM_WARNING: No floating ice neighbors to calculate Hav, set ice thickness to vHnew=%.2e at %d,%d \n",vHnew(i,j),i,j); CHKERRQ(ierr);
			  	} 
		   	}
		}
	}

    PetscScalar gHcut; //check, if redistribution should be run once more
    ierr = PetscGlobalSum(&Hcut, &gHcut, grid.com); CHKERRQ(ierr);
	if (gHcut>0.0) { repeatRedist=PETSC_TRUE;}
    else {			repeatRedist=PETSC_FALSE;}
	//ierr = verbPrintf(3, grid.com,"!!! Hcut=%f \n",gHcut); CHKERRQ(ierr);

	ierr = vH.end_access(); CHKERRQ(ierr);
  	ierr = vHnew.end_access(); CHKERRQ(ierr);
  	// finally copy vHnew into vH and communicate ghosted values
  	ierr = vHnew.beginGhostComm(vH); CHKERRQ(ierr);
  	ierr = vHnew.endGhostComm(vH); CHKERRQ(ierr);

  	ierr = vHref.end_access(); CHKERRQ(ierr); 
	ierr = vbed.end_access(); CHKERRQ(ierr);
	
	ierr = vHresidual.end_access(); CHKERRQ(ierr);
	ierr = vHresidualnew.end_access(); CHKERRQ(ierr);
  	ierr = vHresidualnew.beginGhostComm(vHresidual); CHKERRQ(ierr);
  	ierr = vHresidualnew.endGhostComm(vHresidual); CHKERRQ(ierr);


  return 0;
}
