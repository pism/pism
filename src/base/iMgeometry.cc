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

  if (config.get_flag("kill_icebergs")) {
    ierr = killIceBergs(); CHKERRQ(ierr);	
  }

  return 0;
}


PetscErrorCode IceModel::update_mask() {
  PetscErrorCode ierr;

  if (ocean == PETSC_NULL) {  SETERRQ(1,"PISM ERROR: ocean == PETSC_NULL");  }
  PetscReal currentSeaLevel;
  ierr = ocean->sea_level_elevation(currentSeaLevel); CHKERRQ(ierr);

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
  ierr = ocean->sea_level_elevation(currentSeaLevel); CHKERRQ(ierr);

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

The PISMSurfaceModel IceModel::surface provides \f$M\f$.  The
PISMOceanModel IceModel::ocean provides \f$S\f$ at locations below
floating ice (ice shelves).

Even if we regard the map-plane flux as defined by the formula
\f$\mathbf{q} = \bar{\mathbf{U}} H\f$, the flow of ice is at least somewhat
diffusive in almost all cases, and in simplified models (the SIA model) it
is exactly true that \f$\mathbf{q} = - D \nabla h\f$.  In the current method the
flux is split into the part from the diffusive non-sliding SIA model
and a part which is a less-diffusive, presumably membrane-stress-dominated
2D advective velocity, which generally describes sliding:
  \f[ \mathbf{q} = - D \nabla h + \mathbf{U}_b H.\f]
Here \f$D\f$ is the (positive, scalar) effective diffusivity of the non-sliding
SIA and \f$\mathbf{U}_b\f$ is the less-diffusive sliding velocity.
We interpret \f$\mathbf{U}_b\f$ as a basal sliding velocity in the hybrid.

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
the PIK upwinding technique [equation (25) in \ref Winkelmannetal2010TCD].
The CFL condition for this advection scheme is checked; see 
computeMax2DSlidingSpeed() and determineTimeStep().  This method implements the
direct-superposition (PIK) hybrid which adds the SSA velocity to the SIA velocity
[equation (15) in \ref Winkelmannetal2010TCD].  The hybrid described by equations
(21) and (22) in \ref BBL is no longer used.

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
    ierr = surface->ice_surface_mass_flux(acab); CHKERRQ(ierr);
  } else { SETERRQ(1,"PISM ERROR: surface == NULL"); }

  if (ocean != NULL) {
    ierr = ocean->shelf_base_mass_flux(shelfbmassflux); CHKERRQ(ierr);
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

  // related to PIK part_grid mechanism; see Albrecht et al 2011
  const bool do_part_grid = config.get_flag("part_grid"),
                do_redist = config.get_flag("part_redist");
  if (do_part_grid) {
    ierr = vHref.begin_access(); CHKERRQ(ierr);
    if (do_redist) {
      ierr = vHresidual.begin_access(); CHKERRQ(ierr);
      // FIXME: next line causes mass loss if max_loopcount in redistResiduals()
      //        was not sufficient to zero-out vHresidual already
      ierr = vHresidual.set(0.0); CHKERRQ(ierr);
    }
  }

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      PetscScalar divQ = 0.0;

      if (vMask.is_grounded(i,j)) {
        // staggered grid Div(Q) for diffusive non-sliding SIA deformation part:
        //    Qdiff = - D grad h
        divQ =   ((*Qdiff)(i,j,0) - (*Qdiff)(i-1,j,0)) / dx
               + ((*Qdiff)(i,j,1) - (*Qdiff)(i,j-1,1)) / dy;
      }

      // get non-diffusive velocities according to old or -part_grid scheme
      const int Mo = vMask.value(i,j),
                Me = vMask.value(i+1,j), Mw = vMask.value(i-1,j),
                Mn = vMask.value(i,j+1), Ms = vMask.value(i,j-1);
      PetscScalar velE, velW, velN, velS;
      if (do_part_grid) {
          ierr = velsPartGrid(Mo, Me, Mw, Mn, Ms,
                              vel(i,j), vel(i+1,j), vel(i-1,j), vel(i,j+1), vel(i,j-1),
                              velE, velW, velN, velS); CHKERRQ(ierr);
      } else {
          // just compute (i,j)-centered "face" velocity components by average
          velE = 0.5 * (vel(i,j).u + vel(i+1,j).u),
          velW = 0.5 * (vel(i-1,j).u + vel(i,j).u),
          velN = 0.5 * (vel(i,j).v + vel(i,j+1).v),
          velS = 0.5 * (vel(i,j-1).v + vel(i,j).v);
      }

      // membrane stress (and/or basal sliding) part: upwind by staggered grid
      // PIK method;  this is   \nabla \cdot [(u,v) H]
      divQ += (  velE * (velE > 0 ? vH(i,j) : vH(i+1,j))
               - velW * (velW > 0 ? vH(i-1,j) : vH(i,j)) ) / dx;
      divQ += (  velN * (velN > 0 ? vH(i,j) : vH(i,j+1))
               - velS * (velS > 0 ? vH(i,j-1) : vH(i,j)) ) / dy;

      // decide whether to apply Albrecht et al 2011 subgrid-scale parameterization (-part_grid)
      const bool adjacenttofloating = (Me == MASK_FLOATING ||
                                       Mw == MASK_FLOATING ||
                                       Mn == MASK_FLOATING ||
                                       Ms == MASK_FLOATING),
                 adjacenttogrounded = (Me < MASK_FLOATING ||
                                       Mw < MASK_FLOATING ||
                                       Mn < MASK_FLOATING ||
                                       Ms < MASK_FLOATING);
      // case where we apply -part_grid
      if ((do_part_grid) && (Mo > MASK_FLOATING) && (adjacenttofloating)) {
        vHref(i,j) -= divQ * dt;
        PetscReal Hav = getHav(do_redist, Me, Mw, Mn, Ms,
                               vH(i+1,j), vH(i-1,j), vH(i,j+1), vH(i,j-1));
        // To calculate the surface balance contribution with respect to the 
        // coverage ratio, let  X = vHref_new  be the new value of Href.  We assume
        //   X = vHref_old + (M - S) * dt * coverageRatio
        // equivalently
        //   X = vHref_old + (M - S) * dt * X / Hav.
        // where M = acab and S = shelfbaseflux for floating ice.  Solving for X we get
        //   X = vHref_old / (1.0 - (M - S) * dt * Hav))
        const PetscReal  MminusS = acab(i,j) - (include_bmr_in_continuity ? shelfbmassflux(i,j) : 0.0);
        vHref(i,j) = vHref(i,j) / (1.0 - MminusS * dt / Hav);
        const PetscScalar coverageRatio = vHref(i,j) / Hav;
        if (coverageRatio > 1.0) { // partially filled grid cell is considered to be full
          if (do_redist) {  vHresidual(i,j) = vHref(i,j) - Hav;  } //residual ice thickness
          vHnew(i,j) = Hav; // gets a "real" ice thickness
          vHref(i,j) = 0.0;
        } else {
          vHnew(i,j) = 0.0; // no change from vH value, actually
          // vHref(i,j) not changed
        }

      // grounded/floating default case, and case of ice-free ocean adjacent to grounded
      } else if ( (Mo <= MASK_FLOATING) ||
                  ((Mo > MASK_FLOATING) && (adjacenttogrounded)) ) {
        vHnew(i,j) += (acab(i,j) - divQ) * dt; // include M
        if (include_bmr_in_continuity) { // include S
          if (vMask.is_floating(i,j)) {
            vHnew(i,j) -= shelfbmassflux(i,j) * dt;
          } else {
            vHnew(i,j) -= bmr_gnded[i][j] * dt;
          }
        }

      // last possibility: ice-free, not adjacent to a "full" cell at all
      } else {  vHnew(i,j)=0.0;  }

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

  if (do_part_grid) {
    ierr = vHref.end_access(); CHKERRQ(ierr);
    if (do_redist) {
      ierr = vHresidual.end_access(); CHKERRQ(ierr);
    }
  }
  
  // flux accounting
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
  } //FIXME: flux reporting not yet adjusted to part_grid scheme

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

  // the remaining calls are new routines adopted from PISM-PIK. The place and order is not clear yet!

  // There is no reporting of single ice fluxes yet in comparison to total ice thickness change.

  // distribute residual ice mass if desired
  if (do_redist) {
    ierr = redistResiduals(); CHKERRQ(ierr);
  }

  // FIXME: maybe calving should be applied *before* the redistribution part?
  if (config.get_flag("do_eigen_calving")) {
    ierr = stress_balance->get_principle_strain_rates( vPrinStrain1,vPrinStrain2); CHKERRQ(ierr);
    ierr = eigenCalving(); CHKERRQ(ierr);
  }
  if (config.get_flag("do_thickness_calving")) { 
    if (config.get_flag("part_grid")==true) { 
      ierr = calvingAtThickness(); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(4,grid.com,
        "PISM-WARNING: calving at certain terminal ice thickness without application\n"
        "              of partially filled grid cell scheme may lead to non-moving\n"
        "              ice shelf front!\n"); CHKERRQ(ierr);
    }
  }
	
  return 0;
}

