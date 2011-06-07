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
#include "Mask.hh"

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

  if (ocean == PETSC_NULL) {  SETERRQ(1, "PISM ERROR: ocean == PETSC_NULL");  }
  PetscReal sea_level;
  ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);

  GeometryCalculator gc(sea_level, *ice, config);
  MaskQuery mask(vMask);

  ierr =    vH.begin_access();    CHKERRQ(ierr);
  ierr =  vbed.begin_access();  CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);

  PetscInt GHOSTS = 2;
  for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
      vMask(i, j) = gc.mask(vbed(i, j), vH(i,j));
    } // inner for loop (j)
  } // outer for loop (i)

  ierr =         vH.end_access(); CHKERRQ(ierr);
  ierr =       vbed.end_access(); CHKERRQ(ierr);
  ierr =      vMask.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::update_surface_elevation() {
  PetscErrorCode ierr;

  MaskQuery mask(vMask);

  if (ocean == PETSC_NULL) {  SETERRQ(1, "PISM ERROR: ocean == PETSC_NULL");  }
  PetscReal sea_level;
  ierr = ocean->sea_level_elevation(sea_level); CHKERRQ(ierr);

  GeometryCalculator gc(sea_level, *ice, config);

  ierr =    vh.begin_access();    CHKERRQ(ierr);
  ierr =    vH.begin_access();    CHKERRQ(ierr);
  ierr =  vbed.begin_access();  CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);

  PetscInt GHOSTS = 2;
  for (PetscInt   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (PetscInt j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
      // take this opportunity to check that vH(i, j) >= 0
      if (vH(i, j) < 0) {
        SETERRQ2(1, "Thickness negative at point i=%d, j=%d", i, j);
      }
      vh(i, j) = gc.surface(vbed(i, j), vH(i, j));
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
  } else { SETERRQ(1, "PISM ERROR: surface == NULL"); }

  if (ocean != NULL) {
    ierr = ocean->shelf_base_mass_flux(shelfbmassflux); CHKERRQ(ierr);
  } else { SETERRQ(2, "PISM ERROR: ocean == NULL"); }

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
  ierr = vel_advective->begin_access(); CHKERRQ(ierr);
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
      //        was not sufficient to zero - out vHresidual already
      ierr = vHresidual.set(0.0); CHKERRQ(ierr);
    }
  }
  const bool dirichlet_bc = config.get_flag("dirichlet_bc");
  if (dirichlet_bc) {
    ierr = vBCMask.begin_access();  CHKERRQ(ierr);
    ierr = vBCvel.begin_access();  CHKERRQ(ierr);
    ierr = vbed.begin_access();  CHKERRQ(ierr);
  }

  if (do_ocean_kill) {
    ierr = ocean_kill_mask.begin_access(); CHKERRQ(ierr);
  }

  MaskQuery mask(vMask);

  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {

      PetscScalar divQ = 0.0;

      if (mask.grounded(i, j)) {
        // staggered grid Div(Q) for diffusive non-sliding SIA deformation part:
        //    Qdiff = - D grad h
        divQ = ((*Qdiff)(i, j, 0) - (*Qdiff)(i - 1, j, 0)) / dx
          + ((*Qdiff)(i, j, 1) - (*Qdiff)(i, j - 1, 1)) / dy;
      }

      // get non-diffusive velocities according to old or -part_grid scheme
      planeStar<int> M = vMask.int_star(i, j);
      planeStar<PetscScalar> v;

      if (dirichlet_bc) {
        //the staggered velocities have to be adjusted to Dirichlet boundary conditions
        if (vBCMask.as_int(i,j) == 0) {
          if (vBCMask.as_int(i+1,j) == 1) v.e = vBCvel(i + 1, j).u;
          if (vBCMask.as_int(i-1,j) == 1) v.w = vBCvel(i - 1, j).u;
          if (vBCMask.as_int(i,j+1) == 1) v.n = vBCvel(i, j + 1).v;
          if (vBCMask.as_int(i,j-1) == 1) v.s = vBCvel(i, j - 1).v;
        }
      } else {
        ierr = cell_interface_velocities(do_part_grid, i, j, v); CHKERRQ(ierr);
      }

      // membrane stress (and/or basal sliding) part: upwind by staggered grid
      // PIK method;  this is   \nabla \cdot [(u, v) H]
      divQ += (  v.e * (v.e > 0 ? vH(i, j) : vH(i + 1, j))
                 - v.w * (v.w > 0 ? vH(i - 1, j) : vH(i, j)) ) / dx;
      divQ += (  v.n * (v.n > 0 ? vH(i, j) : vH(i, j + 1))
                 - v.s * (v.s > 0 ? vH(i, j - 1) : vH(i, j)) ) / dy;


      PetscReal S = 0.0;
      if (include_bmr_in_continuity) {
        if (mask.ocean(i, j))
          S = shelfbmassflux(i,j);
        else
          S = bmr_gnded[i][j];
      }

      // decide whether to apply Albrecht et al 2011 subgrid-scale parameterization (-part_grid)

      // case where we apply -part_grid
	  //if (do_part_grid && mask.next_to_floating_ice(i, j) && !mask.next_to_grounded_ice(i, j)) {
      if (do_part_grid && mask.next_to_floating_ice(i, j)) {
        vHref(i, j) -= divQ * dt;
		if (vHref(i, j) < 0.0) { 
			my_nonneg_rule_flux += ( - vHref(i, j));
			vHref(i, j) = 0.0;
			ierr = verbPrintf(2, grid.com,"!!! PISM_WARNING: vHref is negative at i=%d, j=%d\n",i,j); CHKERRQ(ierr);
		}

        PetscReal H_average = get_average_thickness(do_redist, M, vH.star(i, j));

        // To calculate the surface balance contribution with respect to the
        // coverage ratio, let  X = vHref_new  be the new value of Href.  We assume
        //   X = vHref_old + (M - S) * dt * coverageRatio
        // equivalently
        //   X = vHref_old + (M - S) * dt * X / H_average.
        // where M = acab and S = shelfbaseflux for floating ice.  Solving for X we get
        //   X = vHref_old / (1.0 - (M - S) * dt * H_average))
		/*
        if ((acab(i, j) - S) * dt < H_average) {
			vHref(i, j) = vHref(i, j) / (1.0 - (acab(i, j) - S) * dt / H_average);
		} else {
			ierr = verbPrintf(4, grid.com,"!!! PISM_WARNING: H_average is smaller than surface mass balance at i=%d, j=%d.\n",i,j); CHKERRQ(ierr);
		}
		*/

        const PetscScalar coverageRatio = vHref(i, j) / H_average;

        if (coverageRatio > 1.0) { // partially filled grid cell is considered to be full
          if (do_redist) {  vHresidual(i, j) = vHref(i, j) - H_average;  } //residual ice thickness
          vHnew(i, j) = H_average; // gets a "real" ice thickness
	   	  vHnew(i, j)+= (acab(i, j) - S) * dt; // no implicit SMB in partially filled cells any more
          vHref(i, j) = 0.0;
        } else {
          vHnew(i, j) = 0.0; // no change from vH value, actually
          // vHref(i, j) not changed
        }

      } else if (mask.grounded(i, j) ||
                 mask.floating_ice(i, j) ||
                 mask.next_to_grounded_ice(i, j) ) {
        // grounded/floating default case, and case of ice-free ocean adjacent to grounded
        vHnew(i, j) += (acab(i, j) - S - divQ) * dt;
      } else {
        // last possibility: ice-free, not adjacent to a "full" cell at all
        vHnew(i, j) = 0.0;
      }
      
      if (dirichlet_bc && vBCMask.as_int(i,j) == 1) {
        vHnew(i, j) = vH(i,j);
      }

      // apply free boundary rule: negative thickness becomes zero
      if (vHnew(i, j) < 0) {
        my_nonneg_rule_flux += ( - vHnew(i, j));
        vHnew(i, j) = 0.0;
      }

      // the following conditionals, both -ocean_kill and -float_kill, are also applied in
      //   IceModel::computeMax2DSlidingSpeed() when determining CFL

      // force zero thickness at points which were originally ocean (if "-ocean_kill");
      //   this is calving at original calving front location
      if ( do_ocean_kill && ocean_kill_mask.as_int(i, j) == 1) {
        my_ocean_kill_flux -= vHnew(i, j);
        vHnew(i, j) = 0.0;
      }

      // force zero thickness at points which are floating (if "-float_kill");
      //   this is calving at grounding line
      if ( floating_ice_killed && mask.ocean(i, j) ) {
        my_float_kill_flux -= vHnew(i, j);
        vHnew(i, j) = 0.0;
      }

    } // end of the inner for loop
  } // end of the outer for loop

  ierr = vbmr.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = Qdiff->end_access(); CHKERRQ(ierr);
  ierr = vel_advective->end_access(); CHKERRQ(ierr);
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

  if (dirichlet_bc) {
    ierr = vBCMask.end_access();  CHKERRQ(ierr);
    ierr = vBCvel.end_access();  CHKERRQ(ierr);
    ierr = vbed.end_access();  CHKERRQ(ierr);
  }

  if (do_ocean_kill) {
    ierr = ocean_kill_mask.end_access(); CHKERRQ(ierr);
  }

  // flux accounting
  {
    ierr = PetscGlobalSum(&my_nonneg_rule_flux, &nonneg_rule_flux, grid.com); CHKERRQ(ierr);
    ierr = PetscGlobalSum(&my_ocean_kill_flux, &ocean_kill_flux, grid.com); CHKERRQ(ierr);
    ierr = PetscGlobalSum(&my_float_kill_flux, &float_kill_flux, grid.com); CHKERRQ(ierr);

    // FIXME: use corrected cell areas (when available)
    PetscScalar ice_density = config.get("ice_density"),
      factor = ice_density * (dx * dy) / dt;
    nonneg_rule_flux *= factor;
    ocean_kill_flux *= factor;
    float_kill_flux *= factor;
  } //FIXME: flux reporting not yet adjusted to part_grid scheme

  // compute dH/dt (thickening rate) for viewing and for saving at end; only diagnostic
  ierr = vHnew.add(-1.0, vH, vdHdt); CHKERRQ(ierr); // vdHdt = vHnew - vH
  ierr = vdHdt.scale(1.0 / dt); CHKERRQ(ierr);        // vdHdt = vdHdt / dt

  // d(volume) / dt
  {
    PetscScalar dvol = 0.0;

    ierr = vdHdt.begin_access(); CHKERRQ(ierr);
    ierr = cell_area.begin_access(); CHKERRQ(ierr);
    for (PetscInt i = grid.xs; i < grid.xs + grid.xm; ++i) {
      for (PetscInt j = grid.ys; j < grid.ys + grid.ym; ++j) {
        dvol += vdHdt(i, j) * cell_area(i, j);
      }
    }
    ierr = cell_area.end_access(); CHKERRQ(ierr);
    ierr = vdHdt.end_access(); CHKERRQ(ierr);

    ierr = PetscGlobalSum(&dvol, &dvoldt, grid.com); CHKERRQ(ierr);
  }

  // average value of dH / dt;
  PetscScalar ice_area;
  ierr = compute_ice_area(ice_area); CHKERRQ(ierr);

  ierr = vdHdt.sum(gdHdtav); CHKERRQ(ierr);
  gdHdtav = gdHdtav / ice_area; // m/s

  // finally copy vHnew into vH and communicate ghosted values
  ierr = vHnew.beginGhostComm(vH); CHKERRQ(ierr);
  ierr = vHnew.endGhostComm(vH); CHKERRQ(ierr);

  // the following calls are new routines adopted from PISM-PIK. The place and
  // order is not clear yet!

  // There is no reporting of single ice fluxes yet in comparison to total ice
  // thickness change.

  // distribute residual ice mass if desired
  if (do_redist) {
    ierr = redistResiduals(); CHKERRQ(ierr);
  }

  // FIXME: maybe calving should be applied *before* the redistribution part?
  if (config.get_flag("do_eigen_calving") && config.get_flag("use_ssa_velocity")) {
    ierr = stress_balance->get_principal_strain_rates(vPrinStrain1, vPrinStrain2); CHKERRQ(ierr);
    ierr = eigenCalving(); CHKERRQ(ierr);
  }

  if (config.get_flag("do_thickness_calving")) {
    if (config.get_flag("part_grid") == true) {
      ierr = calvingAtThickness(); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(4, grid.com,
                        "PISM - WARNING: calving at certain terminal ice thickness without application\n"
                        "              of partially filled grid cell scheme may lead to non-moving\n"
                        "              ice shelf front!\n"); CHKERRQ(ierr);
    }
  }

  // Check if the ice thickness exceeded the height of the computational box
  // and extend the grid if necessary:
  ierr = check_maximum_thickness(); CHKERRQ(ierr);

  return 0;
}

