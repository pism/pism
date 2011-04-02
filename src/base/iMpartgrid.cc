// Copyright (C) 2011 Torsten Albrecht and Ed Bueler
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


// methods implementing PIK logic for -part_grid; see Albrecht et al 2011


//! Compute staggered grid velocities according to mask and regular grid velocities.
/*!
In the finite volume interpretation, these are normal velocities at the faces
of the cell.  The method avoids differencing velocities from ice free ocean locations.

The outputs are velE, velW, velN, velS.
 */
PetscErrorCode IceModel::velsPartGrid(PetscReal Mo, 
                                      PetscReal Me, PetscReal Mw,
                                      PetscReal Mn, PetscReal Ms,
                                      PISMVector2 vrego,
                                      PISMVector2 vrege, PISMVector2 vregw,
                                      PISMVector2 vregn, PISMVector2 vregs,
                                      PetscReal &velE, PetscReal &velW,
                                      PetscReal &velN, PetscReal &velS) {
  const bool oneneighboricefree = (Me > MASK_FLOATING ||
                                   Mw > MASK_FLOATING ||
                                   Mn > MASK_FLOATING ||
                                   Ms > MASK_FLOATING),
             oneneighboricefilled = (Me <= MASK_FLOATING ||
                                     Mw <= MASK_FLOATING ||
                                     Mn <= MASK_FLOATING ||
                                     Ms <= MASK_FLOATING);

  //case1: [i][j] in the middle of ice or bedrock: default scheme
  if (Mo <= MASK_FLOATING && (!oneneighboricefree)) {
    // compute (i,j)-centered "face" velocity components by average
    velE = 0.5 * (vrego.u + vrege.u);
    velW = 0.5 * (vregw.u + vrego.u);
    velN = 0.5 * (vrego.v + vregn.v);
    velS = 0.5 * (vregs.v + vrego.v);
  //case2: [i][j] on floating or grounded ice, but next to a ice-free ocean grid cell
  } else if (Mo <= MASK_FLOATING && (oneneighboricefree)) {
    velE = (Me > MASK_FLOATING ? vrego.u : 0.5 * (vrego.u + vrege.u));
    velW = (Mw > MASK_FLOATING ? vrego.u : 0.5 * (vregw.u + vrego.u));
    velN = (Mn > MASK_FLOATING ? vrego.v : 0.5 * (vrego.v + vregn.v));
    velS = (Ms > MASK_FLOATING ? vrego.v : 0.5 * (vregs.v + vrego.v));
  //case3: [i][j] on ice-free ocean (or partially filled), but next to ice grid cell
  } else if (Mo > MASK_FLOATING && (oneneighboricefilled)) {
    velE = (Me <= MASK_FLOATING ? vrege.u : 0.0);
    velW = (Mw <= MASK_FLOATING ? vregw.u : 0.0);
    velN = (Mn <= MASK_FLOATING ? vregn.v : 0.0);
    velS = (Ms <= MASK_FLOATING ? vregs.v : 0.0);
  //case4: [i][j] on ice-free ocean, and no ice neighbors, and else
  } else {		
    velE = 0.0;
    velW = 0.0;
    velN = 0.0;
    velS = 0.0;
  }

  return 0;
}


//! For ice-free (or partially-filled) cells adjacent to "full" floating cells, update Href.
/*!
Should only be called if one of the neighbors is floating, i.e. only if at
least one of Me,Mw,Mn,Ms is MASK_FLOATING.

Variable Href is modified by this procedure; its input value matters.
Variable Hav is an output.
 */
PetscReal IceModel::getHav(bool do_redist,
                           PetscReal Me, PetscReal Mw, PetscReal Mn, PetscReal Ms,
                           PetscReal He, PetscReal Hw, PetscReal Hn, PetscReal Hs) {
  // FIXME:  in this form does not account for grounded tributaries: thin
  //         ice shelves may evolve from grounded tongue
  // get mean ice thickness over adjacent floating ice shelf neighbors
  PetscReal Hav = 0.0;
  PetscInt countIceNeighbors=0;
  if (Me == MASK_FLOATING) { Hav += He; countIceNeighbors++; } 
  if (Mw == MASK_FLOATING) { Hav += Hw; countIceNeighbors++; }
  if (Mn == MASK_FLOATING) { Hav += Hn; countIceNeighbors++; }
  if (Ms == MASK_FLOATING) { Hav += Hs; countIceNeighbors++; }
  if (countIceNeighbors == 0) {
    SETERRQ(1,"countIceNeighbors == 0;  call me only if a neighbor is floating!\n");
  }
  Hav = Hav / countIceNeighbors;
  // reduces the guess at the front
  if (do_redist) {	
    const PetscReal  mslope = 2.4511e-18*grid.dx/(300*600/secpera);
    // for declining front C/Q0 according to analytical flowline profile in
    //   vandeveen with v0=300m/yr and H0=600m	    
    Hav -= 0.8*mslope*pow(Hav,5);
  }
  return Hav;
}


//! Redistribute residual ice mass from subgrid-scale parameterization, when using -part_redist option.
/*!
See [\ref Albrechtetal2011].  Manages the loop.

FIXME: Reporting!

FIXME: repeatRedist should be config flag?

FIXME: resolve fixed number (=3) of loops issue
 */
PetscErrorCode IceModel::redistResiduals() {
  PetscErrorCode ierr;
  ierr = calculateRedistResiduals(); CHKERRQ(ierr); //while loop?
  PetscInt loopcount=0;
  const PetscInt max_loopcount = 3;
  while ((repeatRedist==PETSC_TRUE) && (loopcount < max_loopcount)) {
    ierr = calculateRedistResiduals(); CHKERRQ(ierr);
    loopcount+=1;
    ierr = verbPrintf(4,grid.com, "redistribution loopcount = %d\n",loopcount); CHKERRQ(ierr);
  }
  return 0;
}


// This routine carries-over the ice mass when using -part_redist option, one step in the loop.
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


// FIXME: following is deprecated and should be removed once its pieces are no longer needed
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
   	   	  if (config.get_flag("part_redist")) {	
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

