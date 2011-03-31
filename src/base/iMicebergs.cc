// Copyright (C) 2004--2011 Torsten Albrecht
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
#include <petscda.h>
#include "iceModel.hh"


//! Identify and eliminate free-floating icebergs, which cause well-posedness problems for stress solvers.
/*!
Icebergs are, in this context, floating regions that are \e not attached, through a chain of positive thickness ice-filled cells, to at least one grounded cell.  They are observed to cause unrealistically large velocities that may (numerically) affect the ice velocities everywhere.  They cause the SSA operator to have a nontrivial null space, or, under approximation errors, they lead to extremely small time steps and can eventually cause a KSP-ERROR.

This method calls the routines is which first identify and then eliminate these icebergs.

FIXME:  a fundamental aspect of the semantics here is not clear to me (bueler), namely how many times the iceberg-eliminate "sweep" might occur, and what properties control that?  for now, should there be some (low-verbosity) indication that it is occurring, such as when more than one sweep happened?

FIXME:  this package of methods *might* appropriately be a class

FIXME:  this package of routines *should* have a regression
*/
PetscErrorCode IceModel::killIceBergs() {
  PetscErrorCode ierr;

  ierr = findIceBergCandidates(); CHKERRQ(ierr);
  ierr = identifyNotAnIceBerg(); CHKERRQ(ierr);
  ierr = killIdentifiedIceBergs(); CHKERRQ(ierr);
  if (config.get_flag("do_eigen_calving") || config.get_flag("do_thickness_calving")) {
    ierr = killEasyIceBergs(); CHKERRQ(ierr);
  }
  return 0;
}


//! This routine comes from PISM-PIK.
/*!
The aim of this routine is to find floating regions that *might* be icebergs. If these regions actually *are* icebergs is checked in identifyNotAnIceBerg().
*/
PetscErrorCode IceModel::findIceBergCandidates() {
  PetscErrorCode ierr;

  ierr = verbPrintf(4,grid.com,"######### findIceBergCandidates is called \n");    CHKERRQ(ierr);

  const PetscInt Mx = grid.Mx, My = grid.My;
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vIcebergMask.begin_access(); CHKERRQ(ierr);
  ierr = vbed.begin_access(); CHKERRQ(ierr);

  //const PetscScalar   a = grid.dx * grid.dy * 1e-3 * 1e-3; // area unit (km^2)
  PetscReal currentSeaLevel;
  if (ocean != NULL) {
    ierr = ocean->sea_level_elevation(grid.year, dt / secpera, currentSeaLevel); CHKERRQ(ierr);
  } else { SETERRQ(2,"PISM ERROR: ocean == NULL"); }
  double 	ocean_rho = config.get("sea_water_density"),
			ice_rho = config.get("ice_density");

	
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

	    const PetscScalar hgrounded = vbed(i,j) + vH(i,j),
		hfloating = currentSeaLevel + (1.0 - ice_rho/ocean_rho) * vH(i,j);
	
      //cut of border of computational domain	
      if (hgrounded<hfloating &&  (i <= 0 || i >= Mx-1 || j <= 0 || j >= My-1)) {
    	//if ((i <= 0 || i >= Mx-1 || j <= 0 || j >= My-1)) {
          vH(i,j) = 0.0;
          vIcebergMask(i,j) = ICEBERGMASK_STOP_OCEAN;
          vMask(i,j) = MASK_ICE_FREE_OCEAN;	
      }else{
        vIcebergMask(i,j) = ICEBERGMASK_NOT_SET;
      }
    }
  }

  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr); 
  ierr = vIcebergMask.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);

  ierr = vMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vMask.endGhostComm(); CHKERRQ(ierr);
  ierr = vIcebergMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vIcebergMask.endGhostComm(); CHKERRQ(ierr);

  // set all floating points to ICEBERGMASK_ICEBERG_CAND
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vIcebergMask.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      if (vIcebergMask(i,j) == ICEBERGMASK_NOT_SET) {
        if (vMask(i,j) == MASK_FLOATING){
          vIcebergMask(i,j) = ICEBERGMASK_ICEBERG_CAND;
        }

      }
    }
  }
  ierr = vMask.end_access(); CHKERRQ(ierr); 
  ierr = vIcebergMask.end_access(); CHKERRQ(ierr);

  // set borders of shelves/icebergs to ICEBERGMASK_STOP_ATTACHED or ICEBERGMASK_STOP_OCEAN respectively.
  ierr = vIcebergMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vIcebergMask.endGhostComm(); CHKERRQ(ierr);

  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vIcebergMask.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      if (vIcebergMask(i,j) == ICEBERGMASK_NOT_SET) {
		
        bool neighbor_is_candidate = (	vIcebergMask(i+1,j)== ICEBERGMASK_ICEBERG_CAND ||
                                        vIcebergMask(i+1,j+1)== ICEBERGMASK_ICEBERG_CAND ||
                                        vIcebergMask(i+1,j-1)== ICEBERGMASK_ICEBERG_CAND ||
                                        vIcebergMask(i,j+1)== ICEBERGMASK_ICEBERG_CAND ||
                                        vIcebergMask(i,j-1)== ICEBERGMASK_ICEBERG_CAND ||
                                        vIcebergMask(i-1,j+1)== ICEBERGMASK_ICEBERG_CAND ||
                                        vIcebergMask(i-1,j)== ICEBERGMASK_ICEBERG_CAND ||
                                        vIcebergMask(i-1,j-1)== ICEBERGMASK_ICEBERG_CAND);
											
        if (vMask(i,j)<MASK_FLOATING && neighbor_is_candidate) vIcebergMask(i,j) = ICEBERGMASK_STOP_ATTACHED;
        else if (vMask(i,j)>MASK_FLOATING && neighbor_is_candidate) vIcebergMask(i,j) = ICEBERGMASK_STOP_OCEAN;		

      }
    }
  }

  ierr = vMask.end_access(); CHKERRQ(ierr); 
  ierr = vIcebergMask.end_access(); CHKERRQ(ierr);

  ierr = vIcebergMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vIcebergMask.endGhostComm(); CHKERRQ(ierr);

  return 0;
}




PetscErrorCode IceModel::identifyNotAnIceBerg() {
  PetscErrorCode ierr;
  
  ierr = verbPrintf(4,grid.com,"######### identifyNotAnIceBerg is called \n");    CHKERRQ(ierr);
  
  // this communication of ghostvalues is done here to make sure that asking about neighbouring values in this routine 
  // doesn't lead to inconsistencies in parallel computation, if the neighbour belongs to another processor domain.
  ierr = vIcebergMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vIcebergMask.endGhostComm(); CHKERRQ(ierr);
  
  PetscTruth checkingNoIceBergs=PETSC_TRUE;
  PetscInt loopcount = 0;
  while(checkingNoIceBergs==PETSC_TRUE){
    ierr = vIcebergMask.begin_access(); CHKERRQ(ierr);
    
    checkingNoIceBergs = PETSC_FALSE;
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
	
        bool attached_to_grounded = (vIcebergMask(i,j+1) == ICEBERGMASK_STOP_ATTACHED ||
                                     vIcebergMask(i,j-1) == ICEBERGMASK_STOP_ATTACHED ||
                                     vIcebergMask(i+1,j) == ICEBERGMASK_STOP_ATTACHED ||
                                     vIcebergMask(i-1,j) == ICEBERGMASK_STOP_ATTACHED),
		
          attached_to_no_iceberg = (vIcebergMask(i,j+1) == ICEBERGMASK_NO_ICEBERG ||
                                    vIcebergMask(i,j-1) == ICEBERGMASK_NO_ICEBERG ||
                                    vIcebergMask(i+1,j) == ICEBERGMASK_NO_ICEBERG ||
                                    vIcebergMask(i-1,j) == ICEBERGMASK_NO_ICEBERG);
			
        if (vIcebergMask(i,j) == ICEBERGMASK_ICEBERG_CAND && (attached_to_grounded || attached_to_no_iceberg)) {
	
          vIcebergMask(i,j) = ICEBERGMASK_NO_ICEBERG;
          checkingNoIceBergs = PETSC_TRUE;	  
        }
	
      }
    }
    ierr = vIcebergMask.end_access(); CHKERRQ(ierr);
    /* xxxGhostComm() are collective operations. They must be invoked if anything changed on any processor.
     * Thus, collect here the checkingNoIceBergs flags from all processors and compute the global logical OR
     */
    MPI_Allreduce(MPI_IN_PLACE,&checkingNoIceBergs,1,MPI_INT,MPI_LOR,grid.com);
    ierr = vIcebergMask.beginGhostComm(); CHKERRQ(ierr);
    ierr = vIcebergMask.endGhostComm(); CHKERRQ(ierr);
    loopcount+=1;
  }
  ierr = verbPrintf(5,grid.com,"!!! %d loop(s) were needed to identify whether there are icebergs \n",loopcount);    CHKERRQ(ierr);

  return 0;
}


/*!
 We have distinguished icebergs from attached floating regions in identifyNotAnIceBerg(). Now we actually eliminate the former (meaning we set the ice thickness to zero and mark the boxes as icefree-ocean) and leave the latter as they are.
*/
PetscErrorCode IceModel::killIdentifiedIceBergs() {
  PetscErrorCode ierr;
  
  ierr = verbPrintf(4,grid.com,"######### killIceBergs is called \n");    CHKERRQ(ierr);
  
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vIcebergMask.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vh.begin_access(); CHKERRQ(ierr);
  ierr = vbed.begin_access(); CHKERRQ(ierr);

  /* //still don't know why it was there...
  if (ocean == PETSC_NULL) {  SETERRQ(1,"PISM ERROR: ocean == PETSC_NULL");  }
  PetscReal currentSeaLevel=0.0; 
  ierr = ocean->sea_level_elevation(grid.year, dt / secpera, currentSeaLevel); CHKERRQ(ierr);

  bool include_bmr_in_continuity = config.get_flag("include_bmr_in_continuity");

  if (surface != NULL) {
    ierr = surface->ice_surface_mass_flux(grid.year, dt / secpera, acab); CHKERRQ(ierr);
  } else { SETERRQ(1,"PISM ERROR: surface == NULL"); }

  if (ocean != NULL) {
    ierr = ocean->shelf_base_mass_flux(grid.year, dt / secpera, shelfbmassflux); CHKERRQ(ierr);
  } else { SETERRQ(2,"PISM ERROR: ocean == NULL"); }

	double 	ocean_rho = config.get("sea_water_density"),
			ice_rho = config.get("ice_density");
  */

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      
      //const PetscScalar hgrounded = vbed(i,j) + vH(i,j),
	  //	hfloating = currentSeaLevel + (1.0 - ice_rho/ocean_rho) * vH(i,j);
      
      if (vIcebergMask(i,j) == ICEBERGMASK_ICEBERG_CAND){ // actually it's not a candidate any more, it is an iceberg!

		vH(i,j) = 0.0;
		//vh(i,j) = hfloating; //why?
		vh(i,j) = 0.0;
		vMask(i,j) = MASK_ICE_FREE_OCEAN;
		PetscPrintf(PETSC_COMM_SELF,"PISMPIK_INFO: [rank %d] killed iceberg at i=%d,j=%d\n",grid.rank,i,j);
	    // ierr = verbPrintf(5,grid.com,"PISMPIK_INFO killed iceberg at %d,%d\n", i,j);    CHKERRQ(ierr);	
        // killIceBergsWorked=PETSC_TRUE;//not needed any mores since killEasyIceBerg() is called before now
      }
      
    }
  }

  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = vIcebergMask.end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);
  ierr = vh.end_access(); CHKERRQ(ierr);
  ierr = vbed.end_access(); CHKERRQ(ierr);
  
  ierr = vMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vMask.endGhostComm(); CHKERRQ(ierr);
  ierr = vH.beginGhostComm(); CHKERRQ(ierr);
  ierr = vH.endGhostComm(); CHKERRQ(ierr);
  ierr = vh.beginGhostComm(); CHKERRQ(ierr);
  ierr = vh.endGhostComm(); CHKERRQ(ierr);
  
  return 0;
}


/*!
This routine is used when strain-rate based calving is applied. It kills single grid cell icebergs 
and one-grid cell wide ice 'noses', because no proper strain rate eigenvalues can be derived there.
*/
PetscErrorCode IceModel::killEasyIceBergs() {
  PetscErrorCode ierr;
  
  ierr = verbPrintf(4,grid.com,"######### killEasyIceBergs is called \n");    CHKERRQ(ierr);

  //ierr = vH.beginGhostComm(); CHKERRQ(ierr); 
  //ierr = vH.endGhostComm(); CHKERRQ(ierr);

  IceModelVec2S vHnew = vWork2d[0];
  ierr = vH.copy_to(vHnew); CHKERRQ(ierr);
  ierr = vHnew.begin_access(); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);

  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vh.begin_access(); CHKERRQ(ierr);
  ierr = vbed.begin_access(); CHKERRQ(ierr);

  /* only for reporting needed
  bool include_bmr_in_continuity = config.get_flag("include_bmr_in_continuity");
  if (surface != NULL) {
    ierr = surface->ice_surface_mass_flux(grid.year, dt / secpera, acab); CHKERRQ(ierr);
  } else { SETERRQ(1,"PISM ERROR: surface == NULL"); }
  */

  PetscReal currentSeaLevel;
  if (ocean != NULL) {
    //ierr = ocean->shelf_base_mass_flux(grid.year, dt / secpera, shelfbmassflux); CHKERRQ(ierr);
    ierr = ocean->sea_level_elevation(grid.year, dt / secpera, currentSeaLevel); CHKERRQ(ierr);
  } else { SETERRQ(2,"PISM ERROR: ocean == NULL"); }
 
	double 	ocean_rho = config.get("sea_water_density"),
			ice_rho = config.get("ice_density");
  

  /////////////////////////////////////////////////////////////////////////////
  // looking for grid-cell wide floating ice noses that have at least six neighbors 
  // of thickness H=0 like this (o ocean, fl and x floating):
  //   o o o			o o o		  o o o
  //  fl x fl    OR     o x o   OR    o x fl
  //   o o o			o o o		  o o o
  /////////////////////////////////////////////////////////////////////////// 

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      
	  // instead of updating surface elevation, counting here floating or icefree neighbors
      const PetscScalar hgrounded = vbed(i,j) + vH(i,j),
	  hfloating = currentSeaLevel + (1.0 - ice_rho/ocean_rho) * vH(i,j);
      
      if (vH(i,j)>0.0 && hgrounded<hfloating) { //is floating ice shelf
	
		const PetscScalar 
		hgrounded_eb = vbed(i+1,j) + vH(i+1,j),
	  	hfloating_eb = currentSeaLevel + (1.0 - ice_rho/ocean_rho) * vH(i+1,j),
		hgrounded_wb = vbed(i-1,j) + vH(i-1,j),
		hfloating_wb = currentSeaLevel + (1.0 - ice_rho/ocean_rho) * vH(i-1,j),
		hgrounded_nb = vbed(i,j+1) + vH(i,j+1),
	  	hfloating_nb = currentSeaLevel + (1.0 - ice_rho/ocean_rho) * vH(i,j+1),
		hgrounded_sb = vbed(i,j-1) + vH(i,j-1),
	  	hfloating_sb = currentSeaLevel + (1.0 - ice_rho/ocean_rho) * vH(i,j-1);
	
		PetscInt jcount=0, icount=0; // grid-cell wide floating ice nose
		if (vH(i+1,j+1)==0.0) {jcount+=1; icount+=1;}
		if (vH(i+1,j-1)==0.0) {jcount+=1; icount+=1;}
		if (vH(i-1,j+1)==0.0) {jcount+=1; icount+=1;}
		if (vH(i-1,j-1)==0.0) {jcount+=1; icount+=1;}
		if (vH(i+1,j)==0.0) jcount+=1;
		if (vH(i-1,j)==0.0) jcount+=1;
		if (vH(i,j+1)==0.0) icount+=1;
		if (vH(i,j-1)==0.0) icount+=1;
	
		if ((icount == 6 && hgrounded_eb < hfloating_eb &&  hgrounded_wb < hfloating_wb) || 
	        (jcount == 6 && hgrounded_nb < hfloating_nb &&  hgrounded_sb < hfloating_sb)) {
		
	  		vHnew(i,j)=0.0;
	  		//vh(i,j) = hfloating; // why?
			vh(i,j) = 0.0;
	  		vMask(i,j) = MASK_ICE_FREE_OCEAN;
          	PetscPrintf(PETSC_COMM_SELF,"PISMPIK_INFO: [rank %d] cut off nose or one-box-iceberg at i=%d,j=%d\n",grid.rank,i,j);
		}
      }
    }
  }
  
  //ierr =  vbed.end_access(); CHKERRQ(ierr);
  //ierr =  vMask.end_access(); CHKERRQ(ierr);
  //ierr =  vh.end_access(); CHKERRQ(ierr);
  ierr =  vH.end_access(); CHKERRQ(ierr);
  ierr =  vHnew.end_access(); CHKERRQ(ierr);
  // finally copy vHnew into vH and communicate ghosted values
  ierr = vHnew.beginGhostComm(vH); CHKERRQ(ierr);
  ierr = vHnew.endGhostComm(vH); CHKERRQ(ierr);


/////////////////////////////////////////////////////////////////////////////
// looking for one-grid-cell icebergs, that have 4 neighbors of thickness H=0  
/////////////////////////////////////////////////////////////////////////////

  ierr = vH.copy_to(vHnew); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHnew.begin_access(); CHKERRQ(ierr);
  //ierr = vMask.begin_access(); CHKERRQ(ierr);
  //ierr = vh.begin_access(); CHKERRQ(ierr);
  //ierr = vbed.begin_access(); CHKERRQ(ierr);
  
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      
	  // instead of updating surface elevation, counting here floating or icefree neighbors
      const PetscScalar hgrounded = vbed(i,j) + vH(i,j),
	  hfloating = currentSeaLevel + (1.0 - ice_rho/ocean_rho) * vH(i,j);
      
	  bool all_4neighbors_iceless = (vH(i+1,j)==0.0 && vH(i-1,j)==0.0 && vH(i,j+1)==0.0 && vH(i,j-1)==0.0);

      if (vH(i,j)>0.0 && hgrounded < hfloating && all_4neighbors_iceless) { 
	    vHnew(i,j)=0.0;
		//vHnew2(i,j)=0.0;
		//vh(i,j) = hfloating; // why?
		vh(i,j) = 0.0;
	  	vMask(i,j) = MASK_ICE_FREE_OCEAN;
		PetscPrintf(PETSC_COMM_SELF,"PISMPIK_INFO: [rank %d] killed isolated one-box-iceberg at i=%d,j=%d\n",grid.rank,i,j);
      }   
    }
  }

  //ierr =  vbed.end_access(); CHKERRQ(ierr);
  //ierr =  vMask.end_access(); CHKERRQ(ierr);
  //ierr =  vh.end_access(); CHKERRQ(ierr);
  ierr =  vH.end_access(); CHKERRQ(ierr);
  ierr =  vHnew.end_access(); CHKERRQ(ierr);
  ierr = vHnew.beginGhostComm(vH); CHKERRQ(ierr);
  ierr = vHnew.endGhostComm(vH); CHKERRQ(ierr);

/////////////////////////////////////////////////////////////////////////////
// looking for one-grid-cell partially filled grid cells, that have 4 neighbors of thickness H=0  
/////////////////////////////////////////////////////////////////////////////

  ierr = vH.copy_to(vHnew); CHKERRQ(ierr);
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vHnew.begin_access(); CHKERRQ(ierr);
  //ierr = vMask.begin_access(); CHKERRQ(ierr);
  //ierr = vh.begin_access(); CHKERRQ(ierr);
  //ierr = vbed.begin_access(); CHKERRQ(ierr);
  ierr = vHref.begin_access(); CHKERRQ(ierr);

  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      
	  // instead of updating surface elevation, counting here floating or icefree neighbors
//FIXME: unused var:      const PetscScalar hgrounded = vbed(i,j) + vH(i,j),
//FIXME: unused var:	  hfloating = currentSeaLevel + (1.0 - ice_rho/ocean_rho) * vH(i,j);
      
	  bool all_4neighbors_iceless = (vH(i+1,j)==0.0 && vH(i-1,j)==0.0 && vH(i,j+1)==0.0 && vH(i,j-1)==0.0);
	  // What about firstStepAfterInit?
      if (vHref(i,j)>0.0 && all_4neighbors_iceless) { 
	    vHref(i,j)=0.0;
	  	//vMask(i,j) = MASK_ICE_FREE_OCEAN;
		PetscPrintf(PETSC_COMM_SELF,"PISMPIK_INFO: [rank %d] killed lonely partially filled grid cell at i=%d,j=%d\n",grid.rank,i,j);
      }   
    }
  }

  ierr =  vbed.end_access(); CHKERRQ(ierr);
  ierr =  vMask.end_access(); CHKERRQ(ierr);
  ierr =  vh.end_access(); CHKERRQ(ierr);
  ierr =  vHref.end_access(); CHKERRQ(ierr);
  ierr =  vH.end_access(); CHKERRQ(ierr);
  ierr =  vHnew.end_access(); CHKERRQ(ierr);

  ierr = vHnew.beginGhostComm(vH); CHKERRQ(ierr);
  ierr = vHnew.endGhostComm(vH); CHKERRQ(ierr);
  ierr = vMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vMask.endGhostComm(); CHKERRQ(ierr);
  ierr = vh.beginGhostComm(); CHKERRQ(ierr);
  ierr = vh.endGhostComm(); CHKERRQ(ierr);
  ierr = vHref.beginGhostComm(); CHKERRQ(ierr);
  ierr = vHref.endGhostComm(); CHKERRQ(ierr);

  return 0;
}


