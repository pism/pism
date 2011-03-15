// Copyright (C) 2004--2011 Torsten Albrecht, Jed Brown, Ed Bueler and Constantine Khroulev
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

/// !!!!!! 
/// Icebergs, i.e. floating regions that are not attached to at least one dragging box, cause unrealistically large velocities that numerically affect the ice velocities everywhere else, lead to extremely small time steps and can eventually cause a KSP-ERROR. The purpose of the following routines is to identify and eliminate these icebergs.
/// !!!!!! 

//! This routine comes from PISM-PIK.
/*!
 The aim of this routine is to find floating regions that *might* be icebergs. If these regions actually *are* icebergs is checked in extIdentifyNotAnIceBerg.
*/


PetscErrorCode IceModel::FindIceBergCandidates() {
  PetscErrorCode ierr;

  ierr = verbPrintf(4,grid.com,"######### FindIceBergCandidates is called \n");    CHKERRQ(ierr);

  const PetscInt Mx = grid.Mx, My = grid.My;
  ierr = vH.begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);
  ierr = vIcebergMask.begin_access(); CHKERRQ(ierr);

  //const PetscScalar   a = grid.dx * grid.dy * 1e-3 * 1e-3; // area unit (km^2)
	
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {

      //cut of border of computational domain
      if(i <= 0 || i >= Mx-1 || j <= 0 || j >= My-1) {
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

  ierr = vMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vMask.endGhostComm(); CHKERRQ(ierr);
  ierr = vIcebergMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vIcebergMask.endGhostComm(); CHKERRQ(ierr);

  //////////////////////////////////////////////////////////////////////////////////
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

  //////////////////////////////////////////////////////////////////////////////////
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




PetscErrorCode IceModel::IdentifyNotAnIceBerg() {
  PetscErrorCode ierr;
  
  ierr = verbPrintf(4,grid.com,"######### IdentifyNotAnIceBerg is called \n");    CHKERRQ(ierr);
  
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
  ierr = verbPrintf(4,grid.com,"!!! %d loop(s) were needed to identify whether there are icebergs \n",loopcount);    CHKERRQ(ierr);

  return 0;
}


/*!
 We have distinguished icebergs from attached floating regions in extIdentifyNotAnIceBerg. Now we actually eliminate the former (meaning we set the ice thickness to zero and mark the boxes as icefree-ocean) and leave the latter as they are.
*/

PetscErrorCode IceModel::killIceBergs() {
  PetscErrorCode ierr;
  
  ierr = verbPrintf(3,grid.com,"######### killIceBergs is called \n");    CHKERRQ(ierr);
  
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
