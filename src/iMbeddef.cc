// Copyright (C) 2004-2007 Ed Bueler
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

#include <petscda.h>
#include <cmath>
#include "iceModel.hh"
#include "extrasGSL.hh"
#include "beddefLC.hh"

/* the following is from the PETSc FAQ page:

How do I collect all the values from a parallel PETSc vector into a vector on the zeroth processor?

    * Create the scatter context that will do the communication
          o VecScatterCreateToZero(v,&ctx,&w);
    * Actually do the communication; this can be done repeatedly as needed
          o VecScatterBegin(ctx,v,w,INSERT_VALUES,SCATTER_FORWARD);
          o VecScatterEnd(ctx,v,w,INSERT_VALUES,SCATTER_FORWARD);
    * Remember to free the scatter context when no longer needed
          o VecScatterDestroy(ctx);

Note that this simply concatenates in the parallel ordering of the vector.  If you are using a vector from DACreateGlobalVector() you likely want to first call DAGlobalToNaturalBegin/End() to scatter the original vector into the natural ordering in a new global vector before calling VecScatterBegin/End() to scatter the natural vector onto process 0.

*/

PetscErrorCode IceModel::createScatterToProcZero(Vec& samplep0) {
  PetscErrorCode ierr;

  if (top0ctx_created == PETSC_FALSE) {
    // note we want a global Vec but reordered in the natural ordering so when it is
    // scattered to proc zero it is not all messed up; see above
    ierr = DACreateNaturalVector(grid.da2, &g2natural); CHKERRQ(ierr);
    // next get context *and* allocate samplep0 (on proc zero only, naturally)
    ierr = VecScatterCreateToZero(g2natural, &top0ctx, &samplep0); CHKERRQ(ierr);
    top0ctx_created = PETSC_TRUE;
  }
  return 0;
}


PetscErrorCode IceModel::destroyScatterToProcZero() {
  PetscErrorCode ierr;

  if (top0ctx_created == PETSC_TRUE) {
    ierr = VecDestroy(g2natural); CHKERRQ(ierr);
    ierr = VecScatterDestroy(top0ctx); CHKERRQ(ierr);
    top0ctx_created = PETSC_FALSE;
  }
  return 0;
}


PetscErrorCode IceModel::putLocalOnProcZero(Vec& vlocal, Vec& onp0) {
  PetscErrorCode ierr;

  // scatter local Vec to proc zero: from a global Vec in the global ordering to 
  //    a global Vec in the natural ordering and then to a Vec on proc zero (i.e. empty on other procs)
  // requires g2, g2natural, and top0ctx to all be set up properly
  ierr = DALocalToGlobal(grid.da2,vlocal,INSERT_VALUES,g2); CHKERRQ(ierr);
  ierr = DAGlobalToNaturalBegin(grid.da2,g2,INSERT_VALUES,g2natural); CHKERRQ(ierr);
  ierr = DAGlobalToNaturalEnd(grid.da2,g2,INSERT_VALUES,g2natural); CHKERRQ(ierr);
  ierr = VecScatterBegin(top0ctx, g2natural,onp0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecScatterEnd(top0ctx, g2natural,onp0,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::getLocalFromProcZero(Vec& onp0, Vec& vlocal) {
  PetscErrorCode ierr;

  // undo scatter to proc zero: put onp0 back into vlocal
  ierr = VecScatterBegin(top0ctx, onp0,g2natural,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = VecScatterEnd(top0ctx, onp0,g2natural,INSERT_VALUES,SCATTER_REVERSE); CHKERRQ(ierr);
  ierr = DANaturalToGlobalBegin(grid.da2,g2natural,INSERT_VALUES,g2); CHKERRQ(ierr);
  ierr = DANaturalToGlobalEnd(grid.da2,g2natural,INSERT_VALUES,g2); CHKERRQ(ierr);
  ierr = DAGlobalToLocalBegin(grid.da2,g2,INSERT_VALUES,vlocal); CHKERRQ(ierr);
  ierr = DAGlobalToLocalEnd(grid.da2,g2,INSERT_VALUES,vlocal); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::bedDefSetup() {
  PetscErrorCode  ierr;
  
  if (doBedDef == PETSC_TRUE) {
    ierr = VecDuplicate(vH,&vHlast); CHKERRQ(ierr);
    ierr = VecCopy(vH,vHlast); CHKERRQ(ierr);
    ierr = VecDuplicate(vbed,&vbedlast); CHKERRQ(ierr);
    ierr = VecCopy(vbed,vbedlast); CHKERRQ(ierr);
    lastBedDefUpdateYear = grid.p->year;
    
    if (doBedIso == PETSC_TRUE) {
      ierr = verbPrintf(2, grid.com,"using pointwise isostasy Earth deformation model\n");
           CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid.com,
          "using Lingle & Clark Earth deformation model; setting up bed def variables ..."); 
          CHKERRQ(ierr);
#if (WITH_FFTW==0)
      ierr = PetscPrintf(grid.com,
          "\n  WARNING: compiled without FFTW.  -bed_def_lc (-bed_def) will not work.\n"); 
          CHKERRQ(ierr);
#endif
      // create scatter context and allocate the Vecs for BedDeformLC
      ierr = createScatterToProcZero(Hp0); CHKERRQ(ierr);
      ierr = VecDuplicate(Hp0,&bedp0);
      ierr = VecDuplicate(Hp0,&Hstartp0);
      ierr = VecDuplicate(Hp0,&bedstartp0);
      ierr = VecDuplicate(Hp0,&upliftp0);
      // store initial H, bed, and uplift values on proc zero
      // note load is  - rho_i g (H-Hstart)
      // and  bed = bedstart + U_elastic + (U_plate - U_prebent)
      ierr = putLocalOnProcZero(vH,Hstartp0); CHKERRQ(ierr);
      ierr = putLocalOnProcZero(vbed,bedstartp0); CHKERRQ(ierr);
      ierr = putLocalOnProcZero(vuplift,upliftp0); CHKERRQ(ierr);

      if (grid.rank == 0) {
        ierr = bdLC.settings(PETSC_FALSE, // turn off elastic model for now
                       grid.p->Mx,grid.p->My,grid.p->dx,grid.p->dy,
                       2,                 // use Z = 2 for now
                       ice.rho, bedrock.rho, bedrock.eta, bedrock.D,
                       &Hstartp0, &bedstartp0, &upliftp0, &Hp0, &bedp0); CHKERRQ(ierr);
        ierr = bdLC.alloc(); CHKERRQ(ierr);
        // *always* initialize plate from uplift, but make sure uplift (=dbed/dt)
        // is zero if it is not actually available from data
        ierr = bdLC.uplift_init(); CHKERRQ(ierr);
      } // if (grid.rank == 0)
      ierr = verbPrintf(2,grid.com," done\n"); CHKERRQ(ierr);
    }
  }
  return 0;
}


PetscErrorCode IceModel::bedDefCleanup() {
  PetscErrorCode  ierr;

  if (doBedDef == PETSC_TRUE) {
    ierr = VecDestroy(vHlast); CHKERRQ(ierr);
    ierr = VecDestroy(vbedlast); CHKERRQ(ierr);

    if (doBedIso == PETSC_FALSE) {
      ierr = VecDestroy(Hp0); CHKERRQ(ierr);
      ierr = VecDestroy(bedp0); CHKERRQ(ierr);
      ierr = VecDestroy(Hstartp0); CHKERRQ(ierr);
      ierr = VecDestroy(bedstartp0); CHKERRQ(ierr);
      ierr = VecDestroy(upliftp0); CHKERRQ(ierr);
      ierr = destroyScatterToProcZero(); CHKERRQ(ierr);
      // note destructor for bdLC (i.e. ~BedDeformLC()) will free bdLC's storage
    }
  }
  return 0;
}


PetscErrorCode IceModel::bedDefStepIfNeeded() {
  PetscErrorCode  ierr;

  // This is a front end to the bed deformation update system.  It updates
  // no more often than bedDefIntervalYears.
  if (doBedDef == PETSC_TRUE) {
    //ierr = PetscPrintf(grid.com, " [lastBedDefUpdateYear = %10.3f]\n",
    //           lastBedDefUpdateYear); CHKERRQ(ierr);
    // If the bed elevations are not expired, exit cleanly.
    const PetscScalar dtBedDefYears = grid.p->year - lastBedDefUpdateYear;
    if (dtBedDefYears >= bedDefIntervalYears) {
      if (doBedIso == PETSC_TRUE) { // pointwise isostasy model: in parallel
        ierr = bed_def_step_iso(); CHKERRQ(ierr);
      } else { // L&C model: only on proc zero
        ierr = putLocalOnProcZero(vH,Hp0); CHKERRQ(ierr);
        ierr = putLocalOnProcZero(vbed,bedp0); CHKERRQ(ierr);
        if (grid.rank == 0) {  // only processor zero does the step
          ierr = bdLC.step(dtBedDefYears, grid.p->year - startYear); CHKERRQ(ierr);
        }
        ierr = getLocalFromProcZero(bedp0, vbed); CHKERRQ(ierr);
        ierr = DALocalToLocalBegin(grid.da2,vbed,INSERT_VALUES,vbed); CHKERRQ(ierr);
        ierr = DALocalToLocalEnd(grid.da2,vbed,INSERT_VALUES,vbed); CHKERRQ(ierr);
      }
      // update uplift: uplift = (bed - bedlast) / dt
      ierr = VecWAXPY(vuplift,-1.0,vbedlast,vbed); CHKERRQ(ierr);  
      ierr = VecScale(vuplift, 1.0 / (dtBedDefYears * secpera)); CHKERRQ(ierr); 
      // copy current values of H, bed in prep for next step
      ierr = VecCopy(vH,vHlast); CHKERRQ(ierr);
      ierr = VecCopy(vbed,vbedlast); CHKERRQ(ierr);
      ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);
      lastBedDefUpdateYear = grid.p->year;
      ierr = PetscPrintf(grid.com, "b"); CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(grid.com, "$"); CHKERRQ(ierr);
    }
  } else {
    ierr = PetscPrintf(grid.com, "$"); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode IceModel::bed_def_step_iso() {
  PetscErrorCode ierr;
  Vec vHdiff = vWork2d[0];

  const PetscScalar  f = ice.rho / bedrock.rho;
  ierr = VecWAXPY(vHdiff,-1.0,vHlast,vH); CHKERRQ(ierr);  // Hdiff = H - Hlast
  ierr = VecWAXPY(vbed,-f,vHdiff,vbedlast); CHKERRQ(ierr);  // bed = bedlast - f (Hdiff)
  return 0;
}

