// Copyright (C) 2007 Ed Bueler
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

static char help[] =
  "Simple testing program for Lingle & Clark bed deformation model.\n";

#include <cmath>
#include <cstdio>
#include <petscvec.h>
#include <petscda.h>
#include "../base/pism_const.hh"
#include "../base/materials.hh"
#include "../base/beddefLC.hh"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;  // won't be used except for rank,size
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);
  
  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    DeformableEarthType bedrock;
    //ThermoGlenIce       ice(com,ICE_PB);  // linker error???
    BedDeformLC bdlc;
    DA          da2;
    Vec         H, bed, Hstart, bedstart, uplift;
    
    PetscTruth  include_elastic = PETSC_FALSE,
                do_uplift = PETSC_FALSE;
    PetscScalar H0 = 1000.0;            // ice disc load thickness

    /* SCENARIOS:  run "./tryLCbd N" where N=1,2,3,4 as follows
       (1) dump ice disc on initially level, non-uplifting land, use only viscous 
           half-space model:
                 include_elastic = FALSE, do_uplift = FALSE, H0 = 1000.0
       (2) dump ice disc on initially level, non-uplifting land, use both viscous 
           half-space model and elastic model
                 include_elastic = TRUE, do_uplift = FALSE, H0 = 1000.0
       (3) never loaded, initially level, uplifting land, use only viscous 
           half-space model (because elastic model gives no additional when no load):
                 include_elastic = FALSE, do_uplift = TRUE, H0 = 0.0
       (4) dump ice disc on initially level, uplifting land, use both viscous 
           half-space model and elastic model:
                 include_elastic = TRUE, do_uplift = TRUE, H0 = 1000.0
    */

    if (argc >= 2) {
      switch (argv[1][0]) {
        case '1':
          include_elastic = PETSC_FALSE;  do_uplift = PETSC_FALSE;  H0 = 1000.0;
          break;
        case '2':
          include_elastic = PETSC_TRUE;  do_uplift = PETSC_FALSE;  H0 = 1000.0;
          break;
        case '3':
          include_elastic = PETSC_FALSE;  do_uplift = PETSC_TRUE;  H0 = 0.0;
          break;
        case '4':
          include_elastic = PETSC_TRUE;  do_uplift = PETSC_TRUE;  H0 = 1000.0;
          break;
        default:
          break; // accept default which is scenario 1
      }
    }
    const PetscScalar R0 = 1000.0e3;          // ice disc load radius
    const PetscScalar tfinalyears = 150.0e3;  // total run time

    const PetscInt    Mx = 97, 
                      My = 65;
    const PetscScalar Lx = 3000.0e3, 
                      Ly = 2000.0e3;
    const PetscInt    Z = 2;
    const PetscScalar dtyears = 100.0;
    
    if (rank == 0) { // only runs on proc 0; all sequential
      // allocate the variables needed before BedDeformLC can work:
      ierr = VecCreateSeq(PETSC_COMM_SELF, Mx*My, &H); CHKERRQ(ierr);
      ierr = VecDuplicate(H, &Hstart); CHKERRQ(ierr);
      ierr = VecDuplicate(H, &bedstart); CHKERRQ(ierr);
      ierr = VecDuplicate(H, &uplift); CHKERRQ(ierr);
      
      // in order to show bed elevation as a picture, create a da 
      ierr = DACreate2d(PETSC_COMM_SELF, DA_XYPERIODIC, DA_STENCIL_STAR,
                    My, Mx, PETSC_DECIDE, PETSC_DECIDE, 1, 0,
                    PETSC_NULL, PETSC_NULL, &da2); CHKERRQ(ierr);
      ierr = DASetUniformCoordinates(da2, -Ly, Ly, -Lx, Lx, 
                                 PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
      ierr = DACreateGlobalVector(da2, &bed); CHKERRQ(ierr);
      
      // create a bed viewer
      PetscViewer viewer;
      PetscDraw   draw;
      const PetscInt  windowx = 500,
                      windowy = (PetscInt) (((float) windowx) * Ly / Lx);
      ierr = PetscViewerDrawOpen(PETSC_COMM_SELF, PETSC_NULL, "bed elevation will go here",
           PETSC_DECIDE, PETSC_DECIDE, windowy, windowx, &viewer);  CHKERRQ(ierr);
      // following should be redundant, but may put up a title even under 2.3.3-p1:3 where
      // there is a no-titles bug
      ierr = PetscViewerDrawGetDraw(viewer,0,&draw); CHKERRQ(ierr);
      ierr = PetscDrawSetDoubleBuffer(draw); CHKERRQ(ierr);  // remove flicker while we are at it
      ierr = PetscDrawSetTitle(draw,"bed elevation will go here"); CHKERRQ(ierr);

      // make disc load
      ierr = PetscPrintf(PETSC_COMM_SELF,"creating disc load\n"); CHKERRQ(ierr);
      // see "Results: Earth deformation only" section of Bueler et al "Fast computation ..."
      const PetscScalar dx = (2.0*Lx)/((PetscScalar) Mx - 1), 
                        dy = (2.0*Ly)/((PetscScalar) My - 1);
      const PetscInt    imid = (Mx-1)/2, jmid = (My-1)/2;
      PetscScalar **HH;
      ierr = VecGetArray2d(H, Mx, My, 0, 0, &HH); CHKERRQ(ierr);
      for (PetscInt i=0; i<Mx; i++) {
        for (PetscInt j=0; j<My; j++) {
          const PetscScalar r = sqrt( PetscSqr(dx * (i - imid)) + PetscSqr(dy * (j - jmid)) );
          if (r < R0) {
            HH[i][j] = H0;
          } else {
            HH[i][j] = 0.0;
          }
        }
      }
      ierr = VecRestoreArray2d(H, Mx, My, 0, 0, &HH); CHKERRQ(ierr);
      ierr = VecSet(Hstart, 0.0); CHKERRQ(ierr);    // load was zero up till t=0
      ierr = VecSet(bedstart, 0.0); CHKERRQ(ierr);       // initially flat bed
      
      // initialize uplift
      if (do_uplift == PETSC_TRUE) {
        PetscScalar **upl;
        ierr = VecGetArray2d(uplift, Mx, My, 0, 0, &upl); CHKERRQ(ierr);
        for (PetscInt i=0; i<Mx; i++) {
          for (PetscInt j=0; j<My; j++) {
            const PetscScalar peak_up = 0.010 / secpera;  // 10 mm/a
            const PetscScalar r = sqrt( PetscSqr(dx * (i - imid)) + PetscSqr(dy * (j - jmid)) );
            if (r < 1.5 * R0) {
              upl[i][j] = peak_up * (cos(pi * (r / (1.5 * R0))) + 1.0) / 2.0; 
            } else {
              upl[i][j] = 0.0;
            }
          }
        }
        ierr = VecRestoreArray2d(uplift, Mx, My, 0, 0, &upl); CHKERRQ(ierr);
      } else {
        ierr = VecSet(uplift, 0.0); CHKERRQ(ierr);
      }

      ierr = PetscPrintf(PETSC_COMM_SELF,"setting BedDeformLC\n"); CHKERRQ(ierr);
      ierr = bdlc.settings(include_elastic,Mx,My,dx,dy,Z, 
//                           ice.rho, bedrock.rho, bedrock.eta, bedrock.D,
                           910.0, bedrock.rho, bedrock.eta, bedrock.D,
                           &Hstart, &bedstart, &uplift, &H, &bed); CHKERRQ(ierr);

      ierr = PetscPrintf(PETSC_COMM_SELF,"allocating BedDeformLC\n"); CHKERRQ(ierr);
      ierr = bdlc.alloc(); CHKERRQ(ierr);
      
      ierr = PetscPrintf(PETSC_COMM_SELF,"initializing BedDeformLC from uplift map\n"); CHKERRQ(ierr);
      ierr = bdlc.uplift_init(); CHKERRQ(ierr);
      
      ierr = PetscPrintf(PETSC_COMM_SELF,"stepping BedDeformLC\n"); CHKERRQ(ierr);
      const PetscInt     KK = (PetscInt) (tfinalyears / dtyears);
      PetscScalar **b;
      ierr = VecGetArray2d(bedstart, Mx, My, 0, 0, &b); CHKERRQ(ierr);
      PetscScalar b0old = b[imid][jmid];
      ierr = VecRestoreArray2d(bedstart, Mx, My, 0, 0, &b); CHKERRQ(ierr);
      for (PetscInt k=0; k<KK; k++) {
        const PetscScalar tyears = k*dtyears;
        ierr = bdlc.step(dtyears, tyears); CHKERRQ(ierr);
        ierr = VecView(bed,viewer); CHKERRQ(ierr);
        char title[100];
        PetscScalar **b;
        ierr = VecGetArray2d(bed, Mx, My, 0, 0, &b); CHKERRQ(ierr);
        const PetscScalar b0new = b[imid][jmid];
        ierr = VecRestoreArray2d(bed, Mx, My, 0, 0, &b); CHKERRQ(ierr);
        const PetscScalar dbdt0 = (b0new - b0old) / (dtyears);
        sprintf(title, "bed elevation (m)   [t = %9.1f]", tyears);
        ierr = PetscPrintf(PETSC_COMM_SELF,
                  "   t=%8.0f (a)   b(0,0)=%11.5f (m)  dbdt(0,0)=%11.7f (m/a)\n",
                  tyears, b0new, dbdt0); CHKERRQ(ierr);
        ierr = PetscDrawSetTitle(draw,title); CHKERRQ(ierr);
        b0old = b0new;
      }
      ierr = PetscPrintf(PETSC_COMM_SELF,"\ndone\n"); CHKERRQ(ierr);

      ierr = VecDestroy(H); CHKERRQ(ierr);
      ierr = VecDestroy(bed); CHKERRQ(ierr);
      ierr = VecDestroy(Hstart); CHKERRQ(ierr);
      ierr = VecDestroy(bedstart); CHKERRQ(ierr);
      ierr = VecDestroy(uplift); CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer); CHKERRQ(ierr);
      ierr = DADestroy(da2); CHKERRQ(ierr);
    }
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
