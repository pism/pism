// Copyright (C) 2007--2011, 2013, 2014 Ed Bueler and Constantine Khroulev
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

// get this message with 'tryLCbd -help':
static char help[] =
  "Simple testing program for Lingle & Clark bed deformation model.\n"
  "  Runs go for 150,000 years on 63.5km grid with 100a time steps and Z=2 in L&C model.\n"
  "  SCENARIOS:  run './tryLCbd N' where N=1,2,3,4 as follows\n"
  "     (1) dump ice disc on initially level, non-uplifting land, use only viscous \n"
  "         half-space model:\n"
  "               include_elastic = FALSE, do_uplift = FALSE, H0 = 1000.0\n"
  "         center depth b(0,0) should eventually equilibriate to near\n"
  "         -1000 * (910/3300) = -275.76 m\n"
  "     (2) dump ice disc on initially level, non-uplifting land, use both viscous \n"
  "         half-space model and elastic model\n"
  "               include_elastic = TRUE, do_uplift = FALSE, H0 = 1000.0\n"
  "     (3) never loaded, initially level, uplifting land, use only viscous \n"
  "         half-space model (because elastic model gives no additional when no load):\n"
  "               include_elastic = FALSE, do_uplift = TRUE, H0 = 0.0\n"
  "     (4) dump ice disc on initially level, uplifting land, use both viscous \n"
  "         half-space model and elastic model:\n"
  "               include_elastic = TRUE, do_uplift = TRUE, H0 = 1000.0\n\n";


#include <cmath>
#include <cstdio>
#include <petscvec.h>
#include <petscdmda.h>
#include "pism_const.hh"
#include "PISMConfig.hh"
#include "deformation.hh"
#include "pism_options.hh"

#include "PetscInitializer.hh"
#include "error_handling.hh"

using namespace pism;

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm com = MPI_COMM_WORLD;  // won't be used except for rank
  int rank;

  PetscInitializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  
  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    UnitSystem unit_system;
    Config config(com, "pism_config", unit_system),
      overrides(com, "pism_overrides", unit_system);
    ierr = init_config(com, config, overrides); CHKERRQ(ierr);

    BedDeformLC bdlc;
    DM          da2;
    Vec         H, bed, Hstart, bedstart, uplift;
    
    PetscBool  include_elastic = PETSC_FALSE,
                do_uplift = PETSC_FALSE;
    double H0 = 1000.0;            // ice disc load thickness

    if (argc >= 2) {
      // FIXME:  should use PETSC-style options
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
    const double R0 = 1000.0e3;          // ice disc load radius
    const double tfinalyears = 150.0e3;  // total run time

    // FIXME: should accept options here
    const int    Mx = 193, 
                      My = 129;

    const double Lx = 3000.0e3, 
                      Ly = 2000.0e3;
    const int    Z = 2;
    const double dtyears = 100.0;
    
    if (rank == 0) { // only runs on proc 0; all sequential
      // allocate the variables needed before BedDeformLC can work:
      ierr = VecCreateSeq(PETSC_COMM_SELF, Mx*My, &H); CHKERRQ(ierr);
      ierr = VecDuplicate(H, &Hstart); CHKERRQ(ierr);
      ierr = VecDuplicate(H, &bedstart); CHKERRQ(ierr);
      ierr = VecDuplicate(H, &uplift); CHKERRQ(ierr);

      // in order to show bed elevation as a picture, create a da 
#if PETSC_VERSION_LT(3,5,0)
      ierr = DMDACreate2d(PETSC_COMM_SELF,
                          DMDA_BOUNDARY_PERIODIC, DMDA_BOUNDARY_PERIODIC,
                          DMDA_STENCIL_STAR,
                          My, Mx, PETSC_DECIDE, PETSC_DECIDE, 1, 0,
                          NULL, NULL, &da2); CHKERRQ(ierr);
#else
      ierr = DMDACreate2d(PETSC_COMM_SELF,
                          DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                          DMDA_STENCIL_STAR,
                          My, Mx, PETSC_DECIDE, PETSC_DECIDE, 1, 0,
                          NULL, NULL, &da2); CHKERRQ(ierr);
#endif
      ierr = DMDASetUniformCoordinates(da2, -Ly, Ly, -Lx, Lx, 0, 0); CHKERRQ(ierr);
      ierr = DMCreateGlobalVector(da2, &bed); CHKERRQ(ierr);

      // create a bed viewer
      PetscViewer viewer;
      PetscDraw   draw;
      const int  windowx = 500,
                      windowy = (int) (((float) windowx) * Ly / Lx);
      ierr = PetscViewerDrawOpen(PETSC_COMM_SELF, NULL, "bed elev (m)",
           PETSC_DECIDE, PETSC_DECIDE, windowy, windowx, &viewer);  CHKERRQ(ierr);
      // following should be redundant, but may put up a title even under 2.3.3-p1:3 where
      // there is a no-titles bug
      ierr = PetscViewerDrawGetDraw(viewer,0,&draw); CHKERRQ(ierr);
      ierr = PetscDrawSetDoubleBuffer(draw); CHKERRQ(ierr);  // remove flicker while we are at it
      ierr = PetscDrawSetTitle(draw,"bed elev (m)"); CHKERRQ(ierr);

      // make disc load
      ierr = PetscPrintf(PETSC_COMM_SELF,"creating disc load\n"); CHKERRQ(ierr);
      // see "Results: Earth deformation only" section of Bueler et al "Fast computation ..."
      const double dx = (2.0*Lx)/((double) Mx - 1), 
                        dy = (2.0*Ly)/((double) My - 1);
      const int    imid = (Mx-1)/2, jmid = (My-1)/2;
      double **HH;
      ierr = VecGetArray2d(H, Mx, My, 0, 0, &HH); CHKERRQ(ierr);
      for (int i=0; i<Mx; i++) {
        for (int j=0; j<My; j++) {
          const double r = sqrt(PetscSqr(dx * (i - imid)) + PetscSqr(dy * (j - jmid)));
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
      
      const double peak_up = unit_system.convert(10, "mm/year", "m/s");  // 10 mm/year
      // initialize uplift
      if (do_uplift == PETSC_TRUE) {
        double **upl;
        ierr = VecGetArray2d(uplift, Mx, My, 0, 0, &upl); CHKERRQ(ierr);
        for (int i=0; i<Mx; i++) {
          for (int j=0; j<My; j++) {
            const double r = sqrt(PetscSqr(dx * (i - imid)) + PetscSqr(dy * (j - jmid)));
            if (r < 1.5 * R0) {
              upl[i][j] = peak_up * (cos(M_PI * (r / (1.5 * R0))) + 1.0) / 2.0; 
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
      ierr = bdlc.settings(config,
                           include_elastic, Mx, My, dx, dy, Z,
                           &Hstart, &bedstart, &uplift, &H, &bed); CHKERRQ(ierr);

      ierr = PetscPrintf(PETSC_COMM_SELF,"allocating BedDeformLC\n"); CHKERRQ(ierr);
      ierr = bdlc.alloc(); CHKERRQ(ierr);
      
      ierr = PetscPrintf(PETSC_COMM_SELF,"initializing BedDeformLC from uplift map\n"); CHKERRQ(ierr);
      ierr = bdlc.init(); CHKERRQ(ierr);
      ierr = bdlc.uplift_init(); CHKERRQ(ierr);
      
      ierr = PetscPrintf(PETSC_COMM_SELF,"stepping BedDeformLC\n"); CHKERRQ(ierr);
      const int     KK = (int) (tfinalyears / dtyears);
      double **b;
      ierr = VecGetArray2d(bedstart, Mx, My, 0, 0, &b); CHKERRQ(ierr);
      double b0old = b[imid][jmid];
      ierr = VecRestoreArray2d(bedstart, Mx, My, 0, 0, &b); CHKERRQ(ierr);
      for (int k=0; k<KK; k++) {
        const double tyears = k*dtyears;
        ierr = bdlc.step(dtyears, tyears); CHKERRQ(ierr);
        ierr = VecView(bed,viewer); CHKERRQ(ierr);
        ierr = VecGetArray2d(bed, Mx, My, 0, 0, &b); CHKERRQ(ierr);
        const double b0new = b[imid][jmid];
        ierr = VecRestoreArray2d(bed, Mx, My, 0, 0, &b); CHKERRQ(ierr);
        const double dbdt0 = (b0new - b0old) / (dtyears);

        ierr = PetscPrintf(PETSC_COMM_SELF,
                  "   t=%8.0f (a)   b(0,0)=%11.5f (m)  dbdt(0,0)=%11.7f (m/year)\n",
                  tyears, b0new, dbdt0); CHKERRQ(ierr);

        char title[100];
        snprintf(title,100, "bed elev (m)  [t = %9.1f]", tyears);
        ierr = PetscDrawSetTitle(draw,title); CHKERRQ(ierr);
        b0old = b0new;
      }
      ierr = PetscPrintf(PETSC_COMM_SELF,"\ndone\n"); CHKERRQ(ierr);

      ierr = VecDestroy(&H); CHKERRQ(ierr);
      ierr = VecDestroy(&bed); CHKERRQ(ierr);
      ierr = VecDestroy(&Hstart); CHKERRQ(ierr);
      ierr = VecDestroy(&bedstart); CHKERRQ(ierr);
      ierr = VecDestroy(&uplift); CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
      ierr = DMDestroy(&da2); CHKERRQ(ierr);
    }
  }
  catch (...) {
    handle_fatal_errors(com);
  }

  return 0;
}
