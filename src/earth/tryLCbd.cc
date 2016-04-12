// Copyright (C) 2007--2011, 2013, 2014, 2015, 2016 Ed Bueler and Constantine Khroulev
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
#include <gsl/gsl_math.h>       // M_PI
#include <cstdio>
#include <petscvec.h>
#include <petscdmda.h>
#include <petscdraw.h>

#include "base/util/pism_const.hh"
#include "base/util/Context.hh"
#include "base/util/PISMConfigInterface.hh"
#include "deformation.hh"
#include "base/util/pism_options.hh"

#include "base/util/petscwrappers/Viewer.hh"
#include "base/util/petscwrappers/PetscInitializer.hh"
#include "base/util/error_handling.hh"
#include "base/util/petscwrappers/Vec.hh"
#include "base/util/petscwrappers/DM.hh"

using namespace pism;

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm com = MPI_COMM_WORLD;  // won't be used except for rank
  int rank;

  petsc::Initializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;
  MPI_Comm_rank(com, &rank);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    Context::Ptr ctx = context_from_options(com, "tryLCbd");
    Config::Ptr config = ctx->config();

    pism::petsc::DM da2;
    petsc::Vec H, bed, Hstart, bedstart, uplift;

    bool
      include_elastic = false,
      do_uplift       = false;
    double H0         = 1000.0; // ice disc load thickness

    options::Keyword scenario("-scenario",
                              "chooses a scenario",
                              "1,2,3,4", "1");

    if (scenario == "1") {
      include_elastic = false;
      do_uplift       = false;
      H0              = 1000.0;
    } else if (scenario == "2") {
      include_elastic = true;
      do_uplift       = false;
      H0              = 1000.0;
    } else if (scenario == "3") {
      include_elastic = false;
      do_uplift       = true;
      H0              = 0.0;
    } else if (scenario == "4") {
      include_elastic = true;
      do_uplift       = true;
      H0              = 1000.0;
    } else {
      // this can't happen (options::Keyword validates its input), but still
      throw RuntimeError::formatted("invalid scenario %s", scenario->c_str());
    }

    const double R0 = 1000.0e3;          // ice disc load radius
    const double tfinalyears = 150.0e3;  // total run time

    options::Integer Mx("-Mx", "grid size in the X direction", 193);
    options::Integer My("-My", "grid size in the Y direction", 129);

    options::Real Lx("-Lx", "grid half-width in the X direction", 3000.0e3);
    options::Real Ly("-Ly", "grid half-width in the Y direction", 2000.0e3);

    const int Z = 2;
    const double dtyears = 100.0;

    if (rank == 0) { // only runs on proc 0; all sequential
      // allocate the variables needed before BedDeformLC can work:
      ierr = VecCreateSeq(PETSC_COMM_SELF, Mx*My, H.rawptr());
      PISM_CHK(ierr, "VecCreateSeq");

      ierr = VecDuplicate(H, Hstart.rawptr());
      PISM_CHK(ierr, "VecDuplicate");

      ierr = VecDuplicate(H, bedstart.rawptr());
      PISM_CHK(ierr, "VecDuplicate");

      ierr = VecDuplicate(H, uplift.rawptr());
      PISM_CHK(ierr, "VecDuplicate");

      // in order to show bed elevation as a picture, create a da
#if PETSC_VERSION_LT(3,5,0)
      ierr = DMDACreate2d(PETSC_COMM_SELF,
                          DMDA_BOUNDARY_PERIODIC, DMDA_BOUNDARY_PERIODIC,
                          DMDA_STENCIL_STAR,
                          My, Mx, PETSC_DECIDE, PETSC_DECIDE, 1, 0,
                          NULL, NULL, da2.rawptr());
      PISM_CHK(ierr, "DMDACreate2d");
#else
      ierr = DMDACreate2d(PETSC_COMM_SELF,
                          DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC,
                          DMDA_STENCIL_STAR,
                          My, Mx, PETSC_DECIDE, PETSC_DECIDE, 1, 0,
                          NULL, NULL, da2.rawptr());
      PISM_CHK(ierr, "DMDACreate2d");
#endif
      ierr = DMDASetUniformCoordinates(da2, -Ly, Ly, -Lx, Lx, 0, 0);
      PISM_CHK(ierr, "DMDASetUniformCoordinates");

      ierr = DMCreateGlobalVector(da2, bed.rawptr());
      PISM_CHK(ierr, "DMCreateGlobalVector");

      // create a bed viewer
      PetscDraw draw = NULL;
      const int
        windowx = 500,
        windowy = (int) (((float) windowx) * Ly / Lx);

      petsc::Viewer viewer;
      ierr = PetscViewerDrawOpen(PETSC_COMM_SELF, NULL, "bed elev (m)",
                                 PETSC_DECIDE, PETSC_DECIDE, windowy, windowx,
                                 viewer.rawptr());
      PISM_CHK(ierr, "PetscViewerDrawOpen");
      // following should be redundant, but may put up a title even under 2.3.3-p1:3 where
      // there is a no-titles bug
      ierr = PetscViewerDrawGetDraw(viewer,0,&draw);
      PISM_CHK(ierr, "PetscViewerDrawGetDraw");

      ierr = PetscDrawSetDoubleBuffer(draw);
      PISM_CHK(ierr, "PetscDrawSetDoubleBuffer");  // remove flicker while we are at it

      ierr = PetscDrawSetTitle(draw,"bed elev (m)");
      PISM_CHK(ierr, "PetscDrawSetTitle");

      // make disc load
      ierr = PetscPrintf(PETSC_COMM_SELF,"creating disc load\n");
      PISM_CHK(ierr, "PetscPrintf");

      // see "Results: Earth deformation only" section of Bueler et al "Fast computation ..."
      const double
        dx = (2.0*Lx)/((double) Mx - 1),
        dy = (2.0*Ly)/((double) My - 1);
      const int
        imid = (Mx-1)/2,
        jmid = (My-1)/2;

      {
        petsc::VecArray2D HH(H, Mx, My);
        for (int j=0; j<My; j++) {
          for (int i=0; i<Mx; i++) {
            const double r = sqrt(PetscSqr(dx * (i - imid)) + PetscSqr(dy * (j - jmid)));
            if (r < R0) {
              HH(i, j) = H0;
            } else {
              HH(i, j) = 0.0;
            }
          }
        }
      }

      ierr = VecSet(Hstart, 0.0);
      PISM_CHK(ierr, "VecSet");

      ierr = VecSet(bedstart, 0.0);
      PISM_CHK(ierr, "VecSet");

      const double peak_up = units::convert(ctx->unit_system(), 10, "mm year-1", "m second-1");  // 10 mm year-1
      // initialize uplift
      if (do_uplift == true) {
        petsc::VecArray2D upl(uplift, Mx, My);
        for (int j=0; j<My; j++) {
          for (int i=0; i<Mx; i++) {
            const double r = sqrt(PetscSqr(dx * (i - imid)) + PetscSqr(dy * (j - jmid)));
            if (r < 1.5 * R0) {
              upl(i, j) = peak_up * (cos(M_PI * (r / (1.5 * R0))) + 1.0) / 2.0;
            } else {
              upl(i, j) = 0.0;
            }
          }
        }
      } else {
        ierr = VecSet(uplift, 0.0);
        PISM_CHK(ierr, "VecSet");
      }

      ierr = PetscPrintf(PETSC_COMM_SELF,"setting BedDeformLC\n");
      PISM_CHK(ierr, "PetscPrintf");

      pism::bed::BedDeformLC bdlc(*config,
                                  include_elastic, Mx, My, dx, dy, Z,
                                  Hstart, bedstart, uplift, H, bed);

      ierr = PetscPrintf(PETSC_COMM_SELF,"allocating BedDeformLC\n");
      PISM_CHK(ierr, "PetscPrintf");

      ierr = PetscPrintf(PETSC_COMM_SELF,"initializing BedDeformLC from uplift map\n");
      PISM_CHK(ierr, "PetscPrintf");

      bdlc.uplift_init();

      ierr = PetscPrintf(PETSC_COMM_SELF,"stepping BedDeformLC\n");
      PISM_CHK(ierr, "PetscPrintf");

      const int KK = (int) (tfinalyears / dtyears);

      double b0old = 0.0;
      {
        petsc::VecArray2D b(bedstart, Mx, My);
        b0old = b(imid, jmid);
      }

      for (int k=0; k<KK; k++) {
        const double tyears = k*dtyears;

        bdlc.step(dtyears, tyears);

        ierr = VecView(bed,viewer);
        PISM_CHK(ierr, "VecView");

        double b0new = 0.0;
        {
          petsc::VecArray2D b(bed, Mx, My);
          b0new = b(imid, jmid);
        }

        const double dbdt0 = (b0new - b0old) / (dtyears);

        ierr = PetscPrintf(PETSC_COMM_SELF,
                  "   t=%8.0f (a)   b(0,0)=%11.5f (m)  dbdt(0,0)=%11.7f (m year-1)\n",
                           tyears, b0new, dbdt0);
        PISM_CHK(ierr, "PetscPrintf");

        char title[100];
        snprintf(title,100, "bed elev (m)  [t = %9.1f]", tyears);

        ierr = PetscDrawSetTitle(draw,title);
        PISM_CHK(ierr, "PetscDrawSetTitle");

        b0old = b0new;
      }

      ierr = PetscPrintf(PETSC_COMM_SELF,"\ndone\n");
      PISM_CHK(ierr, "PetscPrintf");
    }
  }
  catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
