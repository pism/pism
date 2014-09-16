// Copyright (C) 2010, 2011, 2012, 2013, 2014 Ed Bueler
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

static char help[] = "\nBEDROUGH_TEST\n"
  "  Simple testing program for Schoof (2003)-type bed smoothing and roughness-\n"
  "  parameterization schemes.  Allows comparison of computed theta to result\n"
  "  from Matlab/Octave code exampletheta.m in src/base/bedroughplay.  Also\n"
  "  used in PISM software (regression) test.\n\n";

#include "PISMConfig.hh"
#include <cmath>
#include <cstdio>
#include "pism_options.hh"
#include "IceGrid.hh"
#include "iceModelVec.hh"
#include "PISMBedSmoother.hh"

using namespace pism;

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, help); CHKERRQ(ierr);

  MPI_Comm com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    UnitSystem unit_system;
    Config config(com, "pism_config", unit_system),
      overrides(com, "pism_overrides", unit_system);
    ierr = init_config(com, config, overrides); CHKERRQ(ierr);

    IceGrid grid(com, config);
    grid.Mx = 81;
    grid.My = 81;
    grid.Lx = 1200e3;
    grid.Ly = grid.Lx;
    grid.periodicity = NOT_PERIODIC;
    grid.compute_nprocs();
    grid.compute_ownership_ranges();
    ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
    ierr = grid.allocate(); CHKERRQ(ierr);

    PetscPrintf(grid.com,"BedSmoother TEST\n");

    bool show;
    ierr = OptionsIsSet("-show", show); CHKERRQ(ierr);

    IceModelVec2S topg, usurf, theta;
    ierr = topg.create(grid, "topg", WITH_GHOSTS, 1); CHKERRQ(ierr);
    ierr = topg.set_attrs(
                          "trybedrough_tool", "original topography",
                          "m", "bedrock_altitude"); CHKERRQ(ierr);
    ierr = usurf.create(grid, "usurf", WITH_GHOSTS, 1); CHKERRQ(ierr);
    ierr = usurf.set_attrs(
                           "trybedrough_tool", "ice surface elevation",
                           "m", "surface_altitude"); CHKERRQ(ierr);
    ierr = theta.create(grid, "theta", WITH_GHOSTS, 1); CHKERRQ(ierr);
    ierr = theta.set_attrs(
                           "trybedrough_tool",
                           "coefficient theta in Schoof (2003) bed roughness parameterization",
                           "", ""); CHKERRQ(ierr);

    // put in bed elevations, a la this Matlab:
    //    topg0 = 400 * sin(2 * pi * xx / 600e3) + ...
    //            100 * sin(2 * pi * (xx + 1.5 * yy) / 40e3);
    ierr = topg.begin_access(); CHKERRQ(ierr);
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      topg(i,j) = 400.0 * sin(2.0 * M_PI * grid.x[i] / 600.0e3) +
        100.0 * sin(2.0 * M_PI * (grid.x[i] + 1.5 * grid.y[j]) / 40.0e3);
    }
    ierr = topg.end_access(); CHKERRQ(ierr);

    ierr = usurf.set(1000.0); CHKERRQ(ierr);  // compute theta for this constant thk

    // actually use the smoother/bed-roughness-parameterizer
    config.set_double("Glen_exponent", 3.0);
    config.set_double("bed_smoother_range", 50.0e3);
    BedSmoother smoother(grid, config, 1);
    ierr = smoother.preprocess_bed(topg); CHKERRQ(ierr);
    int Nx,Ny;
    ierr = smoother.get_smoothing_domain(Nx,Ny); CHKERRQ(ierr);
    PetscPrintf(grid.com,"  smoothing domain:  Nx = %d, Ny = %d\n",Nx,Ny);
    ierr = smoother.get_theta(usurf, &theta); CHKERRQ(ierr);

    if (show) {
      const int  window = 400;
      ierr = topg.view(window);  CHKERRQ(ierr);
      ierr = smoother.topgsmooth.view(window);  CHKERRQ(ierr);
      ierr = theta.view(window);  CHKERRQ(ierr);
      printf("[showing topg, smoother.topgsmooth, theta in X windows for 10 seconds ...]\n");
      ierr = PetscSleep(10); CHKERRQ(ierr);
    }

    double topg_min, topg_max, topgs_min, topgs_max, theta_min, theta_max;
    ierr = topg.min(topg_min); CHKERRQ(ierr);
    ierr = topg.max(topg_max); CHKERRQ(ierr);
    ierr = smoother.topgsmooth.min(topgs_min); CHKERRQ(ierr);
    ierr = smoother.topgsmooth.max(topgs_max); CHKERRQ(ierr);
    ierr = theta.min(theta_min); CHKERRQ(ierr);
    ierr = theta.max(theta_max); CHKERRQ(ierr);
    PetscPrintf(grid.com,
                "  original bed    :  min elev = %12.6f m,  max elev = %12.6f m\n",
                topg_min, topg_max);
    PetscPrintf(grid.com,
                "  smoothed bed    :  min elev = %12.6f m,  max elev = %12.6f m\n",
                topgs_min, topgs_max);
    PetscPrintf(grid.com,
                "  Schoof's theta  :  min      = %12.9f,    max      = %12.9f\n",
                theta_min, theta_max);

  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
