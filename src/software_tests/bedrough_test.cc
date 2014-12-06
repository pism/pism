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

#include "PetscInitializer.hh"
#include "error_handling.hh"

using namespace pism;

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;
  MPI_Comm com = MPI_COMM_WORLD;

  PetscInitializer petsc(argc, argv, help);

  com = PETSC_COMM_WORLD;

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  try {
    UnitSystem unit_system;
    Config config(com, "pism_config", unit_system),
      overrides(com, "pism_overrides", unit_system);
    ierr = init_config(com, config, overrides); CHKERRQ(ierr);

    IceGrid grid(com, config);
    grid.set_Mx(81);
    grid.set_My(81);
    grid.Lx = 1200e3;
    grid.Ly = grid.Lx;
    grid.periodicity = NOT_PERIODIC;
    ierr = grid.allocate(); CHKERRQ(ierr);

    PetscPrintf(grid.com,"BedSmoother TEST\n");

    bool show;
    ierr = OptionsIsSet("-show", show); CHKERRQ(ierr);

    IceModelVec2S topg, usurf, theta;
    topg.create(grid, "topg", WITH_GHOSTS, 1);
    topg.set_attrs("trybedrough_tool", "original topography",
                   "m", "bedrock_altitude");
    usurf.create(grid, "usurf", WITH_GHOSTS, 1);
    usurf.set_attrs("trybedrough_tool", "ice surface elevation",
                    "m", "surface_altitude");
    theta.create(grid, "theta", WITH_GHOSTS, 1);
    theta.set_attrs("trybedrough_tool",
                    "coefficient theta in Schoof (2003) bed roughness parameterization",
                    "", "");

    // put in bed elevations, a la this Matlab:
    //    topg0 = 400 * sin(2 * pi * xx / 600e3) + ...
    //            100 * sin(2 * pi * (xx + 1.5 * yy) / 40e3);
    IceModelVec::AccessList list(topg);
    for (Points p(grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      topg(i,j) = 400.0 * sin(2.0 * M_PI * grid.x(i) / 600.0e3) +
        100.0 * sin(2.0 * M_PI * (grid.x(i) + 1.5 * grid.y(j)) / 40.0e3);
    }

    usurf.set(1000.0);  // compute theta for this constant thk

    // actually use the smoother/bed-roughness-parameterizer
    config.set_double("Glen_exponent", 3.0);
    config.set_double("bed_smoother_range", 50.0e3);
    BedSmoother smoother(grid, config, 1);
    smoother.preprocess_bed(topg);
    int Nx,Ny;
    smoother.get_smoothing_domain(Nx,Ny);
    PetscPrintf(grid.com,"  smoothing domain:  Nx = %d, Ny = %d\n",Nx,Ny);
    smoother.get_theta(usurf, &theta);

    const IceModelVec2S &topg_smoothed = smoother.get_smoothed_bed();
    if (show) {
      const int  window = 400;
      topg.view(window);
      topg_smoothed.view(window);
      theta.view(window);
      printf("[showing topg, topg_smoothed, theta in X windows for 10 seconds ...]\n");
      ierr = PetscSleep(10);
      PISM_PETSC_CHK(ierr, "PetscSleep");
    }

    double topg_min, topg_max, topgs_min, topgs_max, theta_min, theta_max;
    topg.min(topg_min);
    topg.max(topg_max);
    topg_smoothed.min(topgs_min);
    topg_smoothed.max(topgs_max);
    theta.min(theta_min);
    theta.max(theta_max);
    PetscPrintf(grid.com,
                "  original bed    :  min elev = %12.6f m,  max elev = %12.6f m\n",
                topg_min, topg_max);
    PetscPrintf(grid.com,
                "  smoothed bed    :  min elev = %12.6f m,  max elev = %12.6f m\n",
                topgs_min, topgs_max);
    PetscPrintf(grid.com,
                "  Schoof's theta  :  min      = %12.9f,    max      = %12.9f\n",
                theta_min, theta_max);

    bool dump = false;
    OptionsIsSet("-dump", dump);
    if (dump) {
      topg.dump("bedrough_test_topg.nc");
      topg_smoothed.dump("bedrough_test_topg_smoothed.nc");
      theta.dump("bedrough_test_theta.nc");
    }

  }
  catch (...) {
    handle_fatal_errors(com);
  }
  return 0;
}
