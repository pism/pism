// Copyright (C) 2007--2010 Ed Bueler and Constantine Khroulev
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

// get this message with 'tryLCbd -help':
static char help[] =
  "Simple testing program for Schoof (2003)-based bed smoothing and bed roughness\n"
  "  parameterization schemes.\n\n";

#include <cmath>
#include <cstdio>
#include <petscvec.h>
#include <petscda.h>
#include "../base/pism_const.hh"
#include "../base/grid.hh"
#include "../base/NCVariable.hh"
#include "../base/PISMBedSmoother.hh"

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
    ierr = setVerbosityLevel(2); CHKERRQ(ierr);

    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    IceGrid grid(com, rank, size, config);
    grid.Mx = 81;
    grid.My = 81;
    grid.Lx = 1200e3;
    grid.Ly = grid.Lx;
    grid.compute_nprocs();
    grid.compute_ownership_ranges();
    ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
    double *x, *y;
    ierr = grid.compute_horizontal_coordinates(x, y); CHKERRQ(ierr);
    //ierr = grid.printInfo(1); CHKERRQ(ierr);
    ierr = grid.createDA(); CHKERRQ(ierr);

    IceModelVec2S topg, thk, theta;
    ierr = topg.create(grid, "topg", false); CHKERRQ(ierr);
    ierr = topg.set_attrs(
      "trybedrough_tool", "original topography",
      "m", "bedrock_altitude"); CHKERRQ(ierr);
    ierr = thk.create(grid, "thk", false); CHKERRQ(ierr);
    ierr = thk.set_attrs(
      "trybedrough_tool", "ice thickness",
      "m", "land_ice_thickness"); CHKERRQ(ierr);
    ierr = theta.create(grid, "theta", false); CHKERRQ(ierr);
    ierr = theta.set_attrs(
      "trybedrough_tool",
      "coefficient theta in Schoof (2003) bed roughness parameterization",
      "", ""); CHKERRQ(ierr);

    // put in bed elevations, a la this Matlab:
    //    topg0 = 400 * sin(2 * pi * xx / 600e3) + ...
    //            100 * sin(2 * pi * (xx + 1.5 * yy) / 40e3);
    ierr = topg.begin_access(); CHKERRQ(ierr);
    for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
        topg(i,j) = 400.0 * sin(2.0 * pi * x[i] / 600.0e3) +
                    100.0 * sin(2.0 * pi * (x[i] + 1.5 * y[j]) / 40.0e3);
      }
    }
    ierr = topg.end_access(); CHKERRQ(ierr);

    ierr = thk.set(1000.0); CHKERRQ(ierr);  // compute theta for this constant thk

    // actually use the smoother/bed-roughness-parameterizer
    PISMBedSmoother smoother(grid, config);
    ierr = smoother.preprocess_bed(topg, 50.0e3, 3.0); CHKERRQ(ierr);
    ierr = smoother.get_theta(thk, 3.0, &theta); CHKERRQ(ierr);

    const PetscInt  window = 400;
    ierr = topg.view(window);  CHKERRQ(ierr);
    ierr = smoother.topgsmooth.view(window);  CHKERRQ(ierr);
    ierr = theta.view(window);  CHKERRQ(ierr);

    ierr = PetscSleep(15); CHKERRQ(ierr);

// FIXME ideas:  write theta in ascii and compare to Matlab/Octave out of same thing

    delete [] x; delete [] y;
  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
