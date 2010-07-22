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
    NCConfigVariable config, overrides;
    ierr = init_config(com, rank, config, overrides); CHKERRQ(ierr);

    IceGrid grid(com, rank, size, config);
    grid.Mx = 81;
    grid.My = 81;
    grid.Lx = 1200e3;
    grid.Ly = grid.Lx;
    ierr = grid.createDA(); CHKERRQ(ierr);
    ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);

    double *x, *y;
    ierr = grid.compute_horizontal_coordinates(x, y); CHKERRQ(ierr);

    PISMBedSmoother smoother(grid, config);

    IceModelVec2S topg, theta;

    const PetscInt  window = 500;
    ierr = topg.view(window);  CHKERRQ(ierr);
    ierr = smoother.topgsmooth.view(window);  CHKERRQ(ierr);
    ierr = theta.view(window);  CHKERRQ(ierr);


    // FIXME do nothing


  }
  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
