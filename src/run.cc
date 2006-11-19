// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

static char help[] =
  "Ice sheet driver for ice sheet simulations initialized from data.\n"
  "Allows use of MacAyeal-Morland equations to compute velocities \n"
  "for ice shelves and dragging ice shelves (i.e. ice streams).\n"
  "Includes Goldsby-Kohlstedt law, age field calculation, age-dependent \n"
  "grain size, spatially-varying geothermal heat flux, interface to NetCDF files,\n"
  "arbitrary grids chosen at the command line, regridding, and many other features.\n"
  "(Jed Kallen-Brown and Ed Bueler, authors, and Craig Lingle, concepts.).\n";

#include <petscbag.h>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"

int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    IceGrid g(com, rank, size);
    PetscInt   flowlawNumber = 0;  // use Paterson-Budd by default
    IceType*   ice;

    ierr = PetscPrintf(com, "PISMR (primary run mode)\n"); CHKERRQ(ierr);
    
    ierr = getFlowLawFromUser(com, ice, flowlawNumber); CHKERRQ(ierr);
    IceModel m(g, *ice);
    ierr = m.setFromOptions(); CHKERRQ(ierr);
    ierr = m.initFromOptions(); CHKERRQ(ierr);
    ierr = m.setSoundingFromOptions(); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "running ...\n"); CHKERRQ(ierr);
    ierr = m.run(); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "... done with run\n"); CHKERRQ(ierr);

    // We provide a default base name.  If the user does not specify an
    // output file name, then we will write one with the following base in
    // the default format.
    // At this time, the default (and only supported output format) is Petsc
    // binary format.  Thus we append `.pb' to the default base name and write
    // a Petsc binary output file.
    ierr = m.writeFiles("unnamed"); CHKERRQ(ierr);

    ierr = PetscPrintf(com, " ... done.\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
