// Copyright (C) 2004-2007 Ed Bueler and Nathan Shemonski
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
  "Ice sheet driver for EISMINT Greenland ice sheet simulations.\n";

#include <cstring>
#include <petscbag.h>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"
#include "iceGRNModel.hh"

int main(int argc, char *argv[]){
  PetscErrorCode ierr;

  MPI_Comm com;
  PetscMPIInt rank, size;
  
  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);
  
  // forces calling of deconstrucors before PetscFinalize() 
  {
    IceGrid g(com, rank, size);
    IceType *ice;
    PetscInt flowlawNumber = 0;
 
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = verbPrintf(1, com, "PGRN (EISMINT Greenland mode)\n"); CHKERRQ(ierr);

    ierr = getFlowLawFromUser(com, ice, flowlawNumber); CHKERRQ(ierr);
 
    IceGRNModel mGRN(g, *ice);
    
    ierr = mGRN.setFromOptions(); CHKERRQ(ierr);
    ierr = mGRN.initFromOptions(); CHKERRQ(ierr);
 
    ierr = mGRN.run(); CHKERRQ(ierr);
    ierr = verbPrintf(2, com, "done with run ... \n"); CHKERRQ(ierr);

    // FIXME: fix data for the netCDF file
    ierr = mGRN.removeBedDiff(); CHKERRQ(ierr);

    ierr = mGRN.writeFiles("grn_exper"); CHKERRQ(ierr);

    ierr = verbPrintf(2, com, " ... done.\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
