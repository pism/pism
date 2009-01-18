// Copyright (C) 2009 Ed Bueler
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
  "Driver which allows testing PISMClimateCoupler without IceModel.\n";

#include <petscda.h>
#include "../base/grid.hh"
#include "../base/LocalInterpCtx.hh"
#include "../base/nc_util.hh"
#include "pccoupler.hh"
// FIXME: #include "pPDDcoupler.hh"


/* example run: pcctest extracts 'mask' and 'usurf' = surface elevation from sheet.nc;
   extracts 'acab' and 'artm' from sheet.nc, because they are stored there by EISMINT
   II choices, and initializes PISMConstAtmosCoupler;
   (FIXME: will allow PISMPDDCoupler instead); "computes", i.e. no actual computation,
   surface mass balance and surface temp at ys:dt:ye and saves these in pccmovie.nc

$ mpiexec -n 2 pisms -eisII A -y 10000 -o sheet.nc

$ pcctest -if sheet.nc -ys 0.0 -ye 1000.0 -dt 100.0 -of pccmovie.nc

*/

PetscErrorCode setupIceGridFromFile(const char *filename, const MPI_Comm com, IceGrid &grid) {
  PetscErrorCode ierr;
  bool file_exists = false;
  grid_info gi;
  double *zlevs, *zblevs;

  NCTool nc(&grid);

  ierr = nc.open_for_reading(filename, file_exists); CHKERRQ(ierr);
  if (!file_exists) {
    ierr = PetscPrintf(com, "PCCTEST ERROR: can't open file '%s'\n", filename); CHKERRQ(ierr);
    PetscEnd();
  }

  ierr = nc.get_grid_info(gi);
  grid.year = gi.time / secpera;
  grid.Mx = gi.x_len;
  grid.My = gi.y_len;
  grid.Mz = gi.z_len;
  grid.Mbz = gi.zb_len;
  // grid.Lx, grid.Ly set from gi.x_max, gi.y_max below in call to grid.rescale_using_zlevels()

  zlevs = new double[grid.Mz];
  zblevs = new double[grid.Mbz];
  ierr = nc.get_vertical_dims(grid.Mz, grid.Mbz, zlevs, zblevs); CHKERRQ(ierr);
  // re-allocate and fill grid.zlevels & zblevels with read values
  delete [] grid.zlevels;  delete [] grid.zblevels;
  grid.zlevels = new PetscScalar[grid.Mz];
  grid.zblevels = new PetscScalar[grid.Mbz];
  for (PetscInt k = 0; k < grid.Mz; k++) {
    grid.zlevels[k] = (PetscScalar) zlevs[k];
  }
  for (PetscInt k = 0; k < grid.Mbz; k++) {
    grid.zblevels[k] = (PetscScalar) zblevs[k];
  }
  delete [] zlevs;  delete [] zblevs;  // done with these
  
  // finally, set DA 
  ierr = grid.createDA(); CHKERRQ(ierr);
  ierr = grid.rescale_using_zlevels(-gi.x_min, -gi.y_min); CHKERRQ(ierr);
  return 0;
}


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
    PetscErrorCode ierr;
    char inname[PETSC_MAX_PATH_LEN], outname[PETSC_MAX_PATH_LEN];

    IceGrid grid(com, rank, size);

    PISMConstAtmosCoupler pcac;
    // PISMPDDCoupler ppddc;  FIXME: test this one too, in same program; user chooses
    
    ierr = PetscPrintf(com, "PCCTEST (test of climate couplers offline from PISM)\n"); CHKERRQ(ierr);
    
    PetscTruth optSet;
    ierr = PetscOptionsGetString(PETSC_NULL, "-if", inname, 
                               PETSC_MAX_PATH_LEN, &optSet); CHKERRQ(ierr);
    if (optSet != PETSC_TRUE) {
      ierr = PetscPrintf(com, "PCCTEST ERROR: no -if to initialize from\n"); CHKERRQ(ierr);
      PetscEnd();
    }

    ierr = PetscPrintf(com, "  initializing grid from NetCDF file '%s'...\n",
                       inname); CHKERRQ(ierr);
    ierr = setupIceGridFromFile(inname,com,grid); CHKERRQ(ierr);
    
    ierr = pcac.initFromOptions(&grid); CHKERRQ(ierr);

    ierr = PetscOptionsGetString(PETSC_NULL, "-of", outname, 
                               PETSC_MAX_PATH_LEN, &optSet); CHKERRQ(ierr);
    if (optSet != PETSC_TRUE) {
      ierr = PetscPrintf(com, "PCCTEST ERROR: no -of to write to\n"); CHKERRQ(ierr);
      PetscEnd();
    }

    NCTool nc(&grid);
    PetscInt NN = 3; // FIXME: read options -ys, -ye, -dt to get when to write
    for (PetscInt k=0; k<NN; k++) {
      ierr = nc.open_for_writing(outname, false); CHKERRQ(ierr); // replace == false
      ierr = nc.append_time(grid.year * secpera); CHKERRQ(ierr);
      //ierr = nc.write_history(tmp); CHKERRQ(ierr); // append the history
      ierr = nc.close(); CHKERRQ(ierr);

      // FIXME: call updateSurfMassFluxAndProvide() and updateSurfTempAndProvide()
      //        here so that nontrivial PISMClimateCoupler can change fields 

      ierr = pcac.writeCouplingFieldsToFile(outname);
    }
    
    ierr = PetscPrintf(com, "... done\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
