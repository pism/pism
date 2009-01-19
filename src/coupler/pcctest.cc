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
#include "pPDDcoupler.hh"


/* example run: here pcctest extracts 'mask' and 'usurf' = surface elevation from sheet.nc;
   extracts 'acab' and 'artm' from sheet.nc, because they are stored there by EISMINT
   II choices, and initializes PISMConstAtmosCoupler;
   (FIXME: will allow PISMPDDCoupler instead); "computes", i.e. no actual computation,
   surface mass balance and surface temp at ys:dt:ye and saves these in pccmovie.nc

$ mpiexec -n 2 pisms -eisII A -y 10000 -o state.nc

$ pcctest -if state.nc -ys 0.0 -ye 300.0 -dt 100.0 -of pccmovie.nc

*/

PetscErrorCode setupIceGridFromFile(const char *filename, const MPI_Comm com, IceGrid &grid) {
  PetscErrorCode ierr;
  bool file_exists = false;
  grid_info gi;
  double *zlevs, *zblevs;

  NCTool nc(&grid);

  // read grid params
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
  zlevs = new double[grid.Mz];
  zblevs = new double[grid.Mbz];
  ierr = nc.get_vertical_dims(grid.Mz, grid.Mbz, zlevs, zblevs); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

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
  delete [] zlevs;  delete [] zblevs;
  
  // finally, set DA 
  ierr = grid.createDA(); CHKERRQ(ierr);
  // sets grid.Lx, grid.Ly:
  ierr = grid.rescale_using_zlevels(-gi.x_min, -gi.y_min); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode writePCCStateAtTimes(
                 PISMClimateCoupler *pcc,
                 const char *filename, const MPI_Comm com, IceGrid* grid,
                 int argc, char *argv[],
                 const PetscReal ys, const PetscReal ye, const PetscReal dt_years) {

  PetscErrorCode ierr;
  NCTool nc(grid);

  // put calling command in history string
  char cmdstr[TEMPORARY_STRING_LENGTH], timestr[TEMPORARY_STRING_LENGTH];
  strcpy(cmdstr, "");
  strncat(cmdstr, argv[0], sizeof(cmdstr)); // Does not null terminate on overflow
  cmdstr[sizeof(cmdstr) - 1] = '\0';
  for (PetscInt i=1; i < argc; i++) {
    PetscInt remaining_bytes = sizeof(cmdstr) - strlen(cmdstr) - 1;
    // strncat promises to null terminate, so we must only make sure that the
    // end of the buffer is not overwritten.
    strncat(cmdstr, " ", remaining_bytes--);
    strncat(cmdstr, argv[i], remaining_bytes);
  }
//  strcat(cmdstr,"\n");
  cmdstr[sizeof(cmdstr) - 1] = '\0';
  
  // compare IceModel::stampHistory() for this way of getting date etc in file
  time_t now;
  tm tm_now;
  now = time(NULL);
  localtime_r(&now, &tm_now);
  char date_str[50], username[50], hostname[100], wwstr[TEMPORARY_STRING_LENGTH];
  strftime(date_str, sizeof(date_str), "%F %T %Z", &tm_now);
  ierr = PetscGetUserName(username, sizeof(username)); CHKERRQ(ierr);
  ierr = PetscGetHostName(hostname, sizeof(hostname)); CHKERRQ(ierr);
  int length = snprintf(wwstr, sizeof(wwstr), "%s@%s %s:  %s\n",
                        username, hostname, date_str, cmdstr);  
  if (length < 0) {
    SETERRQ(3, "PCCTEST ERROR: snprintf() is not C99 compliant?");
  }
  if (length > (int)sizeof(wwstr)) {
    ierr = PetscPrintf(com,
       "PCCTEST WARNING: command line truncated to %d chars in history.\n",
       length + 1 - sizeof(wwstr)); CHKERRQ(ierr);
    wwstr[sizeof(wwstr) - 2] = '\n';
    wwstr[sizeof(wwstr) - 1] = '\0';
  }

  ierr = nc.open_for_writing(filename, true); CHKERRQ(ierr);
  ierr = nc.write_history(wwstr); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  // get NN = number of times at which PCC state is written
  PetscInt NN;
  if (dt_years < 0.0001) {
    ierr = PetscPrintf(com,
      "PCCTEST WARNING: dt_years less than 10^-4 year so just writing state for year %f\n",
      ys); CHKERRQ(ierr);
    NN = 1;
  } else {
    NN = 1 + (int) ceil((ye - ys) / dt_years);
  }
  if (NN > 1000)  SETERRQ(2,"PCCTEST ERROR: refuse to write more than 1000 times");
  if (NN > 50) {
    ierr = PetscPrintf(com, "\n\nPCCTEST WARNING: writing more than 50 times to '%s'!!\n\n\n",
                       filename); CHKERRQ(ierr);
  }

  // write the states
  for (PetscInt k=0; k < NN; k++) {
    const PetscReal pccyear = ys + k * dt_years;
    ierr = nc.open_for_writing(filename, false); CHKERRQ(ierr);
    ierr = nc.append_time(pccyear * secpera); CHKERRQ(ierr);
    snprintf(timestr, sizeof(timestr), "  pcc state at year %11.3f ...\n", pccyear);
    ierr = nc.write_history(timestr); CHKERRQ(ierr); // append the history
    ierr = nc.close(); CHKERRQ(ierr);

    // FIXME:  need to fill struct* with ice info needed for non-trivial update; here NULL
    ierr = pcc->updateClimateFields(pccyear, dt_years, NULL); CHKERRQ(ierr);

    ierr = pcc->writeCouplingFieldsToFile(filename); CHKERRQ(ierr);
  }

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
    PISMPDDCoupler        ppddc;
    PISMClimateCoupler*   PCC;
    
    ierr = PetscPrintf(com, "PCCTEST (test of climate couplers offline from PISM)\n"); CHKERRQ(ierr);
    
    PetscTruth optSet;
    ierr = PetscOptionsGetString(PETSC_NULL, "-if", inname, 
                               PETSC_MAX_PATH_LEN, &optSet); CHKERRQ(ierr);
    if (optSet != PETSC_TRUE) { SETERRQ(1,"PCCTEST ERROR: no -if to initialize from\n"); }

    ierr = PetscPrintf(com, "  initializing grid from NetCDF file '%s'...\n",
                       inname); CHKERRQ(ierr);
    ierr = setupIceGridFromFile(inname,com,grid); CHKERRQ(ierr);
    
    // FIXME:  set PCC from options
    PCC = (PISMClimateCoupler*) &pcac;
    
    ierr = PCC->initFromOptions(&grid); CHKERRQ(ierr);

    ierr = PetscOptionsGetString(PETSC_NULL, "-of", outname, 
                               PETSC_MAX_PATH_LEN, &optSet); CHKERRQ(ierr);
    if (optSet != PETSC_TRUE) { SETERRQ(2,"PCCTEST ERROR: no -of to write to\n"); }

    PetscReal ys = 0.0, ye = 0.0, dt_years = 0.0;
    ierr = PetscOptionsGetReal(PETSC_NULL, "-ys", &ys, PETSC_NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL, "-ye", &ye, PETSC_NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL, "-dt", &dt_years, PETSC_NULL); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "  writing PISMClimateCoupler states to NetCDF file '%s'...\n",
                       outname); CHKERRQ(ierr);
    ierr = writePCCStateAtTimes(PCC,outname,com,&grid,
                                argc,argv,ys,ye,dt_years); CHKERRQ(ierr);

    ierr = PetscPrintf(com, "... done\n"); CHKERRQ(ierr);
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
