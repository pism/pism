// Copyright (C) 2009 Ed Bueler and Constantine Khroulev
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
  "Driver for testing PISMClimateCoupler without IceModel.\n";

#include <ctime>
#include <petscda.h>
#include "base/grid.hh"
#include "base/LocalInterpCtx.hh"
#include "base/nc_util.hh"
#include "coupler/pccoupler.hh"


static PetscErrorCode setupIceGridFromFile(const char *filename, IceGrid &grid) {
  PetscErrorCode ierr;

  NCTool nc(&grid);
  ierr = nc.get_grid(filename); CHKERRQ(ierr);
  ierr = grid.createDA(); CHKERRQ(ierr);  
  return 0;
}


static PetscErrorCode readAtmosInfoFromFile(char *filename, IceGrid *grid,
					    LocalInterpCtx* &lic, 
					    IceInfoNeededByAtmosphereCoupler &info) {
  PetscErrorCode ierr;

  info.lat = new IceModelVec2;
  info.lon = new IceModelVec2;
  info.mask = new IceModelVec2;
  info.surfelev = new IceModelVec2;

  ierr = info.lat->create(*grid, "lat", true); CHKERRQ(ierr);
  ierr = info.lat->set_attrs("mapping", "latitude", "degrees_north", "latitude"); CHKERRQ(ierr);

  ierr = info.lon->create(*grid, "lon", true); CHKERRQ(ierr);
  ierr = info.lon->set_attrs("mapping", "longitude", "degrees_east", "longitude"); CHKERRQ(ierr);

  ierr = info.mask->create(*grid, "mask", true); CHKERRQ(ierr);
  ierr = info.mask->set_attrs(NULL, "grounded_dragging_floating integer mask",
			      NULL, NULL); CHKERRQ(ierr);

  ierr = info.surfelev->create(*grid, "usurf", true); CHKERRQ(ierr);
  ierr = info.surfelev->set_attrs(NULL, "ice upper surface elevation",
		                  "m", "surface_altitude"); CHKERRQ(ierr);

  ierr = info.lat->regrid(filename, *lic, true); CHKERRQ(ierr);
  ierr = info.lon->regrid(filename, *lic, true); CHKERRQ(ierr);
  ierr = info.mask->regrid(filename, *lic, true); CHKERRQ(ierr);
  ierr = info.surfelev->regrid(filename, *lic, true); CHKERRQ(ierr);
  return 0;
}


static PetscErrorCode readOceanInfoFromFile(char *filename, MPI_Comm, IceGrid *grid,
                                     LocalInterpCtx* &lic, 
                                     IceInfoNeededByOceanCoupler &info) {
  PetscErrorCode ierr;

  info.lat = new IceModelVec2;
  info.lon = new IceModelVec2;
  info.mask = new IceModelVec2;
  info.thk = new IceModelVec2;

  ierr = info.lat->create(*grid, "lat", true); CHKERRQ(ierr);
  ierr = info.lat->set_attrs("mapping", "latitude", "degrees_north", "latitude"); CHKERRQ(ierr);

  ierr = info.lon->create(*grid, "lon", true); CHKERRQ(ierr);
  ierr = info.lon->set_attrs("mapping", "longitude", "degrees_east", "longitude"); CHKERRQ(ierr);

  ierr = info.mask->create(*grid, "mask", true); CHKERRQ(ierr);
  ierr = info.mask->set_attrs(NULL, "grounded_dragging_floating integer mask",
			      NULL, NULL); CHKERRQ(ierr);

  ierr = info.thk->create(*grid, "thk", true); CHKERRQ(ierr);
  ierr = info.thk->set_attrs(NULL, "land ice thickness",
		             "m", "land_ice_thickness"); CHKERRQ(ierr);

  ierr = info.lat->regrid(filename, *lic, true); CHKERRQ(ierr);
  ierr = info.lon->regrid(filename, *lic, true); CHKERRQ(ierr);
  ierr = info.mask->regrid(filename, *lic, true); CHKERRQ(ierr);
  ierr = info.thk->regrid(filename, *lic, true); CHKERRQ(ierr);
  return 0;
}



static PetscErrorCode doneWithAtmosInfo(IceInfoNeededByAtmosphereCoupler &info) {
  PetscErrorCode ierr;
  ierr = info.lat->destroy(); CHKERRQ(ierr);
  ierr = info.lon->destroy(); CHKERRQ(ierr);
  ierr = info.mask->destroy(); CHKERRQ(ierr);
  ierr = info.surfelev->destroy(); CHKERRQ(ierr);
  delete info.lat;
  delete info.lon;
  delete info.mask;
  delete info.surfelev;
  return 0;
}


static PetscErrorCode doneWithOceanInfo(IceInfoNeededByOceanCoupler &info) {
  PetscErrorCode ierr;
  ierr = info.lat->destroy(); CHKERRQ(ierr);
  ierr = info.lon->destroy(); CHKERRQ(ierr);
  ierr = info.mask->destroy(); CHKERRQ(ierr);
  ierr = info.thk->destroy(); CHKERRQ(ierr);
  delete info.lat;
  delete info.lon;
  delete info.mask;
  delete info.thk;
  return 0;
}


static PetscErrorCode writePCCStateAtTimes(
                 PISMClimateCoupler *pcc,
                 const char *filename, const MPI_Comm com, IceGrid* grid,
                 int argc, char *argv[],
                 PetscReal ys, PetscReal ye, PetscReal dt_years,
                 void* iceInfoNeeded,
		 PolarStereoParams &psparams) {

  PetscErrorCode ierr;
  NCTool nc(grid);

  // put calling command in history string
  char cmdstr[TEMPORARY_STRING_LENGTH], timestr[TEMPORARY_STRING_LENGTH];
  strcpy(cmdstr, "");
  strncat(cmdstr, argv[0], sizeof(cmdstr)); // Does not null terminate on overflow
  cmdstr[sizeof(cmdstr) - 1] = '\0';
  for (PetscInt i=1; i < argc; i++) {
    size_t remaining_bytes = sizeof(cmdstr) - strlen(cmdstr) - 1;
    // strncat promises to null terminate, so we must only make sure that the
    // end of the buffer is not overwritten.
    strncat(cmdstr, " ", remaining_bytes--);
    strncat(cmdstr, argv[i], remaining_bytes);
  }
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
  ierr = nc.write_global_attrs(false, "CF-1.4"); CHKERRQ(ierr);
  ierr = nc.write_polar_stereographic(psparams); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  PetscInt NN;  // get number of times at which PCC state is written
  if (dt_years < 0.0001) {
    ierr = PetscPrintf(com,
      "PCCTEST WARNING: dt_years less than 10^-4 year so just writing state for year %f\n",
      ys); CHKERRQ(ierr);
    NN = 1;
  } else {
    NN = 1 + (int) ceil((ye - ys) / dt_years);
  }
  if (NN > 1000)  SETERRQ(2,"PCCTEST ERROR: refuse to write more than 1000 times!");
  if (NN > 50) {
    ierr = PetscPrintf(com, "\n\nPCCTEST WARNING: writing more than 50 times to '%s'!!\n\n\n",
                       filename); CHKERRQ(ierr);
  }

  PISMPDDCoupler* pdd_pcc = dynamic_cast<PISMPDDCoupler*>(pcc);
  PetscScalar use_dt_years = dt_years;
  if ((pdd_pcc != NULL) && (dt_years > 1.0)) {
    ierr = PetscPrintf(com,
      "PCCTEST ATTENTION: PISMPDDCoupler will be asked for results from one year periods at\n"
      "  the start of each desired time subinterval; full subinterval evaluation is too slow ...\n");
      CHKERRQ(ierr);
    use_dt_years = 1.0;
  }
  
  // write the states
  for (PetscInt k=0; k < NN; k++) {
    const PetscReal pccyear = ys + k * dt_years; // use original dt_years to get correct subinterval starts
    ierr = nc.open_for_writing(filename, false); CHKERRQ(ierr);
    ierr = nc.append_time(pccyear); CHKERRQ(ierr);
    snprintf(timestr, sizeof(timestr), 
             "  coupler updated for [%11.3f a,%11.3f a] ...\n", pccyear, pccyear + use_dt_years);
    ierr = nc.write_history(timestr); CHKERRQ(ierr); // append the history
    ierr = nc.close(); CHKERRQ(ierr);

    ierr = pcc->updateClimateFields(pccyear, use_dt_years, iceInfoNeeded); CHKERRQ(ierr);
    ierr = pcc->writeCouplingFieldsToFile(filename); CHKERRQ(ierr);
    ierr = PetscPrintf(com, "  coupler updated for [%11.3f a,%11.3f a]; result written to %s ...\n",
             pccyear, pccyear + use_dt_years, filename); CHKERRQ(ierr);
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
    char inname[PETSC_MAX_PATH_LEN], outname[PETSC_MAX_PATH_LEN];
    IceGrid grid(com, rank, size);
    PolarStereoParams psparams = {0, 90, -71};
    
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = PetscPrintf(com, 
		       "PCCTEST %s (test of PISMClimateCoupler offline from IceModel)\n",
		       PISM_Revision); CHKERRQ(ierr);
    
    PetscTruth i_set;
    ierr = PetscOptionsGetString(PETSC_NULL, "-i", inname, 
                               PETSC_MAX_PATH_LEN, &i_set); CHKERRQ(ierr);
    if (!i_set) { SETERRQ(1,"PCCTEST ERROR: no -i file to initialize from\n"); }

    ierr = PetscPrintf(com, 
             "  initializing grid from NetCDF file '%s'...\n", inname); CHKERRQ(ierr);
    ierr = setupIceGridFromFile(inname,grid); CHKERRQ(ierr);

    // Process -ys, -ye, -dt. This should happen *before*
    // PCC->initFromOptions() is called.
    PetscReal ys = 0.0, ye = 0.0, dt_years = 0.0;
    PetscTruth ysSet, yeSet, dtSet;
    ierr = PetscOptionsGetReal(PETSC_NULL, "-ys", &ys, &ysSet); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL, "-ye", &ye, &yeSet); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL, "-dt", &dt_years, &dtSet); CHKERRQ(ierr);

    if (!ysSet || !yeSet || !dtSet) {
      ierr = PetscPrintf(com, "PCCTEST ERROR: All three of -ys, -ye, -dt are required.\n");
      CHKERRQ(ierr);
      PetscEnd();
    }
    grid.year = ys;		// this value is used in PCC->initFromOptions()

    // set PCC from options
    PetscTruth caSet, pddSet, coSet;
    ierr = check_option("-ca",  caSet); CHKERRQ(ierr);
    ierr = check_option("-pdd", pddSet); CHKERRQ(ierr);
    ierr = check_option("-co",  coSet); CHKERRQ(ierr);
    int  choiceSum = (int) caSet + (int) pddSet + (int) coSet;
    if (choiceSum == 0) {
      ierr = PetscPrintf(com,"PCCTEST ERROR: called with no derived class\n");
         CHKERRQ(ierr);    PetscEnd();
    } else if (choiceSum > 1) {
      ierr = PetscPrintf(com,"PCCTEST ERROR: called with more than one derived class\n");
         CHKERRQ(ierr);    PetscEnd();
    }

    PISMConstAtmosCoupler pcac;
    PISMPDDCoupler        ppdd;
    PISMConstOceanCoupler pcoc;
    PISMClimateCoupler*   PCC;
    if (caSet == PETSC_TRUE) { 
      PCC = (PISMClimateCoupler*) &pcac;
    } else if (pddSet == PETSC_TRUE) { 
      PCC = (PISMClimateCoupler*) &ppdd;
    } else if (coSet == PETSC_TRUE) { 
      PCC = (PISMClimateCoupler*) &pcoc;
    } else {
      PCC = PETSC_NULL;
      ierr = PetscPrintf(com,"PCCTEST ERROR: how did I get here?  111\n"); CHKERRQ(ierr);
      PetscEnd();
    }
    
    ierr = PCC->initFromOptions(&grid); CHKERRQ(ierr);

    LocalInterpCtx* lic;
    ierr = PCC->findPISMInputFile(inname, lic); CHKERRQ(ierr); // allocates lic

    // Get the polar stereographic projection parameters.
    NCTool nc(&grid);
    ierr = nc.open_for_reading(inname); CHKERRQ(ierr);
    ierr = nc.read_polar_stereographic(psparams); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);

    IceInfoNeededByAtmosphereCoupler info_atmos;
    IceInfoNeededByOceanCoupler      info_ocean;
    void*                            info;
    if ((caSet == PETSC_TRUE) || (pddSet == PETSC_TRUE)) { 
      ierr = PetscPrintf(com, 
             "  reading fields 'usurf', 'mask', 'lat', 'lon' from NetCDF file '%s' to fill fields\n"
             "    in IceInfoNeededByAtmosphereCoupler ...\n",
             inname); CHKERRQ(ierr);
      ierr = readAtmosInfoFromFile(inname,&grid,lic,info_atmos); CHKERRQ(ierr);
      info = (void*) &info_atmos;
    } else if (coSet == PETSC_TRUE) {
      ierr = PetscPrintf(com, 
             "  reading fields 'thk', 'mask', 'lat', 'lon' from NetCDF file '%s' to fill fields\n"
             "    in IceInfoNeededByOceanCoupler ...\n",
             inname); CHKERRQ(ierr);
      ierr = readOceanInfoFromFile(inname,com,&grid,lic,info_ocean); CHKERRQ(ierr);
      info = (void*) &info_ocean;
    } else {
      info = PETSC_NULL;
      ierr = PetscPrintf(com,"PCCTEST ERROR: how did I get here?  222\n"); CHKERRQ(ierr);
      PetscEnd();
    }
    
    PetscTruth oSet;
    ierr = PetscOptionsGetString(PETSC_NULL, "-o", outname, 
                               PETSC_MAX_PATH_LEN, &oSet); CHKERRQ(ierr);
    if (oSet != PETSC_TRUE) { SETERRQ(2,"PCCTEST ERROR: no -o file to write to\n"); }

    ierr = PetscPrintf(com, "  writing PISMClimateCoupler states to NetCDF file '%s'...\n",
                       outname); CHKERRQ(ierr);
    ierr = writePCCStateAtTimes(PCC,outname,com,&grid, argc,argv, ys,ye,dt_years,
                                info, psparams); CHKERRQ(ierr);

    if ((caSet == PETSC_TRUE) || (pddSet == PETSC_TRUE)) { 
      ierr = doneWithAtmosInfo(info_atmos); CHKERRQ(ierr);
    } else if (coSet == PETSC_TRUE) {
      ierr = doneWithOceanInfo(info_ocean); CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(com,"PCCTEST ERROR: how did I get here?  333\n"); CHKERRQ(ierr);
      PetscEnd();
    }

    ierr = PetscPrintf(com, "... done\n"); CHKERRQ(ierr);
    
    delete lic;
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
