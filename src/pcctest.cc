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
  "Driver for testing PISMClimateCoupler and its derived classes without IceModel.\n";

#include <ctime>
#include <string>
#include <vector>
#include <petscda.h>
#include "base/pism_const.hh"
#include "base/grid.hh"
#include "base/LocalInterpCtx.hh"
#include "base/nc_util.hh"
#include "coupler/pccoupler.hh"
#include "coupler/pGreenlandAtmosCoupler.hh"
#include "base/NCVariable.hh"


static PetscErrorCode setupIceGridFromFile(const char *filename, IceGrid &grid) {
  PetscErrorCode ierr;

  NCTool nc(&grid);
  ierr = nc.get_grid(filename); CHKERRQ(ierr);
  ierr = grid.createDA(); CHKERRQ(ierr);  
  return 0;
}

static PetscErrorCode createVecs(IceGrid &grid, PISMVars &variables) {
  
  PetscErrorCode ierr;
  IceModelVec2 *lat, *lon, *mask, *thk, *surfelev, *topg;

  lat = new IceModelVec2;
  lon = new IceModelVec2;
  mask = new IceModelVec2;
  thk = new IceModelVec2;
  surfelev = new IceModelVec2;
  topg = new IceModelVec2;

  ierr = lat->create(grid, "lat", true); CHKERRQ(ierr);
  ierr = lat->set_attrs("mapping", "latitude", "degrees_north", "latitude"); CHKERRQ(ierr);
  ierr = variables.add(*lat); CHKERRQ(ierr);

  ierr = lon->create(grid, "lon", true); CHKERRQ(ierr);
  ierr = lon->set_attrs("mapping", "longitude", "degrees_east", "longitude"); CHKERRQ(ierr);
  ierr = variables.add(*lon); CHKERRQ(ierr);

  ierr = mask->create(grid, "mask", true); CHKERRQ(ierr);
  ierr = mask->set_attrs("", "grounded_dragging_floating integer mask",
			      "", ""); CHKERRQ(ierr);
  ierr = variables.add(*mask); CHKERRQ(ierr);

  ierr = thk->create(grid, "thk", true); CHKERRQ(ierr);
  ierr = thk->set_attrs("", "land ice thickness",
		             "m", "land_ice_thickness"); CHKERRQ(ierr);
  ierr = variables.add(*thk); CHKERRQ(ierr);

  ierr = surfelev->create(grid, "usurf", true); CHKERRQ(ierr);
  ierr = surfelev->set_attrs("", "ice upper surface elevation",
		                  "m", "surface_altitude"); CHKERRQ(ierr);
  ierr = variables.add(*surfelev); CHKERRQ(ierr);

  ierr = topg->create(grid, "topg", true); CHKERRQ(ierr);
  ierr = topg->set_attrs("", "bedrock surface elevation",
			"m", "bedrock_altitude"); CHKERRQ(ierr);
  ierr = variables.add(*topg); CHKERRQ(ierr);

  return 0;
}

static PetscErrorCode readIceInfoFromFile(char *filename, LocalInterpCtx* &lic, int start,
                                          PISMVars &variables) {
  PetscErrorCode ierr;

  if (lic != NULL) {
    ierr = variables.get("lat")->regrid(filename, *lic, true); CHKERRQ(ierr);
    ierr = variables.get("lon")->regrid(filename, *lic, true); CHKERRQ(ierr);
    ierr = variables.get("mask")->regrid(filename, *lic, true); CHKERRQ(ierr);
    ierr = variables.get("thk")->regrid(filename, *lic, true); CHKERRQ(ierr);
    ierr = variables.get("usurf")->regrid(filename, *lic, true); CHKERRQ(ierr);
    ierr = variables.get("topg")->regrid(filename, *lic, true); CHKERRQ(ierr);
  } else {
    ierr = variables.get("lat")->read(filename, start); CHKERRQ(ierr);
    ierr = variables.get("lon")->read(filename, start); CHKERRQ(ierr);
    ierr = variables.get("mask")->read(filename, start); CHKERRQ(ierr);
    ierr = variables.get("thk")->read(filename, start); CHKERRQ(ierr);
    ierr = variables.get("usurf")->read(filename, start); CHKERRQ(ierr);
    ierr = variables.get("topg")->read(filename, start); CHKERRQ(ierr);
  }

  return 0;
}


static PetscErrorCode doneWithIceInfo(PISMVars &variables) {

  delete variables.get("lat");
  delete variables.get("lon");
  delete variables.get("mask");
  delete variables.get("thk");
  delete variables.get("usurf");
  delete variables.get("topg");
  
  return 0;
}


static PetscErrorCode writePCCStateAtTimes(
                 PISMClimateCoupler *pcc,
                 const char *filename, const MPI_Comm com, IceGrid* grid,
                 int argc, char *argv[],
                 PetscReal ys, PetscReal ye, PetscReal dt_years,
		 NCConfigVariable &mapping) {

  PetscErrorCode ierr;
  NCTool nc(grid);
  NCGlobalAttributes global_attrs;

  global_attrs.init("global_attributes", com, grid->rank);
  global_attrs.set_string("Conventions", "CF-1.4");
  global_attrs.set_string("source", string("pcctest ") + PISM_Revision);

  // Create a string with space-separated command-line arguments:
  string cmdstr;
  for (int j = 0; j < argc; j++)
    cmdstr += string(" ") + argv[j];
  cmdstr += "\n";

  string history = username_prefix() + cmdstr;

  global_attrs.prepend_history(history);

  ierr = nc.open_for_writing(filename, false, true); CHKERRQ(ierr);
  // append == false, check_dims == true
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = mapping.write(filename); CHKERRQ(ierr);
  ierr = global_attrs.write(filename); CHKERRQ(ierr);

  PetscInt NN;  // get number of times at which PCC state is written
  NN = (int) ceil((ye - ys) / dt_years);
  if (NN > 1000)
    SETERRQ(2,"PCCTEST ERROR: refuse to write more than 1000 times!");
  if (NN > 50) {
    ierr = PetscPrintf(com,
        "\nPCCTEST ATTENTION: writing more than 50 times to '%s'!!\n\n",
        filename); CHKERRQ(ierr);
  }

  PISMGreenlandAtmosCoupler* pdd_pcc = dynamic_cast<PISMGreenlandAtmosCoupler*>(pcc);
  PetscScalar use_dt_years = dt_years;
  if ((pdd_pcc != NULL) && (dt_years > 1.0)) {
    ierr = verbPrintf(1,com,
        "PCCTEST ATTENTION: PISMGreenlandAtmosCoupler will be asked for results\n"
        "  from one year periods at the start of each desired time subinterval;\n"
        "  full subinterval evaluation is too slow ...\n\n");
      CHKERRQ(ierr);
    use_dt_years = 1.0;
  }
  
  // write the states
  for (PetscInt k=0; k < NN; k++) {
    // use original dt_years to get correct subinterval starts:
    const PetscReal pccyear = ys + k * dt_years; 
    ierr = nc.open_for_writing(filename, true, false); CHKERRQ(ierr); // append=true,check_dims=false
    ierr = nc.append_time(pccyear); CHKERRQ(ierr);
    
    PetscScalar dt_update_years = PetscMin(use_dt_years, ye - pccyear);

    ierr = pcc->updateClimateFields(pccyear, dt_update_years); CHKERRQ(ierr);

    char timestr[TEMPORARY_STRING_LENGTH];
    snprintf(timestr, sizeof(timestr), 
        "  coupler updated for [%11.3f a,%11.3f a] ...", 
        pccyear, pccyear + dt_update_years);
    ierr = verbPrintf(2,com,"."); CHKERRQ(ierr);
    ierr = verbPrintf(3,com,"\n%s writing result to %s ..",timestr,filename); CHKERRQ(ierr);
    strncat(timestr,"\n",1);
    ierr = nc.write_history(timestr); CHKERRQ(ierr); // append the history
    ierr = nc.close(); CHKERRQ(ierr);

    ierr = pcc->writeCouplingFieldsToFile(pccyear, filename); CHKERRQ(ierr);

  }
  ierr = verbPrintf(2,com,"\n"); CHKERRQ(ierr);

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
    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);

    vector<string> required;
    required.push_back("-i");
    required.push_back("-o");
    required.push_back("-ys");
    required.push_back("-ye");
    required.push_back("-dt");
    ierr = show_usage_check_req_opts(com, "pcctest", required,
      "  pcctest -i IN.nc -o OUT.nc -ys A -ye B -dt C [-ca|-co|-greenland] [OTHER PISM & PETSc OPTIONS]\n\n"
      "where:\n"
      "  -i         input file in NetCDF format\n"
      "  -o         output file in NetCDF format\n"
      "  -ys        start time A (= float) in years\n"
      "  -ye        end time B (= float), B > A, in years\n"
      "  -dt        time step C (= positive float) in years\n"
      "and choose instance of PISMClimateCoupler:\n"
      "  -ca        apply PISMConstAtmosCoupler (default)\n"
      "  -co        apply PISMConstOceanCoupler\n"
      "  -greenland apply PISMGreenlandAtmosCoupler\n"
      ); CHKERRQ(ierr);

    ierr = verbPrintf(2,
      com,"PCCTEST %s (test of PISMClimateCoupler offline from IceModel)\n",
      PISM_Revision); CHKERRQ(ierr);
    
    char inname[PETSC_MAX_PATH_LEN], outname[PETSC_MAX_PATH_LEN];
    IceGrid grid(com, rank, size);
    NCConfigVariable mapping;
    
    ierr = PetscOptionsGetString(PETSC_NULL, "-i", inname, 
                                 PETSC_MAX_PATH_LEN, NULL); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, 
             "  initializing grid from NetCDF file %s...\n", inname); CHKERRQ(ierr);
    ierr = setupIceGridFromFile(inname,grid); CHKERRQ(ierr);

    mapping.init("mapping", com, rank);

    // Process -ys, -ye, -dt *before* PCC->initFromOptions() is called.
    PetscReal ys = 0.0, ye = 0.0, dt_years = 0.0;
    ierr = PetscOptionsGetReal(PETSC_NULL, "-ys", &ys, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL, "-ye", &ye, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetReal(PETSC_NULL, "-dt", &dt_years, NULL); CHKERRQ(ierr);
    // These values are used in PCC->initFromOptions()
    grid.year = ys;		
    grid.start_year = ys;
    grid.end_year = ye;

    // set PCC from options
    PetscTruth caSet, greenlandSet, coSet;
    ierr = check_option("-ca",  caSet); CHKERRQ(ierr);
    ierr = check_option("-co",  coSet); CHKERRQ(ierr);
    ierr = check_option("-greenland", greenlandSet); CHKERRQ(ierr);

    PISMConstAtmosCoupler     pcac;
    PISMConstOceanCoupler     pcoc;
    PISMGreenlandAtmosCoupler pgac;

    PISMClimateCoupler*       PCC;
    if (coSet == PETSC_TRUE) { 
      PCC = (PISMClimateCoupler*) &pcoc;
    } else if (greenlandSet == PETSC_TRUE) { 
      PCC = (PISMClimateCoupler*) &pgac;
    } else  { 
      PCC = (PISMClimateCoupler*) &pcac;
    }

    // allocate IceModelVecs needed by couplers and put them in a dictionary:
    PISMVars variables;
    ierr = createVecs(grid, variables); CHKERRQ(ierr);

    ierr = PCC->initFromOptions(&grid, variables); CHKERRQ(ierr);

    LocalInterpCtx* lic;
    bool regrid;
    int start;
    // allocates lic (if necessary)
    ierr = PCC->findPISMInputFile(inname, lic, regrid, start); CHKERRQ(ierr);

    // Get the projection parameters.
    NCTool nc(&grid);
    ierr = nc.open_for_reading(inname); CHKERRQ(ierr);
    bool mapping_exists;
    ierr = nc.find_variable("mapping", NULL, mapping_exists); CHKERRQ(ierr);
    ierr = nc.close(); CHKERRQ(ierr);
    if (mapping_exists) {
      ierr = mapping.read(inname); CHKERRQ(ierr);
      ierr = mapping.print(); CHKERRQ(ierr);
    }

    ierr = verbPrintf(2,com, 
             "  reading fields lat,lon,mask,thk,topg,usurf from NetCDF file %s\n"
             "    to fill fields in PISMVars ...\n",
             inname); CHKERRQ(ierr);

    ierr = readIceInfoFromFile(inname, lic, start, variables); CHKERRQ(ierr);
    
    ierr = PetscOptionsGetString(PETSC_NULL, "-o", outname, 
                                 PETSC_MAX_PATH_LEN, NULL); CHKERRQ(ierr);

    ierr = verbPrintf(2,
      com, "  writing PISMClimateCoupler states to NetCDF file '%s'...\n",
      outname); CHKERRQ(ierr);
    ierr = writePCCStateAtTimes(PCC,outname,com,&grid, argc,argv, ys,ye,dt_years,
                                mapping); CHKERRQ(ierr);

    ierr = doneWithIceInfo(variables); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "... done\n"); CHKERRQ(ierr);
    
    delete lic;
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}

