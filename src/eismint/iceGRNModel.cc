// Copyright (C) 2007-2009 Nathan Shemonski, Ed Bueler and Constantine Khroulev
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

#include <petscda.h>
#include <cstring>
#include <netcdf.h>
#include "../base/nc_util.hh"
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"
#include "../coupler/forcing.hh"
#include "../coupler/pccoupler.hh"
#include "iceGRNModel.hh"


const PetscScalar EISMINT_G_geothermal = 0.050;      // J/m^2 s; geothermal flux


PISMEISGREENPDDCoupler::PISMEISGREENPDDCoupler() : PISMPDDCoupler() {
}


PetscScalar PISMEISGREENPDDCoupler::getSummerWarming(
     const PetscScalar elevation, const PetscScalar latitude, const PetscScalar Tma) {
  // this is virtual in IceModel
  // EISMINT-Greenland summer surface temperature model (expressed as

  //OLD:  // warming above mean annual); Tma,Ts in deg C; Tma is mean annual, Ts is summer peak
  //OLD:  const PetscScalar Ts = 30.38 - 0.006277 * elevation - 0.3262 * latitude;

  // warming above mean annual); Tma,Ts in K; Tma is mean annual, Ts is summer peak
  const PetscScalar Ts = 273.15 + 30.38 - 0.006277 * elevation - 0.3262 * latitude;
  return Ts - Tma;
}


IceGRNModel::IceGRNModel(IceGrid &g, IceType *i) : IceModel(g, i) {
  // only call parent's constructor; do all classes need constructors?
  pddPCC = PETSC_NULL;
}


PetscErrorCode IceGRNModel::setFromOptions() {
  PetscErrorCode ierr;
  PetscTruth ssl2Set, ssl3Set, ccl3Set, gwl3Set;

  expernum = 1;  // SSL2 is the default
  ierr = PetscOptionsHasName(PETSC_NULL, "-ssl2", &ssl2Set); CHKERRQ(ierr);
  if (ssl2Set == PETSC_TRUE)   expernum = 1;
  ierr = PetscOptionsHasName(PETSC_NULL, "-ccl3", &ccl3Set); CHKERRQ(ierr);
  if (ccl3Set == PETSC_TRUE)   expernum = 3;
  ierr = PetscOptionsHasName(PETSC_NULL, "-gwl3", &gwl3Set); CHKERRQ(ierr);
  if (gwl3Set == PETSC_TRUE)   expernum = 4;

  ierr = PetscOptionsHasName(PETSC_NULL, "-ssl3", &ssl3Set); CHKERRQ(ierr);
  if (ssl3Set == PETSC_TRUE) {
    ierr = PetscPrintf(grid.com,
       "EISMINT-Greenland experiment SSL3 (-ssl3) is not implemented.\n"
       "Choose parameters yourself, by runtime options.\n"); CHKERRQ(ierr);
    PetscEnd();
  }

  ierr = verbPrintf(2, grid.com, 
     "EISMINT-Greenland mode (IceGRNModel) setting flags equivalent to:\n"
     "  '-e 3 -ocean_kill -skip 20', but user options will override\n"); CHKERRQ(ierr);
  enhancementFactor = 3;  
  doOceanKill = PETSC_TRUE;
  doSkip = PETSC_TRUE;
  skipMax = 20;

  if (expernum == 1) { // no bed deformation for steady state (SSL2)
    doBedDef = PETSC_FALSE;
  } else { // use Lingle-Clark bed deformation model for CCL3 and GWL3
    doBedDef = PETSC_TRUE;
    doBedIso = PETSC_FALSE;
  }

  muSliding = 0.0;  // no SIA-type sliding!

  // these flags turn off parts of the EISMINT-Greenland specification;
  //   use when extra/different data is available
  ierr = PetscOptionsHasName(PETSC_NULL, "-have_artm", &haveSurfaceTemp);
     CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-have_bheatflx", &haveGeothermalFlux);
     CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-no_EI_delete", &noEllesmereIcelandDelete);
     CHKERRQ(ierr);
  
  // note: user value for -e, and -gk, and so on, will override settings above
  ierr = IceModel::setFromOptions(); CHKERRQ(ierr);  
    
  return 0;
}


PetscErrorCode IceGRNModel::attachEISGREENPDDPCC(PISMEISGREENPDDCoupler &p) {
  pddPCC = &p;
  return 0;
}


PetscErrorCode IceGRNModel::initFromOptions(PetscTruth doHook) {
  PetscErrorCode ierr;
  char inFile[PETSC_MAX_PATH_LEN];
  PetscTruth inFileSet, bootFileSet;

  pddPCC = (PISMEISGREENPDDCoupler*) atmosPCC;
  
  pddPCC->initialize_vsnowaccum_FromFile = true;

  if (haveSurfaceTemp == PETSC_FALSE) { // default case: EISMINT-Greenland provides formulas
    pddPCC->initialize_vsurftemp_FromFile = false;
  }

  if (doOceanKill == PETSC_TRUE) {
    oceanPCC->reportInitializationToStdOut = false;
  }
  
  ierr = IceModel::initFromOptions(); CHKERRQ(ierr);  

  ierr = PetscOptionsGetString(PETSC_NULL, "-i", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-boot_from", inFile,
                               PETSC_MAX_PATH_LEN, &bootFileSet); CHKERRQ(ierr);

  // PDD is already set by atmosPCC->initFromOptions() in IceModel::initFromOptions();
  //   here we just warn if nondefault values are used
  PetscTruth pddSummerWarmingSet, pddStdDevSet;
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_summer_warming",
              &pddSummerWarmingSet); CHKERRQ(ierr);
  if (pddSummerWarmingSet == PETSC_TRUE) { 
    ierr = verbPrintf(1, grid.com, 
       "WARNING: -pdd_summer_warming option ignored.\n"
       "  Using EISMINT-GREENLAND summer temperature formula\n"); CHKERRQ(ierr);
  }
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_std_dev", &pddStdDevSet); 
     CHKERRQ(ierr);
  if (pddStdDevSet == PETSC_FALSE) {
    pddPCC->pddStdDev = 5.0;  // EISMINT-GREENLAND default; user may override
  }
  
  // set up surface temperature, geothermal flux, and Ellesmere/Iceland delete
  if (inFileSet == PETSC_TRUE) {
    if (bootFileSet) {
      ierr = verbPrintf(1, grid.com, "WARNING: both -boot_from and -i given; using -i\n");
         CHKERRQ(ierr);
    }
  } else if (bootFileSet == PETSC_TRUE) {
    // though default bootstrapping has set the new temperatures, we usually need to set 
    // the surface temp and geothermal flux at base and then set 3D temps again
    if (haveGeothermalFlux == PETSC_FALSE) {
      ierr = verbPrintf(2, grid.com,
         "geothermal flux set to EISMINT-Greenland value %f W/m^2\n",
         EISMINT_G_geothermal); CHKERRQ(ierr);
      ierr = vGhf.set(EISMINT_G_geothermal); CHKERRQ(ierr);
    }
    if (haveSurfaceTemp == PETSC_FALSE) {
      ierr = verbPrintf(2, grid.com, 
         "computing surface temps by EISMINT-Greenland elevation-latitude rule\n");
         CHKERRQ(ierr);
      ierr = updateTs(); CHKERRQ(ierr);
    }
    if ((haveGeothermalFlux == PETSC_FALSE) || (haveSurfaceTemp == PETSC_FALSE)) {
      ierr = verbPrintf(2, grid.com, 
         "filling in temperatures AGAIN at depth using quartic guess (for EISMINT-Greenland)\n");
         CHKERRQ(ierr);
      ierr = putTempAtDepth(); CHKERRQ(ierr);
    }
    if (noEllesmereIcelandDelete == PETSC_FALSE) {
      ierr = verbPrintf(2, grid.com, 
         "removing extra land (Ellesmere and Iceland) using EISMINT-Greenland rule\n");
         CHKERRQ(ierr);
      ierr = cleanExtraLand(); CHKERRQ(ierr);
    }
  } else {
    ierr = PetscPrintf(grid.com,"ERROR: IceGRNModel needs an input file\n"); CHKERRQ(ierr);
    PetscEnd();
  }

  if (!isInitialized()) {
    SETERRQ(1, "IceGRNModel has not been initialized.\n");
  }

  return 0;
}


PetscErrorCode IceGRNModel::additionalAtStartTimestep() {
  PetscErrorCode ierr;

  // for all experiments, at each time step we need to recompute
  // surface temperatures from surface elevation and latitude, unless the 
  // user supplies an additional map of mean annual surface temps
  if (haveSurfaceTemp == PETSC_FALSE) {
    ierr = updateTs(); CHKERRQ(ierr);
  }

  if (expernum == 4) {  // for GWL3 apply global warming temperature forcing
    PetscScalar temp_increase;
    PetscScalar age = grid.year - startYear;
    if (age <= 80.0) {
      temp_increase = age * 0.035;
    } else if (age <= 500.0) {
      temp_increase = 2.8 + (age - 80.0) * 0.0017;
    } else {
      temp_increase = 3.514;
    }
    if (pddPCC == PETSC_NULL) { SETERRQ(1,"pddPCC == PETSC_NULL"); }
    ierr = pddPCC->vsurftemp.shift(temp_increase); CHKERRQ(ierr);
  }
  return 0;
}


// Used below in IceGRNModel::updateTs().  Returns value in K.
PetscScalar IceGRNModel::calculateMeanAnnual(PetscScalar h, PetscScalar lat) {
  // EISMINT-Greenland surface temperature model
  PetscScalar Z = PetscMax(h, 20 * (lat - 65));
  return 49.13 - 0.007992 * Z - 0.7576 * (lat) + ice->meltingTemp;
}


PetscErrorCode IceGRNModel::updateTs() {
  PetscErrorCode ierr;
  PetscScalar **Ts, **lat, **h;
  
  ierr = verbPrintf(4, grid.com, 
     "recomputing surface temperatures according to EISMINT-Greenland rule"
     " and setting TsOffset=0.0\n");
     CHKERRQ(ierr);

  // PDD will use vsurftemp to determine yearly cycle from parameters
  ierr = pddPCC->vsurftemp.get_array(Ts); CHKERRQ(ierr);
  ierr = vLatitude.get_array(lat); CHKERRQ(ierr);
  ierr = vh.get_array(h); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      Ts[i][j] = calculateMeanAnnual(h[i][j], lat[i][j]);
    }
  }
  ierr = pddPCC->vsurftemp.end_access(); CHKERRQ(ierr);
  ierr = vLatitude.end_access(); CHKERRQ(ierr);
  ierr =        vh.end_access(); CHKERRQ(ierr);

//FIXME: additional kludge to deal with hosed initialization sequence
  // also set pddPCC->vsurftemp because it is used at initialization

  PetscScalar **foo;
  ierr = pddPCC->vsurftemp.get_array(Ts); CHKERRQ(ierr);
  ierr = vWork2d[0].get_array(foo); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      foo[i][j] = Ts[i][j];
    }
  }
  ierr = pddPCC->vsurftemp.end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = vWork2d[0].beginGhostComm(); CHKERRQ(ierr);
  ierr = vWork2d[0].endGhostComm(); CHKERRQ(ierr);
//FIXME: end kludge  

  return 0;
}


PetscErrorCode IceGRNModel::ellePiecewiseFunc(PetscScalar lon, PetscScalar *lat) {
  float l1_x1 = -68.18, l1_y1 = 80.1;
  float l1_x2 = -62, l1_y2 = 82.24;
  float m, b;  // piecewise boundaries

  m = (l1_y1 - l1_y2) / (l1_x1 - l1_x2);
  b = (l1_y2) - m * (l1_x2);
  *lat = m * lon + b;
  return 0;
}


PetscErrorCode IceGRNModel::cleanExtraLand(){
  PetscErrorCode ierr;
  PetscScalar lat_line;
  // make all mask points southeast of the following point into FLOATING_OCEAN0
  float ice_lon = 30, ice_lat = 67;
  PetscScalar **lat, **lon, **mask;

  ierr = vLongitude.get_array(lon); CHKERRQ(ierr);
  ierr =  vLatitude.get_array(lat); CHKERRQ(ierr);
  ierr = vMask.get_array(mask); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      ellePiecewiseFunc(lon[i][j], &lat_line);
      if (lat[i][j]>lat_line) {	// Ellesmere case
          mask[i][j] = MASK_FLOATING_OCEAN0;
      } else if (lat[i][j] < ice_lat && lon[i][j] > -ice_lon) {
        mask[i][j] = MASK_FLOATING_OCEAN0; // Iceland case
      }
    }
  } 
  ierr = vLongitude.end_access(); CHKERRQ(ierr);
  ierr =  vLatitude.end_access(); CHKERRQ(ierr);
  ierr = vMask.end_access(); CHKERRQ(ierr);
  
  // when mask is changed we must communicate the ghosted values
  //   because neighbor's mask matters
  ierr = vMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vMask.endGhostComm(); CHKERRQ(ierr);
  return 0;
}

