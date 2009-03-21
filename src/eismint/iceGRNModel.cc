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


PISMEISGREENPDDCoupler::PISMEISGREENPDDCoupler() : PISMPDDCoupler() {
  doGWL3 = PETSC_FALSE;
  startGWL3Year = 0.0;
}


PetscErrorCode PISMEISGREENPDDCoupler::startGWL3AtYear(PetscScalar year) {
  doGWL3 = PETSC_TRUE;
  startGWL3Year = year;
  return 0;
}


PetscScalar PISMEISGREENPDDCoupler::getSummerWarming(
     const PetscScalar elevation, const PetscScalar latitude, const PetscScalar Tma) {
  // EISMINT-Greenland summer surface temperature model (expressed as
  // warming above mean annual); Tma,Ts in K; Tma is mean annual, Ts is summer peak
  const PetscScalar Ts = 273.15 + 30.38 - 0.006277 * elevation - 0.3262 * latitude;
  return Ts - Tma;
}


//! Used by updateSurfTempAndProvide().  Returns value Tma in K.
PetscScalar PISMEISGREENPDDCoupler::calculateMeanAnnual(PetscScalar h, PetscScalar lat) {
  // following EISMINT-Greenland formulas
  PetscScalar Z = PetscMax(h, 20 * (lat - 65));
  return  49.13 - 0.007992 * Z - 0.7576 * (lat) + 273.15;
}


//! Updates forcing and provides access to vsurftemp.  Updates vsurftemp, the mean annual surface temperature, according EISMINT-Greenland choices.
PetscErrorCode PISMEISGREENPDDCoupler::updateSurfTempAndProvide(
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  void *iceInfoNeeded, IceModelVec2* &pvst) {
  PetscErrorCode ierr;

  ierr = PISMPDDCoupler::updateSurfTempAndProvide(
              t_years, dt_years, iceInfoNeeded, pvst); CHKERRQ(ierr);

  if (initialize_vsurftemp_FromFile == PETSC_FALSE) {
    // update mean annual temp; note getSummerWarming() is only used in
    //    PISMPDDCoupler::updateSurfMassFluxAndProvide()
    ierr = verbPrintf(4, grid->com, 
      "computing surface temperatures by elevation,latitude-dependent EISMINT-Greenland formulas ...\n");
      CHKERRQ(ierr);
    IceInfoNeededByAtmosphereCoupler* info = (IceInfoNeededByAtmosphereCoupler*) iceInfoNeeded;
    if (info == PETSC_NULL) { SETERRQ(1,"info == PETSC_NULL"); }
    PetscScalar **Ts, **lat, **h;
    ierr = vsurftemp.get_array(Ts); CHKERRQ(ierr);
    ierr = info->lat->get_array(lat); CHKERRQ(ierr);
    ierr = info->surfelev->get_array(h); CHKERRQ(ierr);
    for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
      for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
        Ts[i][j] = calculateMeanAnnual(h[i][j], lat[i][j]);
      }
    }
    ierr = vsurftemp.end_access(); CHKERRQ(ierr);
    ierr = info->lat->end_access(); CHKERRQ(ierr);
    ierr = info->surfelev->end_access(); CHKERRQ(ierr);

    if (doGWL3 == PETSC_TRUE) {
      // for GWL3 apply global warming temperature forcing
      // compute age back to start of GWL3 by midpoint rule
      PetscScalar age = 0.5 * (2.0 * t_years + dt_years) - startGWL3Year;
      if (age > 0.0) {
        PetscScalar temp_increase;
        if (age <= 80.0) {
          temp_increase = age * 0.035;
        } else if (age <= 500.0) {
          temp_increase = 2.8 + (age - 80.0) * 0.0017;
        } else { // after 500 years
          temp_increase = 3.514;
        }
        ierr = verbPrintf(4, grid->com, 
          "   adding global warming GWL3 temperature increase of %8.4f deg; age = %8.4f; startGWL3Year = %8.4f; ...\n",
          temp_increase, age, startGWL3Year); CHKERRQ(ierr);
        ierr = vsurftemp.shift(temp_increase); CHKERRQ(ierr);
      }
    }
  }

  return 0;
}


IceGRNModel::IceGRNModel(IceGrid &g) : IceModel(g) {
  pddPCC = PETSC_NULL;
}


PetscErrorCode IceGRNModel::attachEISGREENPDDPCC(PISMEISGREENPDDCoupler &p) {
  pddPCC = &p;
  return 0;
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



PetscErrorCode IceGRNModel::init_couplers() {
  PetscErrorCode ierr;

  pddPCC = (PISMEISGREENPDDCoupler*) atmosPCC;
  if (pddPCC == PETSC_NULL) { SETERRQ(1,"pddPCC == PETSC_NULL\n"); }
  
  pddPCC->initialize_vsnowaccum_FromFile = true;

  pddPCC->initialize_vsurftemp_FromFile = (haveSurfaceTemp == PETSC_TRUE);

  PetscTruth gwl3StartSet;
  PetscReal  gwl3StartYear = grid.year;
  ierr = PetscOptionsGetReal(PETSC_NULL, "-gwl3_start_year", &gwl3StartYear, &gwl3StartSet); CHKERRQ(ierr);
  if (expernum == 4) {  // do GWL3, and set start year
    ierr = pddPCC->startGWL3AtYear(gwl3StartYear); CHKERRQ(ierr);
  } else if (gwl3StartSet == PETSC_TRUE) {
    ierr = verbPrintf(1, grid.com, 
       "WARNING: -gwl3_start_year option ignored because expernum != 4 (-gwl3 not set?).\n"); CHKERRQ(ierr);
  }

  // only report initialization on PISMOceanClimateCoupler if allowing ice shelves
  oceanPCC->reportInitializationToStdOut = (doOceanKill == PETSC_FALSE);

  ierr = IceModel::init_couplers(); CHKERRQ(ierr);

  // PDD is already set by atmosPCC->initFromOptions() in IceModel::init_couplers()
  //   here we just warn if nondefault values are used
  PetscTruth pddSummerWarmingSet, pddStdDevSet;
  // FIXME: exposed to -pdd_summer_warming 0 bug?
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_summer_warming",
              &pddSummerWarmingSet); CHKERRQ(ierr);
  if (pddSummerWarmingSet == PETSC_TRUE) { 
    ierr = verbPrintf(1, grid.com, 
       "WARNING: -pdd_summer_warming option ignored.\n"
       "  Using EISMINT-GREENLAND summer temperature formula\n"); CHKERRQ(ierr);
  }

  // FIXME:  NOT TRUE:  -pdd_std_dev 0 does not make sense, so using PetscOptionsHasName is OK
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_std_dev", &pddStdDevSet); 
     CHKERRQ(ierr);
  if (pddStdDevSet == PETSC_FALSE) {
    pddPCC->pddStdDev = 5.0;  // EISMINT-GREENLAND default; user may override
  }

  return 0;
}


PetscErrorCode IceGRNModel::set_vars_from_options() {
  PetscErrorCode ierr;

  const PetscScalar EISMINT_G_geothermal = 0.050;      // J/m^2 s; geothermal flux

  // Let the base class handle bootstrapping:
  ierr = IceModel::set_vars_from_options(); CHKERRQ(ierr);

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
    // init couplers needs to have already occurred, but it has (see IceModel::init()):
    if (pddPCC == PETSC_NULL) { SETERRQ(1,"pddPCC == PETSC_NULL\n"); }
    IceModelVec2* dummy;
    ierr = pddPCC->updateSurfTempAndProvide(0.0, 0.0, &info_atmoscoupler, dummy); CHKERRQ(ierr);
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

  return 0;
}


PetscErrorCode IceGRNModel::ellePiecewiseFunc(PetscScalar lon, PetscScalar *lat) {
  double l1_x1 = -68.18, l1_y1 = 80.1;
  double l1_x2 = -62, l1_y2 = 82.24;
  double m, b;  // piecewise boundaries

  m = (l1_y1 - l1_y2) / (l1_x1 - l1_x2);
  b = (l1_y2) - m * (l1_x2);
  *lat = m * lon + b;
  return 0;
}


PetscErrorCode IceGRNModel::cleanExtraLand(){
  PetscErrorCode ierr;
  PetscScalar lat_line;
  // make all mask points southeast of the following point into FLOATING_OCEAN0
  double ice_lon = 30, ice_lat = 67;
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
  ierr = vMask.beginGhostComm(); CHKERRQ(ierr);
  ierr = vMask.endGhostComm(); CHKERRQ(ierr);
  return 0;
}

