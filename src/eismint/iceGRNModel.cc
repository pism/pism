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


EISGREENAtmosCoupler::EISGREENAtmosCoupler() : PISMSnowModelAtmosCoupler() {
  doGWL3 = PETSC_FALSE;
  startGWL3Year = 0.0;
}


PetscErrorCode EISGREENAtmosCoupler::startGWL3AtYear(PetscScalar year) {
  doGWL3 = PETSC_TRUE;
  startGWL3Year = year;
  return 0;
}

PetscErrorCode EISGREENAtmosCoupler::initFromOptions(IceGrid* g) {
  PetscErrorCode ierr;
  ierr = PISMSnowModelAtmosCoupler::initFromOptions(g); CHKERRQ(ierr);
  EISGREENMassBalance* eisg_scheme = dynamic_cast<EISGREENMassBalance*>(mbscheme);
  if (eisg_scheme == NULL) {
    ierr = verbPrintf(1, g->com, 
      "EISGREENAtmosCoupler ERROR:  attachment of EISGREENMassBalance failed ... ending ...\n"
      ); CHKERRQ(ierr);
    PetscEnd();
  }
  ierr = verbPrintf(2, g->com, 
      "  special climate coupler for EISMINT-Greenland\n"
      "    -- non-default snow and ice surface temperature parameterizations\n"
      "    -- non-default interpretation of PDD factors\n"
      "    -- can add greenhouse warming if -gwl3 chosen\n"
      ); CHKERRQ(ierr);
  return 0;
}


//! Used for both ice surface temperature (boundary condition) and in snow temperature yearly cycle.
PetscScalar EISGREENAtmosCoupler::meanAnnualTemp(PetscScalar h, PetscScalar lat) {
  PetscScalar Z = PetscMax(h, 20 * (lat - 65));
  return  49.13 - 0.007992 * Z - 0.7576 * (lat) + 273.15;
}


//! Compute the temperature shift (increase) appropriate to global warming experiment.
PetscScalar EISGREENAtmosCoupler::shiftForGWL3(PetscScalar t_years, PetscScalar dt_years) {
  // compute age back to start of GWL3; use midpoint of interval as the time
  PetscScalar age = (t_years + 0.5 * dt_years) - startGWL3Year;
  if (age <= 0.0) {
    return 0.0; // before time 0.0, no warming
  } else {
    if (age <= 80.0) {
      return age * 0.035;
    } else if (age <= 500.0) {
      return 2.8 + (age - 80.0) * 0.0017;
    } else { // after 500 years, constant amount
      return 3.514;
    }
  }
}


//! EISMINT-Greenland snow-surface temperature parameterization for use in PDD.
/*! 
Formulas from \ref RitzEISMINT.  Overwrites method of same name in 
PISMSnowModelAtmosCoupler.
*/
PetscErrorCode EISGREENAtmosCoupler::parameterizedUpdateSnowSurfaceTemp(
                       PetscScalar t_years, PetscScalar dt_years, IceInfoNeededByCoupler* info) {
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid->com, 
      "entering EISGREENAtmosCoupler::parameterizedUpdateSnowTemp()\n"); CHKERRQ(ierr);
  PetscScalar **lat, **h, **T_ma, **T_mj;
  if (info == PETSC_NULL) { SETERRQ(1,"info == PETSC_NULL"); }
  ierr = info->surfelev->get_array(h);   CHKERRQ(ierr);
  ierr = info->lat->get_array(lat); CHKERRQ(ierr);
  ierr = vsnowtemp_ma.get_array(T_ma);  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.get_array(T_mj);  CHKERRQ(ierr);
  for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
      // T_ma, T_mj in K; T_ma is mean annual, T_mj is summer peak;
      //   note coupler builds sine wave yearly cycle from these
      T_ma[i][j] = meanAnnualTemp(h[i][j], lat[i][j]);
      T_mj[i][j] = 273.15 + 30.38 - 0.006277 * h[i][j] - 0.3262 * lat[i][j];
    }
  }  
  ierr = info->surfelev->end_access();   CHKERRQ(ierr);
  ierr = info->lat->end_access(); CHKERRQ(ierr);
  ierr = vsnowtemp_ma.end_access();  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.end_access();  CHKERRQ(ierr);

  if (doGWL3 == PETSC_TRUE) {
    ierr = vsnowtemp_ma.shift(shiftForGWL3(t_years,dt_years)); CHKERRQ(ierr);
    ierr = vsnowtemp_mj.shift(shiftForGWL3(t_years,dt_years)); CHKERRQ(ierr);
  }
  return 0;
}


//! EISMINT-Greenland ice-surface temperature parameterization for use in conservation of energy equation.
PetscErrorCode EISGREENAtmosCoupler::updateSurfTempAndProvide(
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  IceInfoNeededByCoupler* info, IceModelVec2* &pvst) {
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid->com, 
      "entering EISGREENAtmosCoupler::updateSurfTempAndProvide();\n"
      "  computing ice surface temperature by elevation,latitude-dependent\n"
      "  EISMINT-Greenland formulas ...\n"); CHKERRQ(ierr);

  ierr = PISMAtmosphereCoupler::updateSurfTempAndProvide(
              t_years, dt_years, info, pvst); CHKERRQ(ierr);

  if (info == PETSC_NULL) { SETERRQ(1,"info == PETSC_NULL"); }
  PetscScalar **Ts, **lat, **h;
  ierr = vsurftemp.get_array(Ts); CHKERRQ(ierr);
  ierr = info->lat->get_array(lat); CHKERRQ(ierr);
  ierr = info->surfelev->get_array(h); CHKERRQ(ierr);
  for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
      Ts[i][j] = meanAnnualTemp(h[i][j], lat[i][j]);
    }
  }
  ierr = vsurftemp.end_access(); CHKERRQ(ierr);
  ierr = info->lat->end_access(); CHKERRQ(ierr);
  ierr = info->surfelev->end_access(); CHKERRQ(ierr);

  if (doGWL3 == PETSC_TRUE) {
    ierr = vsurftemp.shift(shiftForGWL3(t_years,dt_years)); CHKERRQ(ierr);
  }
  return 0;
}


EISGREENMassBalance::EISGREENMassBalance() : PDDMassBalance() {
  // ignor configuration and just set EISMINT-Greenland defaults; user may override anyway
  pddFactorIce    = 0.008;
  pddFactorSnow   = 0.003;
  pddRefreezeFrac = 0.6;
  pddStdDev       = 5.0;

  // degree-day factors in \ref RitzEISMINT are water-equivalent
  //   thickness per degree day; ice-equivalent thickness melted per degree
  //   day is slightly larger; for example, iwfactor = 1000/910
  const PetscScalar iwfactor = config.get("fresh_water_rho") / config.get("ice_rho");
  pddFactorSnow *= iwfactor;
  pddFactorIce  *= iwfactor;
};


PetscErrorCode EISGREENMassBalance::setDegreeDayFactorsFromSpecialInfo(
                                  PetscScalar latitude, PetscScalar T_mj) {
  //PetscPrintf(PETSC_COMM_WORLD, 
  //    "EISGREENMassBalance::setDegreeDayFactorsFromSpecialInfo();\n"
  //    "  pddFactorIce,pddFactorSnow = %f,%f\n",pddFactorIce,pddFactorSnow);
  return 0;
}


PetscErrorCode IceGRNModel::setFromOptions() {
  PetscErrorCode ierr;
  PetscTruth ssl2Set, ssl3Set, ccl3Set, gwl3Set;

  expernum = 1;  // SSL2 is the default
  ierr = check_option("-ssl2", ssl2Set); CHKERRQ(ierr);
  if (ssl2Set == PETSC_TRUE)   expernum = 1;
  ierr = check_option("-ccl3", ccl3Set); CHKERRQ(ierr);
  if (ccl3Set == PETSC_TRUE)   expernum = 3;
  ierr = check_option("-gwl3", gwl3Set); CHKERRQ(ierr);
  if (gwl3Set == PETSC_TRUE)   expernum = 4;

  ierr = check_option("-ssl3", ssl3Set); CHKERRQ(ierr);
  if (ssl3Set == PETSC_TRUE) {
    ierr = PetscPrintf(grid.com,
       "experiment SSL3 (-ssl3) is not implemented\n"
       "  (choose parameters yourself, by runtime options)\n"); CHKERRQ(ierr);
    PetscEnd();
  }

  ierr = verbPrintf(2, grid.com, 
     "  setting flags equivalent to '-e 3 -ocean_kill'; user options will override\n");
     CHKERRQ(ierr);
  enhancementFactor = 3;  
  doOceanKill = PETSC_TRUE;

  if (expernum != 1) { 
    // no bed deformation for steady state (SSL2)
    // otherwise use Lingle-Clark bed deformation model for CCL3 and GWL3
    ierr = verbPrintf(2, grid.com, 
      "  setting flags equivalent to: '-bed_def_lc'; user options will override\n"); CHKERRQ(ierr);
    doBedDef = PETSC_TRUE;
    doBedIso = PETSC_FALSE;
  }

  muSliding = 0.0;  // no SIA-type sliding!

  // these flags turn off parts of the EISMINT-Greenland specification;
  //   use when extra/different data is available
  ierr = check_option("-have_bheatflx", haveGeothermalFlux);
     CHKERRQ(ierr);
  ierr = check_option("-no_EI_delete", noEllesmereIcelandDelete);
     CHKERRQ(ierr);
  
  // note: user value for -e, and -gk, and so on, will override settings above
  ierr = IceModel::setFromOptions(); CHKERRQ(ierr);  
    
  return 0;
}


PetscErrorCode IceGRNModel::init_couplers() {
  PetscErrorCode ierr;

  EISGREENAtmosCoupler* pddPCC = dynamic_cast<EISGREENAtmosCoupler*>(atmosPCC);
  if (pddPCC == PETSC_NULL) { SETERRQ(1,"dynamic_cast fail; pddPCC == PETSC_NULL\n"); }

  // the coupler (*pddPCC) will delete this at the end
  EISGREENMassBalance *eisgreen_mbscheme = new EISGREENMassBalance;
  ierr = pddPCC->setLMBScheme(eisgreen_mbscheme); CHKERRQ(ierr);

  PetscTruth gwl3StartSet;
  PetscReal  gwl3StartYear = grid.year;
  ierr = PetscOptionsGetReal(PETSC_NULL, "-gwl3_start_year", &gwl3StartYear, &gwl3StartSet);
    CHKERRQ(ierr);
  if (expernum == 4) {  // do GWL3, and set start year
    ierr = pddPCC->startGWL3AtYear(gwl3StartYear); CHKERRQ(ierr);
  } else if (gwl3StartSet == PETSC_TRUE) {
    ierr = verbPrintf(1, grid.com, 
       "WARNING: -gwl3_start_year option ignored because expernum != 4 (-gwl3 not set?).\n");
       CHKERRQ(ierr);
  }

  // only report initialization on PISMOceanClimateCoupler if allowing ice shelves
  oceanPCC->reportInitializationToStdOut = (doOceanKill == PETSC_FALSE);

  ierr = IceModel::init_couplers(); CHKERRQ(ierr);

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

