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

#include "iceGRNModel.hh"

EISGREENAtmosCoupler::EISGREENAtmosCoupler() : PISMSnowModelAtmosCoupler() {
  doGreenhouse = PETSC_FALSE;
  startYearGreenhouse = 0.0;
}


PetscErrorCode EISGREENAtmosCoupler::startGreenhouseAtYear(PetscScalar year) {
  doGreenhouse = PETSC_TRUE;
  startYearGreenhouse = year;
  return 0;
}


PetscErrorCode EISGREENAtmosCoupler::initFromOptions(IceGrid* g, const PISMVars &variables) {
  PetscErrorCode ierr;
  ierr = PISMSnowModelAtmosCoupler::initFromOptions(g, variables); CHKERRQ(ierr);
  EISGREENMassBalance* eisg_scheme = dynamic_cast<EISGREENMassBalance*>(mbscheme);
  if (eisg_scheme == NULL) {
    ierr = verbPrintf(1, g->com, 
      "EISGREENAtmosCoupler ERROR:  attachment of EISGREENMassBalance failed ... ending ...\n"
      ); CHKERRQ(ierr);
    PetscEnd();
  }
  // because ice surface temperature is coming from a parameterization,
  //   the variable 'artm' in the output file is merely diagnostic
  ierr = vsurftemp.set_attr("pism_intent","climate_diagnostic"); CHKERRQ(ierr);
  // a bit of info to user
  ierr = verbPrintf(2, g->com, 
      "  special climate coupler for EISMINT-Greenland\n"
      "    -- non-default snow and ice surface temperature parameterizations\n"
      "    -- non-default interpretation of PDD factors\n"
      "    -- can add greenhouse warming if -gwl3 chosen\n"
      ); CHKERRQ(ierr);
  return 0;
}


/*!
Used for both ice surface temperature (boundary condition) and in snow temperature yearly cycle.
 */
PetscScalar EISGREENAtmosCoupler::meanAnnualTemp(PetscScalar h, PetscScalar lat) {
  PetscScalar Z = PetscMax(h, 20 * (lat - 65));
  return  49.13 - 0.007992 * Z - 0.7576 * (lat) + 273.15;
}


PetscScalar EISGREENAtmosCoupler::shiftForGreenhouse(PetscScalar t_years, PetscScalar dt_years) {
  // compute age back to start of GWL3; use midpoint of interval as the time
  PetscScalar age = (t_years + 0.5 * dt_years) - startYearGreenhouse;
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


PetscErrorCode EISGREENAtmosCoupler::parameterizedUpdateSnowSurfaceTemp(
                       PetscScalar t_years, PetscScalar dt_years) {
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid->com, 
      "entering EISGREENAtmosCoupler::parameterizedUpdateSnowTemp()\n"); CHKERRQ(ierr);
  PetscScalar **lat_degN, **h, **T_ma, **T_mj;

  ierr = surfelev->get_array(h);   CHKERRQ(ierr);
  ierr = lat->get_array(lat_degN); CHKERRQ(ierr);
  ierr = vsnowtemp_ma.get_array(T_ma);  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.get_array(T_mj);  CHKERRQ(ierr);
  for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
      // T_ma, T_mj in K; T_ma is mean annual, T_mj is summer peak;
      //   note coupler builds sine wave yearly cycle from these
      T_ma[i][j] = meanAnnualTemp(h[i][j], lat_degN[i][j]);
      T_mj[i][j] = 273.15 + 30.38 - 0.006277 * h[i][j] - 0.3262 * lat_degN[i][j];
    }
  }  
  ierr = surfelev->end_access();   CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);
  ierr = vsnowtemp_ma.end_access();  CHKERRQ(ierr);
  ierr = vsnowtemp_mj.end_access();  CHKERRQ(ierr);

  if (doGreenhouse == PETSC_TRUE) {
    const PetscScalar shift = shiftForGreenhouse(t_years,dt_years);
    ierr = vsnowtemp_ma.shift(shift); CHKERRQ(ierr);
    ierr = vsnowtemp_mj.shift(shift); CHKERRQ(ierr);
  }
  return 0;
}


PetscErrorCode EISGREENAtmosCoupler::updateSurfTempAndProvide(
                  const PetscScalar t_years, const PetscScalar dt_years, 
                  IceModelVec2* &pvst) {
  PetscErrorCode ierr;
  ierr = verbPrintf(4, grid->com, 
      "entering EISGREENAtmosCoupler::updateSurfTempAndProvide();\n"
      "  computing ice surface temperature by elevation,latitude-dependent\n"
      "  EISMINT-Greenland formulas ...\n"); CHKERRQ(ierr);

  ierr = PISMAtmosphereCoupler::updateSurfTempAndProvide(
              t_years, dt_years, pvst); CHKERRQ(ierr);

  PetscScalar **Ts, **lat_degN, **h;
  ierr = vsurftemp.get_array(Ts); CHKERRQ(ierr);
  ierr = lat->get_array(lat_degN); CHKERRQ(ierr);
  ierr = surfelev->get_array(h); CHKERRQ(ierr);
  for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
      Ts[i][j] = meanAnnualTemp(h[i][j], lat_degN[i][j]);
    }
  }
  ierr = vsurftemp.end_access(); CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);
  ierr = surfelev->end_access(); CHKERRQ(ierr);

  if (doGreenhouse == PETSC_TRUE) {
    ierr = vsurftemp.shift(shiftForGreenhouse(t_years,dt_years)); CHKERRQ(ierr);
  }
  return 0;
}


EISGREENMassBalance::EISGREENMassBalance(NCConfigVariable* myconfig) : PDDMassBalance(myconfig) {
  // ignor configuration and just set EISMINT-Greenland defaults; user may override anyway
  pddFactorIce    = 0.008;
  pddFactorSnow   = 0.003;
  pddRefreezeFrac = 0.6;
  pddStdDev       = 5.0;

  // degree-day factors in \ref RitzEISMINT are water-equivalent
  //   thickness per degree day; ice-equivalent thickness melted per degree
  //   day is slightly larger; for example, iwfactor = 1000/910
  const PetscScalar iwfactor = config->get("fresh_water_density") / config->get("ice_density");
  pddFactorSnow *= iwfactor;
  pddFactorIce  *= iwfactor;
}


PetscErrorCode EISGREENMassBalance::setDegreeDayFactorsFromSpecialInfo(
		PetscScalar /*latitude*/, PetscScalar /*T_mj*/) {
  // a message like this may be useful for debugging; note task #6216
  //PetscPrintf(PETSC_COMM_WORLD,"EISGREENMassBalance::setDegreeDayFactorsFromSpecialInfo();\n"
  //    "  pddFactorIce,pddFactorSnow = %f,%f\n",pddFactorIce,pddFactorSnow);
  return 0;
}


PetscErrorCode IceGRNModel::setFromOptions() {
  PetscErrorCode ierr;
  PetscTruth ssl2Set, ssl3Set, ccl3Set, gwl3Set;

  exper = SSL2;  // default
  ierr = check_option("-ssl2", ssl2Set); CHKERRQ(ierr);
  if (ssl2Set == PETSC_TRUE)   exper = SSL2;
  ierr = check_option("-ccl3", ccl3Set); CHKERRQ(ierr);
  if (ccl3Set == PETSC_TRUE)   exper = CCL3;
  ierr = check_option("-gwl3", gwl3Set); CHKERRQ(ierr);
  if (gwl3Set == PETSC_TRUE)   exper = GWL3;

  ierr = check_option("-ssl3", ssl3Set); CHKERRQ(ierr);
  if (ssl3Set == PETSC_TRUE) {
    ierr = PetscPrintf(grid.com,
       "experiment SSL3 (-ssl3) is not implemented ... ENDING\n"
       "  (choose parameters yourself, by runtime options)\n"); CHKERRQ(ierr);
    PetscEnd();
  }

  ierr = verbPrintf(2, grid.com, 
     "  setting flags equivalent to '-e 3 -ocean_kill'; user options may override ...\n");
     CHKERRQ(ierr);
  config.set("enhancement_factor", 3.0);
  config.set_flag("ocean_kill", true);

  if (exper != SSL2) { 
    // use Lingle-Clark bed deformation model for CCL3 and GWL3 but not SSL2
    ierr = verbPrintf(2, grid.com, 
      "  setting flags equivalent to: '-bed_def_lc'; user options may override ...\n");
      CHKERRQ(ierr);
    config.set_flag("do_bed_deformation", true);
    config.set_flag("do_bed_iso", false);
  }

  config.set("mu_sliding", 0.0);  // no SIA-type sliding!; see [\ref RitzEISMINT]

  // use the EISMINT-Greenland value if no value in -boot_from file
  config.set("bootstrapping_geothermal_flux_value_no_var", 0.050);
  
  // note: user value for -e, and -gk, and so on, will override settings above
  ierr = IceModel::setFromOptions(); CHKERRQ(ierr);  
    
  return 0;
}


PetscErrorCode IceGRNModel::init_couplers() {
  PetscErrorCode ierr;

  EISGREENAtmosCoupler* pddPCC = dynamic_cast<EISGREENAtmosCoupler*>(atmosPCC);
  if (pddPCC == PETSC_NULL) { SETERRQ(1,"dynamic_cast fail; pddPCC == PETSC_NULL\n"); }

  // the coupler (*pddPCC) will delete this at the end
  EISGREENMassBalance *eisgreen_mbscheme = new EISGREENMassBalance(&config);
  ierr = pddPCC->setLMBScheme(eisgreen_mbscheme); CHKERRQ(ierr);

  PetscTruth gwl3StartSet;
  PetscReal  gwl3StartYear = grid.year;
  ierr = PetscOptionsGetReal(PETSC_NULL, "-gwl3_start_year", &gwl3StartYear, &gwl3StartSet);
    CHKERRQ(ierr);
  if (exper == GWL3) {  // do GWL3, and set start year
    ierr = pddPCC->startGreenhouseAtYear(gwl3StartYear); CHKERRQ(ierr);
  } else if (gwl3StartSet == PETSC_TRUE) {
    ierr = verbPrintf(1, grid.com, 
       "WARNING: -gwl3_start_year option ignored;  experiment != GWL3;  option -gwl3 not set?\n");
       CHKERRQ(ierr);
  }

  // only report initialization on PISMOceanClimateCoupler if allowing ice shelves
  bool do_ocean_kill = config.get_flag("ocean_kill");
  oceanPCC->reportInitializationToStdOut = (do_ocean_kill == PETSC_FALSE);

  ierr = IceModel::init_couplers(); CHKERRQ(ierr);

  return 0;
}

