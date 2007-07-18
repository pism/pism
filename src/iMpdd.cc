// Copyright (C) 2007 Nathan Shemonski and Ed Bueler
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

#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <petscda.h>
#include "iceModel.hh"

// implementation of default positive degree day model for melting
// note defaults "DEFAULT_PDD_..." are defined in iMdefaults.cc

// FIXME:
//   3. for EIMINT-GREENLAND, summer_warming = Ts - Ta is *not* constant;
//      need alternative determined from flag?

PetscErrorCode IceModel::initPDDFromOptions() {
  PetscErrorCode  ierr;
  PetscTruth      pddSet, pddFactorSnowSet, pddFactorIceSet, 
                  pddRefreezeFracSet, pddStdDevSet, pddSummerWarmingSet,
                  pddRepeatableSet;

  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd", &pddSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_factor_snow", &pddFactorSnow, 
                               &pddFactorSnowSet); CHKERRQ(ierr);
  if (pddFactorSnowSet == PETSC_FALSE)    pddFactorSnow = DEFAULT_PDD_FACTOR_SNOW;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_factor_ice", &pddFactorIce, 
                               &pddFactorIceSet); CHKERRQ(ierr);
  if (pddFactorIceSet == PETSC_FALSE)    pddFactorIce = DEFAULT_PDD_FACTOR_ICE;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_refreeze", &pddRefreezeFrac, 
                               &pddRefreezeFracSet); CHKERRQ(ierr);
  if (pddRefreezeFracSet == PETSC_FALSE)    pddRefreezeFrac = DEFAULT_PDD_REFREEZE_FRAC;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_summer_warming", &pddSummerWarming,
                               &pddSummerWarmingSet); CHKERRQ(ierr);
  if (pddSummerWarmingSet == PETSC_FALSE)    pddSummerWarming = DEFAULT_PDD_SUMMER_WARMING;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_std_dev", &pddStdDev, &pddStdDevSet); CHKERRQ(ierr);
  if (pddStdDevSet == PETSC_FALSE)    pddStdDev = DEFAULT_PDD_STD_DEV;
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_repeatable", &pddRepeatableSet); CHKERRQ(ierr);

  // if any of the above options are turned on, then a PDD model will be used
  if ((doPDD == PETSC_TRUE) ||
      (pddSet == PETSC_TRUE) || (pddFactorSnowSet == PETSC_TRUE) || 
      (pddFactorIceSet == PETSC_TRUE) || (pddRefreezeFracSet == PETSC_TRUE) ||
      (pddStdDevSet == PETSC_TRUE) || (pddRepeatableSet == PETSC_TRUE) ||
      (pddSummerWarmingSet == PETSC_TRUE)                                       ) {
    doPDD = PETSC_TRUE;
  } else
    doPDD = PETSC_FALSE;
 
  if (doPDD == PETSC_TRUE) {
    ierr = VecDuplicate(vh, &vAccumSnow); CHKERRQ(ierr);
    ierr = VecCopy(vAccum, vAccumSnow); CHKERRQ(ierr);  // so assume vAccum *is* snow
    // initialize the random number generator: use GSL's recommended default random
    // number generator, which seems to be "mt19937" and is DIEHARD
    pddRandGen = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(pddRandGen, pddRepeatableSet ? 0 : time(0));  // number of seconds since 1/1/1970 (?)
    pddStuffCreated = PETSC_TRUE;
  } else {  
    pddStuffCreated = PETSC_FALSE;
  }
  return 0;
}


bool IceModel::IsIntegralYearPDD() {

  if (doPDD == PETSC_TRUE) {  // for PDD, go to next integral year
    if (grid.p->year < 0) {
      maxdt_temporary = fmod(fabs(grid.p->year), 1.0) * secpera;
    } else {
      maxdt_temporary = (1.0 - fmod(fabs(grid.p->year), 1.0)) * secpera;
    }
  }
  return ((fmod(grid.p->year, 1.0) + 5e-6) < 1e-5); // always return whether current is integral
}


PetscErrorCode IceModel::putBackSnowAccumPDD() {
  PetscErrorCode ierr;
  ierr = VecCopy(vAccumSnow, vAccum); CHKERRQ(ierr);  // so when vAccum is written it is snow accum
  return 0;
}


PetscErrorCode IceModel::PDDCleanup() {
  PetscErrorCode ierr;
  if (pddStuffCreated == PETSC_TRUE) {
    ierr = VecDestroy(vAccumSnow); CHKERRQ(ierr);
    gsl_rng_free(pddRandGen);
    pddStuffCreated = PETSC_FALSE;
  }  
  return 0;
}


PetscScalar IceModel::getSummerWarming(
       const PetscScalar elevation, const PetscScalar latitude, const PetscScalar Ta) const {
  // THIS PROCEDURE IS virtual AND REPLACEABLE
  // version here ignors elevation, latitude, and mean annual temperature (Ta)
  // and instead uses -pdd_summer_warming setable constant
  // see IceGRNModel for alternate implementation
  return pddSummerWarming;
}


PetscScalar IceModel::getTemperatureFromYearlyCycle(
       const PetscScalar summer_warming, const PetscScalar Ta, const int day) const {
  // THIS PROCEDURE IS virtual AND REPLACEABLE
  // NOTE: Ta = mean annual temperature (at location on surface) in degrees C
  //       summer_warming = difference between peak summer temperature and Ta
  //       day = day on 365.0 day calendar
  const int          summer_offset  = 243; // approximately August 1st
  const PetscScalar  rad_per_day    = 2.0 * PETSC_PI / 365.0;  // treat year as exactly 365 days
  const PetscScalar  days_from_summer = (PetscScalar) (day - summer_offset);
  return summer_warming * cos(rad_per_day * days_from_summer) + Ta;
}


PetscErrorCode IceModel::updateNetAccumFromPDD() {
  PetscErrorCode ierr;
  PetscScalar **accum, **Ts, **snow_accum, **h, **lat;

  if (doPDD == PETSC_FALSE)  return 0;

  // Calculate and store a random number for each day in the year. 
  // The same numbers will be used for all points on the grid
  // Note option -repeatable_pdd sets pdd_rand_gen to fixed value at start of run
  PetscScalar rand_values[365];
  if (pddStdDev > 0.0) {
    for (int day=0; day<365; day++) {
      rand_values[day] = (PetscScalar) gsl_ran_gaussian(pddRandGen, pddStdDev);
    }
  } else { // no need to waste time generating random numbers
    for (int day=0; day<365; day++) {
      rand_values[day] = 0.0;
    }
  }
  
//  verbPrintf(1,grid.com," rand_values[0..3] = %f %f %f %f\n",
//               rand_values[0],rand_values[1],rand_values[2],rand_values[3]);
               
  // compute net accumulation vAccum using surface temperature vTs and PDD model
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vAccumSnow, &snow_accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
  PetscScalar temp, summer_warming;
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      PetscScalar pdd_sum = 0.0;  // units of day^-1 K-1
      for (PetscInt day=0; day<365; day++){ // compute # of pos deg day
        summer_warming = getSummerWarming(h[i][j],lat[i][j],Ts[i][j]-ice.meltingTemp);
        temp = getTemperatureFromYearlyCycle(summer_warming,Ts[i][j]-ice.meltingTemp,day);
        temp += rand_values[day];
        if (temp > 0.0) {
          pdd_sum += temp;
        }
      }
      const PetscScalar snow = snow_accum[i][j] * secpera;  // units of m of ice-equivalent (in a year)
      const PetscScalar snow_melted = pdd_sum * pddFactorSnow;  // ditto
      if (snow_melted <= snow) {  // note this case is never active if snow_accum[i][j] < 0; e.g. rain
        accum[i][j] = ((snow - snow_melted) + (snow_melted * pddRefreezeFrac))/secpera;
      } else {
//        accum[i][j] = snow_accum[i][j] * DEFAULT_PDD_REFREEZE_FACTOR;
//        ice_melt = ((snow_melted - snow) / DEFAULT_PDD_FACTOR_SNOW) * DEFAULT_PDD_FACTOR_ICE;
        const PetscScalar excess_pdd_sum = pdd_sum - (snow / pddFactorSnow);
        const PetscScalar ice_melted = excess_pdd_sum * pddFactorIce;
        accum[i][j] = snow_accum[i][j] - (ice_melted / secpera);
      }

/*
      if ((i == id) && (j == jd)) {
        verbPrintf(1,grid.com,
              "at id,jd:  h = %f, Ts = %f, pdd_sum = %f, snow_accum = %f, accum = %f\n",
              h[id][jd],Ts[id][jd],pdd_sum,snow_accum[id][jd]*secpera,accum[id][jd]*secpera);
      }
*/
    }
  }
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vAccumSnow, &snow_accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
  return 0;
}

