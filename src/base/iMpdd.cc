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
#include <gsl/gsl_sf.h>       // for erfc()
#include <petscda.h>
#include "iceModel.hh"

// implementation of default positive degree day model for melting
// note defaults "DEFAULT_PDD_..." are defined in iMdefaults.cc

// "Calov-Greve" method from R Calov and R Greve (2005), "Correspondence: A semi-analytical 
// solution for the positive degree-day model with stochastic temperature variations,"
// J. Glaciol. 51 (172), pp 173--175.

// FIXME:
//   3. for EISMINT-GREENLAND, summer_warming = Ts - Ta is *not* constant;
//      need alternative determined from flag?

PetscErrorCode IceModel::initPDDFromOptions() {
  PetscErrorCode  ierr;

  PetscTruth      pddSet, pddRandSet, pddRepeatableSet;
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd", &pddSet); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_rand", &pddRandSet); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_rand_repeatable", &pddRepeatableSet); CHKERRQ(ierr);

  // if any of the above options are turned on, then a PDD model will be used
  doPDDTrueRand = PETSC_FALSE;
  if ( (pddRandSet == PETSC_TRUE) || (pddRepeatableSet == PETSC_TRUE) ) {
    doPDDTrueRand = PETSC_TRUE;
  }
  if ( (doPDD == PETSC_TRUE) || (pddSet == PETSC_TRUE) || (doPDDTrueRand == PETSC_TRUE) ) {
    doPDD = PETSC_TRUE;
  } else {
    doPDD = PETSC_FALSE;
  }


  if (doPDD == PETSC_TRUE) {

    if (pddStuffCreated == PETSC_FALSE) {
      ierr = VecDuplicate(vh, &vAccumSnow); CHKERRQ(ierr);
      ierr = VecCopy(vAccum, vAccumSnow); CHKERRQ(ierr);  // assume vAccum *is* snow
      pddStuffCreated = PETSC_TRUE;
    } else {
      SETERRQ(1, "ERROR: initPDDFromOptions() called with pddStuffCreated == TRUE\n");
    }

    if ( (pddRandSet == PETSC_TRUE) || (pddRepeatableSet == PETSC_TRUE) ) {
      if (pddRandStuffCreated == PETSC_TRUE) {
        SETERRQ(2, 
          "ERROR: initPDDFromOptions() called with -pdd_rand[_repeatable] pddRandStuffCreated == TRUE\n");
      }
      // initialize the random number generator: use GSL's recommended default random
      // number generator, which seems to be "mt19937" and is DIEHARD
      pddRandGen = gsl_rng_alloc(gsl_rng_default);
      gsl_rng_set(pddRandGen, pddRepeatableSet ? 0 : time(0));  // number of seconds since 1/1/1970 (?)
      pddRandStuffCreated = PETSC_TRUE;
    } else {
      pddRandStuffCreated = PETSC_FALSE;
    }

    pddFactorSnow = DEFAULT_PDD_FACTOR_SNOW;
    pddFactorIce = DEFAULT_PDD_FACTOR_ICE;
    pddRefreezeFrac = DEFAULT_PDD_REFREEZE_FRAC;
    pddSummerWarming = DEFAULT_PDD_SUMMER_WARMING;
    pddStdDev = DEFAULT_PDD_STD_DEV;
    PetscTruth   pSet;
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_factor_snow", &pddFactorSnow, &pSet); CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_factor_ice", &pddFactorIce, &pSet); CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_refreeze", &pddRefreezeFrac, &pSet); CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_summer_warming", &pddSummerWarming, &pSet); CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_std_dev", &pddStdDev, &pSet); CHKERRQ(ierr);
  }
  
  return 0;
}


PetscErrorCode IceModel::setMaxdtTempPDD() {

  if (doPDDTrueRand == PETSC_TRUE) {  // for PDD based on randomness, go to next integral year
    PetscScalar   mydt;
    if (grid.p->year < 0) {
      mydt = fmod(fabs(grid.p->year), 1.0) * secpera;
    } else {
      mydt = (1.0 - fmod(fabs(grid.p->year), 1.0)) * secpera;
    }
    maxdt_temporary = (maxdt_temporary < 0.0) ? mydt : PetscMin(maxdt_temporary,mydt);
  }
  return 0;
}


bool IceModel::IsIntegralYearPDD() {

  if (doPDDTrueRand == PETSC_TRUE) {  // for PDD based on randomness,
    return ((fmod(grid.p->year, 1.0) + 5e-6) < 1e-5); // return boolean for whether current year is integral
  } else {
    return true; // if default Calov-Greve model, will always compute PDD from integral
  }
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
    pddStuffCreated = PETSC_FALSE;
  }  
  if (pddRandStuffCreated == PETSC_TRUE) {
    gsl_rng_free(pddRandGen);
    pddRandStuffCreated = PETSC_FALSE;
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
       const PetscScalar summer_warming, const PetscScalar Ta, const PetscScalar day) const {
  // THIS PROCEDURE IS virtual AND REPLACEABLE
  // NOTE: Ta = mean annual temperature (at location on surface) in degrees C
  //       summer_warming = difference between peak summer temperature and Ta
  //       day = day on 365.0 day calendar
  const PetscScalar  summer_offset  = 243.0; // approximately August 1st
  const PetscScalar  rad_per_day    = 2.0 * PETSC_PI / 365.0;  // treat year as exactly 365 days
  const PetscScalar  days_from_summer = day - summer_offset;
  return summer_warming * cos(rad_per_day * days_from_summer) + Ta;
}


double IceModel::CalovGreveIntegrand(const double T) {
  double z;
  const double sigSqr = pddStdDev * pddStdDev;
  z = (pddStdDev / sqrt(2.0 * PETSC_PI)) * exp(- T * T / (2.0 * sigSqr)) 
      + (T / 2.0) * gsl_sf_erfc(- T / (sqrt(2.0) * pddStdDev));
  return z;
}


PetscErrorCode IceModel::updateNetAccumFromPDD() {
  PetscErrorCode ierr;

  // do nothing if PDD not turned on
  if (doPDD == PETSC_FALSE)  return 0;

  // do nothing if PDD random model turned on but its not an integral year
//  if ((doPDDTrueRand == PETSC_TRUE) && (IsIntegralYearPDD() == false))  return 0;

  if (doPDDTrueRand == PETSC_TRUE) {

    //  verbPrintf(1,grid.com," true=%d, false=%d, doPDD=%d, doPDDTrueRand=%d, IsIntegralYearPDD()=%d\n",
    //             PETSC_TRUE, PETSC_FALSE, doPDD, doPDDTrueRand, IsIntegralYearPDD());

/*
    // Calculate and store a random number for each day in the year. 
    // The same numbers will be used for all points on the grid
    // Note option -pdd_rand_repeatable sets pddRandGen to fixed value at start of run
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
*/
               
    // compute net accumulation vAccum using surface temperature vTs and PDD model
    PetscScalar **accum, **Ts, **snow_accum, **h, **lat;
    ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vAccumSnow, &snow_accum); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);

    const PetscScalar   start = grid.p->year - dt / secpera;  // note grid.p->year has *end* of step
    const PetscInt      intstartday = (int) ceil(365.24 * (start - floor(start)));
    const PetscInt      num_days = (int) ceil(365.24 * (dt / secpera));

    PetscScalar temp, summer_warming;
    for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
        PetscScalar pdd_sum = 0.0;  // units of day^-1 K-1
        
//        for (PetscInt day=0; day<365; day++){ // compute # of pos deg day
        for (PetscInt day = intstartday; day < num_days; day++){ // compute # of pos deg day
          summer_warming = getSummerWarming(h[i][j],lat[i][j],Ts[i][j]-ice.meltingTemp);
          temp = getTemperatureFromYearlyCycle(summer_warming,Ts[i][j]-ice.meltingTemp,
                                               (PetscScalar) day);
          temp += (PetscScalar) gsl_ran_gaussian(pddRandGen, pddStdDev);
          //temp += rand_values[day];
          if (temp > 0.0) {
            pdd_sum += temp;
          }
        }
        //const PetscScalar dt_pdd = secpera;  // this random model only runs once per year
        const PetscScalar dt_pdd = dt;
        const PetscScalar snow = snow_accum[i][j] * dt_pdd;  // units of m of ice-equivalent (in a year)
        const PetscScalar snow_melted = pdd_sum * pddFactorSnow;  // ditto
        if (snow_melted <= snow) {  // this case is never active if snow_accum[i][j] < 0; e.g. rain
          accum[i][j] = ((snow - snow_melted) + (snow_melted * pddRefreezeFrac)) / dt_pdd;
        } else {
          const PetscScalar excess_pdd_sum = pdd_sum - (snow / pddFactorSnow);
          const PetscScalar ice_melted = excess_pdd_sum * pddFactorIce;
          accum[i][j] = - ice_melted / dt_pdd;
        }


/*
      if ((i == id) && (j == jd)) {
        verbPrintf(1,grid.com,
              "at id,jd:  h = %8.2f, Ts = %8.2f, pdd_sum = %8.2f, snow_accum = %8.4f, accum = %8.4f\n"
              "           num_days = %d, intstartday = %d,\n",
              h[id][jd],Ts[id][jd],pdd_sum,snow_accum[id][jd]*secpera,accum[id][jd]*secpera,
              num_days, intstartday);
      }
*/

      }
    }
    ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vAccumSnow, &snow_accum); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
    
  } else {
  
    // usual case: use Calov-Greve method
    PetscScalar **accum, **Ts, **snow_accum, **h, **lat;
    ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vAccumSnow, &snow_accum); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
    ierr = DAVecGetArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
    // use Simpson's rule to do integral; therefore number of evaluations of integrand always odd; 
    // at least 13 evals per year; at least 3 in any case:
    PetscInt              pddsumcount = 12 * ((int) ceil(dt / secpera)) + 1;
    if (pddsumcount < 3)  pddsumcount = 3;
    const PetscScalar     pddsumstep = dt / pddsumcount;     // seconds
    const PetscScalar     pddsumstepdays = (365.24 / secpera) * pddsumstep; // days
    const PetscScalar     start = grid.p->year - dt / secpera;  // note grid.p->year has *end* of step
    const PetscScalar     startday = 365.24 * (start - floor(start));
    
//    ierr = verbPrintf(1,grid.com,
//            "\npddsumcount = %d, pddsumstep = %6.0f s = %6.5f a = %6.1f days, startday = %6.1f\n",
//            pddsumcount,pddsumstep,pddsumstep/secpera,pddsumstepdays,startday); CHKERRQ(ierr);

    for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
      for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {

        const PetscScalar mean_annual = Ts[i][j] - ice.meltingTemp;  // in deg C
        const PetscScalar summer_warming = getSummerWarming(h[i][j],lat[i][j],mean_annual);
        PetscScalar pdd_sum = 0.0;  // units of day^-1 (deg C)^-1
        for (PetscInt m = 0; m < pddsumcount; ++m) {
          PetscScalar  coeff = ((m % 2) == 1) ? 4.0 : 2.0;
          if (m == 0)  coeff = 1.0;   if (m == (pddsumcount - 1))  coeff = 1.0;
          const PetscScalar day = startday + m * pddsumstepdays;
          const PetscScalar temp = getTemperatureFromYearlyCycle(summer_warming,mean_annual,day);
          pdd_sum += coeff * CalovGreveIntegrand(temp);
        }
        pdd_sum = pdd_sum * (pddsumstepdays / 3.0);

        const PetscScalar snow = snow_accum[i][j] * dt;  // units of m of ice-equivalent
        const PetscScalar snow_melted = pdd_sum * pddFactorSnow;  // m of ice-equivalent
        if (snow_melted <= snow) {  // note this case is never active if snow_accum[i][j] < 0; e.g. rain
          accum[i][j] = ((snow - snow_melted) + (snow_melted * pddRefreezeFrac)) / dt;
        } else { // in this case more snow melted than fell, so some additional (old) ice also melted
          const PetscScalar excess_pdd_sum = pdd_sum - (snow / pddFactorSnow);
          const PetscScalar ice_melted = excess_pdd_sum * pddFactorIce;  // m of ice-equivalent
          accum[i][j] = - ice_melted / dt;
        }

/*
        if ((i == id) && (j == jd)) {
          verbPrintf(1,grid.com,
                "at id,jd:  h = %8.2f, Ts = %8.2f, pdd_sum = %8.2f, snow_accum = %8.4f,\n"
                "           accum = %8.4f, snow = %8.4f, snow_melted = %8.4f\n",
                h[id][jd],Ts[id][jd],pdd_sum,snow_accum[id][jd]*secpera,accum[id][jd]*secpera,
                snow, snow_melted);
        }
*/

      }
    }
    ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vAccumSnow, &snow_accum); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
    ierr = DAVecRestoreArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);

  }

  return 0;
}

