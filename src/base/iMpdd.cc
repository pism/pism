// Copyright (C) 2007 Ed Bueler and Nathan Shemonski
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
#include <gsl/gsl_sf.h>       // for erfc() in CalovGreveIntegrand()
#include <petscda.h>
#include "iceModel.hh"

//! Initialize the positive degree day (PDD) model, either deterministic (the default) or actually stochastic.
PetscErrorCode IceModel::initPDDFromOptions() {
  PetscErrorCode  ierr;

  PetscTruth      pddSet, pddRandSet, pddRepeatableSet;
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd", &pddSet); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_rand", &pddRandSet); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_rand_repeatable", &pddRepeatableSet); CHKERRQ(ierr);

  // if any of the above options are turned on, then a PDD model will be used
  if ( (doPDD == PETSC_TRUE) || (pddSet == PETSC_TRUE) || (doPDDTrueRand == PETSC_TRUE) ) {
    doPDD = PETSC_TRUE;
  } else {
    doPDD = PETSC_FALSE;
  }

  doPDDTrueRand = PETSC_FALSE;

  if (doPDD == PETSC_TRUE) {

    if (pddStuffCreated == PETSC_FALSE) {
      ierr = VecDuplicate(vh, &vAccumSnow); CHKERRQ(ierr);
      ierr = VecCopy(vAccum, vAccumSnow); CHKERRQ(ierr);  // assume vAccum *is* snow
      pddStuffCreated = PETSC_TRUE;
    } else {
      SETERRQ(1, "ERROR: initPDDFromOptions() called with pddStuffCreated == TRUE\n");
    }

    if ( (pddRandSet == PETSC_TRUE) || (pddRepeatableSet == PETSC_TRUE) ) {
      doPDDTrueRand = PETSC_TRUE;
      if (pddRandStuffCreated == PETSC_TRUE) {
        SETERRQ(2, "ERROR: initPDDFromOptions() called with -pdd_rand and pddRandStuffCreated == TRUE\n");
      }
      // initialize the random number generator: use GSL's recommended default random
      // number generator, which seems to be "mt19937" and is DIEHARD
      pddRandGen = gsl_rng_alloc(gsl_rng_default);
      gsl_rng_set(pddRandGen, pddRepeatableSet ? 0 : time(0));  // number of seconds since 1/1/1970 (?)
      pddRandStuffCreated = PETSC_TRUE;
    } else {
      pddRandStuffCreated = PETSC_FALSE;
    }

    // see iMdefaults.cc for values of DEFAULT_PDD_...
    pddFactorSnow = DEFAULT_PDD_FACTOR_SNOW;
    pddFactorIce = DEFAULT_PDD_FACTOR_ICE;
    pddRefreezeFrac = DEFAULT_PDD_REFREEZE_FRAC;
    pddSummerPeakDay = DEFAULT_PDD_SUMMER_PEAK_DAY;
    pddSummerWarming = DEFAULT_PDD_SUMMER_WARMING;
    pddStdDev = DEFAULT_PDD_STD_DEV;

    PetscTruth   pSet;
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_factor_snow", &pddFactorSnow, &pSet); CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_factor_ice", &pddFactorIce, &pSet); CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_refreeze", &pddRefreezeFrac, &pSet); CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_std_dev", &pddStdDev, &pSet); CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_summer_peak_day", &pddSummerPeakDay, &pSet); CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_summer_warming", &pddSummerWarming, &pSet); CHKERRQ(ierr);
  }
  
  return 0;
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
  const PetscScalar  rad_per_day    = 2.0 * PETSC_PI / 365.24;
  return summer_warming * cos(rad_per_day * (day - pddSummerPeakDay)) + Ta;
}


double IceModel::getSurfaceBalanceFromSnowAndPDD(
                     const double snowrate, const double mydt, const double pdds) {

  if (snowrate < 0.0) {
    return snowrate - (pdds * pddFactorIce / mydt);  // snowrate interpreted as ablation
  } else {
    const PetscScalar snow = snowrate * mydt;  // units of m of ice-equivalent (in time mydt)
    const PetscScalar snow_melted = pdds * pddFactorSnow;  // m of ice-equivalent
    if (snow_melted <= snow) {  // this case is never active if snowrate < 0
      return ((snow - snow_melted) + (snow_melted * pddRefreezeFrac)) / mydt;
    } else { // it is snowing, but all the snow melts and refreezes; this ice is
             // then removed, plus possibly more of the underlying ice
      const PetscScalar ice_created = snow * pddRefreezeFrac;
      const PetscScalar excess_pdds = pdds - (snow / pddFactorSnow); // definitely positive
      const PetscScalar ice_melted = excess_pdds * pddFactorIce;
      return (ice_created - ice_melted) / mydt;
    }
  }
}


double IceModel::CalovGreveIntegrand(const double T) {
  double z;
  const double sigSqr = pddStdDev * pddStdDev;
  z = (pddStdDev / sqrt(2.0 * PETSC_PI)) * exp(- T * T / (2.0 * sigSqr)) 
      + (T / 2.0) * gsl_sf_erfc(- T / (sqrt(2.0) * pddStdDev));
  return z;
}


//! Compute the surface mass balance from the PDD model and the stored map of snow rate. 
/*!
PISM implements two positive degree day models for computing surface mass balance from a stored 
map of snow fall rate.  There is a deterministic default method and an actually random (stochastic) method.  
Both methods conceptually include daily variability according to a normally distributed random temperature 
change whose standard deviation can be controlled by option <tt>-pdd_std_dev</tt>.  The deterministic method
computes only the expected amount of melting, while the stochastic method uses pseudo-random numbers to
simulate the melting.

Both models implement the basic positive degree day scheme described in C. Ritz, 
"EISMINT Intercomparison Experiment: Comparison of existing Greenland models,"
<tt>http://homepages.vub.ac.be/~phuybrec/eismint/greenland.html</tt>, 1997.

PISM assumes a yearly temperature cycle with average surface temperature determined by the
2D Vec <tt>vTs</tt>.  The default has a constant choice for summer warming, the difference
between the peak summer temperature and the mean of the yearly cycle.  That is, the daily mean 
temperature is
   \f[T(t) = Ts[i][j] + A \cos(2\pi (t - B) / 365.24)\f]
where \f$t\f$ gives the Julian day, \f$A\f$ is the number of degrees of summer warming and \f$B\f$ is the
Julian day of the peak summer temperature.

The amount of summer warming can also be a function of elevation, latitude, and mean annual surface 
temperature.  The derived class <tt>IceGRNModel</tt> includes an example.

The date of peak summer warming can be controlled by option <tt>-pdd_summer_peak</tt> with a given Julian day.
The magnitude of the yearly cycle by <tt>-pdd_summer_warming</tt> with a temperature change in degrees C.

PISM uses a constant rate of melting for snow, per positive degree day.  It is controlled by option
<tt>-pdd_factor_snow</tt>.  A fraction of the melted snow refreezes; the fraction is controlled by
<tt>-pdd_refreeze</tt>.  If the rate of snowfall is negative then the rate is interpreted as an ice-equivalent 
(direct) ablation rate.  If the number of positive degree days exceeds those needed to melt all of the
snow then the excess number are used to melt both the ice that came from refreeze and perhaps ice which is
already present.  In either case, ice melts at a constant rate per positive degree day, and this rate can be
controlled by option <tt>-pdd_factor_ice</tt>.

The default model only computes the expected number of positive degree days, so it is deterministic.  
It is chosen by option <tt>-pdd</tt>.  For more detail, see R. Calov and R. Greve (2005),
"Correspondence: A semi-analytical solution for the positive degree-day model with stochastic 
temperature variations," J. Glaciol. 51 (172), pp 173--175.

The alternative method, chosen by either <tt>-pdd_rand</tt> or <tt>-pdd_rand_repeatable</tt>, computes
a pseudorandom amount of temperature change for each day at each grid point.  The random amount of daily
variablility is independent and normally distributed.  (Clearly a more realistic pattern for the 
variability would have correlation with some spatial range and some temporal range.)
 */
PetscErrorCode IceModel::updateNetAccumFromPDD() {
  PetscErrorCode ierr;

  PetscScalar **accum, **Ts, **snow_accum, **h, **lat;
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vAccumSnow, &snow_accum); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);

  const PetscScalar     start = grid.p->year - dt / secpera;  // note grid.p->year has *end* of step
  const PetscScalar     startday = 365.24 * (start - floor(start));

  // set up for Calov-Greve method; use Simpson's rule to do integral; therefore number
  // of evaluations of integrand always odd; at least 13 evals per year; at least 3 in any case:
  PetscInt              pddsumcount = 12 * ((int) ceil(dt / secpera)) + 1;
  if (pddsumcount < 3)  pddsumcount = 3;
  const PetscScalar     pddsumstep = dt / pddsumcount;     // seconds
  const PetscScalar     pddsumstepdays = (365.24 / secpera) * pddsumstep; // days

  // set up for actually stochastic model
  const PetscInt        intstartday = (int) ceil(startday);
  const PetscInt        num_days = (int) ceil(365.24 * (dt / secpera));

  // run through (processor owned) grid and compute PDDs and then surface balance at each point
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar mean_annual = Ts[i][j] - ice.meltingTemp;  // in deg C
      const PetscScalar summer_warming = getSummerWarming(h[i][j],lat[i][j],mean_annual);

      PetscScalar pdd_sum = 0.0;  // units of day^-1 (deg C)-1
      
      // use one of the two methods for computing the number of positive degree days
      // at the given i,j grid point for the duration of time step (=dt)
      if (doPDDTrueRand == PETSC_TRUE) { // actually stochastic
        for (PetscInt day = intstartday; day < num_days; day++){ // compute # of pos deg day
          PetscScalar temp = getTemperatureFromYearlyCycle(summer_warming, mean_annual, 
                                                           (PetscScalar) day);
          temp += (PetscScalar) gsl_ran_gaussian(pddRandGen, pddStdDev);
          if (temp > 0.0)   pdd_sum += temp;
        }
      } else { // default Calov-Greve method; apply Simpson's rule for integral
        for (PetscInt m = 0; m < pddsumcount; ++m) {
          PetscScalar  coeff = ((m % 2) == 1) ? 4.0 : 2.0;
          if (m == 0)  coeff = 1.0;   if (m == (pddsumcount - 1))  coeff = 1.0;
          const PetscScalar day = startday + m * pddsumstepdays;
          const PetscScalar temp = getTemperatureFromYearlyCycle(summer_warming,mean_annual,day);
          pdd_sum += coeff * CalovGreveIntegrand(temp);
        }
        pdd_sum = pdd_sum * (pddsumstepdays / 3.0);
      }

      // now that we have number of PDDs, compute mass balance from snow rate
      accum[i][j] = getSurfaceBalanceFromSnowAndPDD(snow_accum[i][j], dt, pdd_sum);

      if ((i == id) && (j == jd)) {
        ierr = verbPrintf(4,grid.com,
              "\nat id,jd:  h = %6.2f, Ts = %5.2f, pdd_sum = %5.2f, snow_accum = %5.4f,\n"
              "           accum = %5.4f, num_days = %d, intstartday = %d\n",
              h[id][jd],Ts[id][jd],pdd_sum,snow_accum[id][jd]*secpera,accum[id][jd]*secpera,
              num_days, intstartday); CHKERRQ(ierr);
      }

    }
  }

  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vAccum, &accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vAccumSnow, &snow_accum); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
  return 0;
}

