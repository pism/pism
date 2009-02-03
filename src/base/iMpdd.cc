// Copyright (C) 2007-2009 Ed Bueler, Nathan Shemonski and Constantine Khroulev
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

//! Initialize the positive degree day (PDD) model
/*!
The default PDD is deterministic, but there is also a monte carlo (random)
implementation.

This procedure reads options
  - "-pdd"
  - "-pdd_rand"
  - "-pdd_rand_repeatable"
  - "-pdd_factor_ice"
  - "-pdd_factor_snow"
  - "-pdd_refreeze"
  - "-pdd_std_dev"
  - "-pdd_summer_peak_day"
  - "-pdd_summer_warming"
  - "-pdd_monthly_temps"
 */
PetscErrorCode IceModel::initPDDFromOptions() {
  PetscErrorCode  ierr;

  PetscTruth      pddSet, pddRandSet, pddRepeatableSet;
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd", &pddSet); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_rand", &pddRandSet); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_rand_repeatable", &pddRepeatableSet);
           CHKERRQ(ierr);

  // if any of the above options are turned on, then a PDD model will be used
  if ( (doPDD == PETSC_TRUE) || (pddSet == PETSC_TRUE) || (doPDDTrueRand == PETSC_TRUE) ) {
    doPDD = PETSC_TRUE;
  } else {
    doPDD = PETSC_FALSE;
  }

  doPDDTrueRand = PETSC_FALSE;

  if (doPDD == PETSC_TRUE) {

    if (pddStuffCreated == PETSC_FALSE) {
      ierr = vAccumSnow.create(grid, "vAccumSnow", true); CHKERRQ(ierr);
      ierr = vAccumSnow.copy_from(vAccum); CHKERRQ(ierr);  // assume vAccum *is* snow
      pddStuffCreated = PETSC_TRUE;
    } else {
      SETERRQ(1, "ERROR: initPDDFromOptions() called with pddStuffCreated == TRUE\n");
    }

    if ( (pddRandSet == PETSC_TRUE) || (pddRepeatableSet == PETSC_TRUE) ) {
      doPDDTrueRand = PETSC_TRUE;
      if (pddRandStuffCreated == PETSC_TRUE) {
        SETERRQ(2, "ERROR: initPDDFromOptions() called with -pdd_rand\n"
                   " and pddRandStuffCreated == TRUE\n");
      }
      // initialize the random number generator: use GSL's recommended default random
      // number generator, which seems to be "mt19937" and is DIEHARD
      pddRandGen = gsl_rng_alloc(gsl_rng_default);
      // seed with number of seconds since 1/1/1970 (?)
      gsl_rng_set(pddRandGen, pddRepeatableSet ? 0 : time(0));  
      pddRandStuffCreated = PETSC_TRUE;
    } else {
      pddRandStuffCreated = PETSC_FALSE;
    }

    PetscTruth   pSet;
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_factor_snow", &pddFactorSnow, &pSet);
             CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_factor_ice", &pddFactorIce, &pSet);
             CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_refreeze", &pddRefreezeFrac, &pSet);
             CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_std_dev", &pddStdDev, &pSet);
             CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_summer_peak_day", 
             &pddSummerPeakDay, &pSet); CHKERRQ(ierr);
    ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_summer_warming", 
             &pddSummerWarming, &pSet); CHKERRQ(ierr);

    ierr = PetscOptionsGetString(PETSC_NULL, "-pdd_monthly_temps", 
             monthlyTempsFile, PETSC_MAX_PATH_LEN, &pddHaveMonthlyTempData); CHKERRQ(ierr);
    if (pddHaveMonthlyTempData == PETSC_TRUE) {
      ierr = readMonthlyTempDataPDD(); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com,
                "using PDD based on standard yearly surface temp cycle\n");
                CHKERRQ(ierr);
    }
  }
  return 0;
}


//! De-allocate various stuff associated to the PDD model.
PetscErrorCode IceModel::PDDCleanup() {
  PetscErrorCode ierr;
  if (pddStuffCreated == PETSC_TRUE) {
    ierr = vAccumSnow.destroy(); CHKERRQ(ierr);
    pddStuffCreated = PETSC_FALSE;
  }  
  if (pddRandStuffCreated == PETSC_TRUE) {
    gsl_rng_free(pddRandGen);
    pddRandStuffCreated = PETSC_FALSE;
  }
  if (pddHaveMonthlyTempData == PETSC_TRUE) {
    if (vmonthlyTs != PETSC_NULL) {
      for (PetscInt j = 0; j < 12; ++j) {
        if (vmonthlyTs[j] != PETSC_NULL) {
          ierr = VecDestroy(vmonthlyTs[j]); CHKERRQ(ierr);
        }
      }
      vmonthlyTs = PETSC_NULL;
    }
    pddHaveMonthlyTempData = PETSC_FALSE;
  }
  return 0;
}


//! Read in monthly temperature data if available (and allocate space for it).
/*!
Called if option "-pdd_monthly_temps" is given with a file name.
 */
PetscErrorCode IceModel::readMonthlyTempDataPDD() {
  PetscErrorCode ierr;
  NCTool nc(&grid);

  ierr = verbPrintf(2, grid.com, 
       "reading monthly surface temperature data for PDD from file %s ...\n", 
       monthlyTempsFile); CHKERRQ(ierr);

  bool file_exists = false;
  ierr = nc.open_for_reading(monthlyTempsFile, file_exists); CHKERRQ(ierr);
  if (!file_exists) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: Can't open file '%s'.\n", monthlyTempsFile);
    CHKERRQ(ierr);
    PetscEnd();
  }

  ierr = verbPrintf(4, grid.com, 
         "  creating local interpolation context for monthly surface temps ...\n");
  CHKERRQ(ierr);
  grid_info g;
  ierr = nc.get_grid_info_2d(g); CHKERRQ(ierr);
  LocalInterpCtx lic(g, NULL, NULL, grid); // 2D only
  //ierr = lic.printGrid(grid.com); CHKERRQ(ierr);
  
  ierr = verbPrintf(4, grid.com, 
         "  allocating space for monthly surface temps ...\n"); CHKERRQ(ierr);
  ierr = VecDuplicateVecs(g2, 12, &vmonthlyTs); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, 
         "  reading (regridding) monthly surface temps ('temp_monN'):\n    month  ");
  CHKERRQ(ierr);

  bool var_exists;
  int varid;
  for (PetscInt j = 0; j < 12; ++j) {
    char monthlyTempName[20];
    snprintf(monthlyTempName, 20, "temp_mon%d", j);
    ierr = nc.find_variable(monthlyTempName, NULL, &varid, var_exists);
    if (!var_exists) {
      ierr = PetscPrintf(grid.com,
			 "PISM ERROR: Variable '%s' not found in file '%s'\n.",
			 monthlyTempName, monthlyTempsFile);
      CHKERRQ(ierr);
      PetscEnd();
    }
    ierr = nc.regrid_global_var(varid, GRID_2D, lic, grid.da2, 
				vmonthlyTs[j], false); CHKERRQ(ierr);
    ierr = verbPrintf(2, grid.com, " %d", j); CHKERRQ(ierr);
  }
  ierr = verbPrintf(2, grid.com, "\n"); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);
  return 0;
}


//! Put the stored map vAccumSnow back in vAccum.
/*!
We need to make sure that when vAccum is written to NetCDF at end, it is snow accumulation
rate and not the surface mass balance (which is vAccum).
 */
PetscErrorCode IceModel::putBackSnowAccumPDD() {
  PetscErrorCode ierr;
  ierr = vAccumSnow.copy_to(vAccum); CHKERRQ(ierr);  
  return 0;
}


//!  Return the amount of summer warming.
/*!
The amount of summer warming can also be a function of elevation, latitude, and mean 
annual surface temperature.  The derived class IceGRNModel includes an example.

The date of peak summer warming is controlled by option <tt>-pdd_summer_peak_day</tt> 
with a given Julian day.  The magnitude of the yearly cycle is controlled by 
<tt>-pdd_summer_warming</tt> with a temperature change in degrees C.

This procedure is \e virtual and replacable.
 */
PetscScalar IceModel::getSummerWarming(
       const PetscScalar elevation, const PetscScalar latitude, const PetscScalar Tma) const {
  // version here ignors elevation, latitude, and mean annual temperature (Tma)
  // and instead uses -pdd_summer_warming setable constant
  // see IceGRNModel for alternate implementation
  return pddSummerWarming;
}


//! Get the surface temperature at a point for a given time (in days) from yearly cycle.
/*!
There are two versions:

  - if pddHaveMonthlyTempData is true then the monthly data read by 
    readMonthlyTempDataPDD() is linearly interpolated, while

  - if pddHaveMonthlyTempData is false then a standard sinusoidal formula is used:
      \f[ T_s(d,i,j) = T_{ma} + S \cos((2\pi/365.24) * (d - P)) \f]
    where \f$d\f$ is the day on a 365.24 day calendar,
    \f$T_{ma}\f$ is the annual mean temperature in degrees C,
    \f$S\f$ is the amplitude of the sinusoid (the ``summer warming'') in degrees C,
    and \f$P\f$ is the day of peak summer warming, pddSummerPeakDay, which is rather
    arbitrarily taken to be August 1st by default.

The first version runs when option "-pdd_monthly_temps" is given with a file name.

EISMINT-Greenland assumes the second version \lo\cite{RitzEISMINT}\elo.  
See <tt>iceGRNModel.cc</tt>.

This procedure is \e virtual and replacable.
 */
PetscScalar IceModel::getTemperatureFromYearlyCycle(
       const PetscScalar summer_warming, const PetscScalar Tma, const PetscScalar day,
       const PetscInt i, const PetscInt j) const {
  
  if (pddHaveMonthlyTempData == PETSC_TRUE) {
    // this case is when "-have_monthly_temps foo.nc"
    PetscScalar month = 12.0 * day / 365.24;
    month = month - static_cast<PetscScalar> (((int) floor(month)) % 12);
    PetscInt  curr = (int) floor(month);
    curr = curr % 12;
    PetscInt  next = curr+1;
    if (next == 12)   next = 0;
    PetscScalar **currTs, **nextTs;
    DAVecGetArray(grid.da2, vmonthlyTs[curr], &currTs);
    DAVecGetArray(grid.da2, vmonthlyTs[next], &nextTs);
    const PetscScalar myTs = (currTs[i][j] - 273.15) 
         + (month - (PetscScalar)curr) * (nextTs[i][j] - 273.15);
    DAVecRestoreArray(grid.da2, vmonthlyTs[curr], &currTs);
    DAVecRestoreArray(grid.da2, vmonthlyTs[next], &nextTs);
    return myTs;
  } else { 
    // this is the default case for EISMINT-Greenland
    const PetscScalar  rad_per_day    = 2.0 * PETSC_PI / 365.24;
    return Tma + summer_warming * cos(rad_per_day * (day - pddSummerPeakDay));
  }
}


//! Compute the surface balance at a location given the number of positive degree days and snowfall there.
/*!  
The surface balance (either net accumulation or ablation) is computed from the
number of positive degree days and the yearly snowfall.  The number of positive degree
days is determined inside updateSurfaceBalanceFromPDD().

We assume a constant rate of melting per positive degree day for snow.  The rate is set
by the option <tt>-pdd_factor_snow</tt>.  A fraction of the melted snow refreezes; 
this fraction is controlled by <tt>-pdd_refreeze</tt>.  If the number of positive 
degree days exceeds those needed to melt all of the snow then the excess number are 
used to melt both the ice that came from refreeze and perhaps ice which is already 
present.  In either case, ice melts at a constant rate per positive degree day, 
and this rate can be controlled by option <tt>-pdd_factor_ice</tt>.

If the rate of snowfall is negative then the rate is interpreted as an ice-equivalent
(direct) ablation rate and the PDD contribution is added (i.e. there is additional 
ablation).  For most uses the snowfall rate is positive everywhere on the ice sheet, 
however.

The default values for the factors come from \lo\cite{RitzEISMINT}\elo.
 */
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


//! Compute the integrand in integral (6) in \lo\cite{CalovGreve05}\elo.
/*!
The integral is
   \f[\mathrm{PDD} = \int_{t_0}^{t_0+\mathtt{dt}} dt\,
         \bigg[\frac{\sigma}{\sqrt{2\pi}}\,\exp\left(-\frac{T_{ac}(t)^2}{2\sigma^2}\right)
               + \frac{T_{ac}(t)}{2}\,\mathrm{erfc}
               \left(-\frac{T_{ac}(t)}{\sqrt{2}\,\sigma}\right)\bigg] \f]
This procedure computes the quantity in square brackets.

The user can choose \f$\sigma\f$ by option <tt>-pdd_std_dev</tt>.

Note that the integral is over a time interval of length \c dt instead of a 
whole year as stated in \lo\cite{CalovGreve05}\elo.
 */
double IceModel::CalovGreveIntegrand(const double Tac) {

  return (pddStdDev / sqrt(2.0 * PETSC_PI))
         * exp(- Tac * Tac / (2.0 * pddStdDev * pddStdDev)) 
         + (Tac / 2.0) * gsl_sf_erfc(- Tac / (sqrt(2.0) * pddStdDev));
}


//! Compute the surface mass balance over the whole grid from the PDD model and a stored map of snow rate. 
/*!
PISM implements two positive degree day models for computing surface mass balance 
from a stored map of snow fall rate.  There is a deterministic default method and an 
monte carlo (stochastic) method.  Both methods conceptually include daily and 
additional weather-related variability according to a normally distributed 
random temperature change for each day and grid point.  The standard deviation of
this temperature change can be controlled by option <tt>-pdd_std_dev</tt>.  
The deterministic method computes only the expected amount of melting, while the 
monte carlo method uses pseudo-random numbers to simulate the melting.

computes the yearly temperature cycle.

The surface mass balance at each point of the grid is computed by calling
getSurfaceBalanceFromSnowAndPDD().

The default model only computes the \e expected number of positive degree days, 
so it is deterministic.  It is chosen by option <tt>-pdd</tt>.  It implements the 
scheme in \lo\cite{CalovGreve05}\elo.  In particular, integral (6) in that paper
is approximated here by Simpson's rule.

The alternative method, chosen by either <tt>-pdd_rand</tt> or 
<tt>-pdd_rand_repeatable</tt>, computes a (pseudo)random amount of temperature 
change for each day at each grid point.  This change is normally distributed 
and is independent at each grid point and each day.  This method can be regarded 
as a monte carlo simulation of a stochastic process (of which the Calov-Greve 
method computes the expected value.

A more realistic pattern for the variability of surface melting would have correlation 
with appropriate spatial and temporal ranges.
 */
PetscErrorCode IceModel::updateSurfaceBalanceFromPDD() {
  PetscErrorCode ierr;

  PetscScalar **accum, **Ts, **snow_accum, **h, **lat;
  ierr =         vh.get_array(h);          CHKERRQ(ierr);
  ierr =        vTs.get_array(Ts);         CHKERRQ(ierr);
  ierr =     vAccum.get_array(accum);      CHKERRQ(ierr);
  ierr = vAccumSnow.get_array(snow_accum); CHKERRQ(ierr);
  ierr =  vLatitude.get_array(lat);        CHKERRQ(ierr);

  const PetscScalar     start = grid.year - dt / secpera; // grid.year has *end* of step
  const PetscScalar     startday = 365.24 * (start - floor(start));

  // set up for Calov-Greve method; use Simpson's rule to do integral, so number
  // of evaluations of integrand always odd; at least 53 evals per year (i.e. approximately
  // weekly); see trials at end of this file;  at least 3 evals in any case:
  PetscInt              CGsumcount = (int) ceil(52 * (dt / secpera) + 1);
//  PetscInt              CGsumcount = (int) ceil(366 * (dt / secpera) + 1);
  if (CGsumcount < 3)   CGsumcount = 3;
  if ((CGsumcount % 2) == 0)  CGsumcount++;  // guarantee it is odd
  const PetscScalar     CGsumstep = dt / CGsumcount;     // seconds
  const PetscScalar     CGsumstepdays = (365.24 / secpera) * CGsumstep; // days

  // set up for monte carlo method
  const PetscInt        intstartday = (int) ceil(startday);
  const PetscInt        num_days = (int) ceil(365.24 * (dt / secpera));

  PetscScalar pdd_sum;  // units of day^-1 (deg C)-1

  // run through grid and compute PDDs and then surface balance at each point
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      const PetscScalar mean_annual = Ts[i][j] - ice->meltingTemp;  // in deg C
      const PetscScalar summer_warming = getSummerWarming(h[i][j],lat[i][j],mean_annual);

      // use one of the two methods for computing the number of positive degree days
      // at the given i,j grid point for the duration of time step (=dt)
      if (doPDDTrueRand == PETSC_TRUE) { // monte carlo
        pdd_sum = 0.0;
        // compute # of pos deg day:
        for (PetscInt day = intstartday; day < intstartday + num_days; day++){ 
          const PetscScalar mytemp
                      = getTemperatureFromYearlyCycle(
                             summer_warming, mean_annual,(PetscScalar) day, i,j);
          const double randadd = gsl_ran_gaussian(pddRandGen, pddStdDev);
          const PetscScalar temp = mytemp + (PetscScalar) randadd;
          if (temp > 0.0)   pdd_sum += temp;
          if ((i == id) && (j == jd)) {
            ierr = verbPrintf(5,grid.com,
              "  day=%d: mytemp=%5.4f, randadd=%5.4f, temp=%5.4f, pdd_sum=%5.4f\n",
              day, mytemp, randadd, temp, pdd_sum); CHKERRQ(ierr);
          }
        }
      } else { // default Calov-Greve method; apply Simpson's rule for integral
        pdd_sum = 0.0;
        for (PetscInt m = 0; m < CGsumcount; ++m) {
          // Simpson's is: (h/3) * sum([1 4 2 4 2 4 ... 4 1] .* [f(x0) f(x1) ... f(xN)])
          PetscScalar  coeff = ((m % 2) == 1) ? 4.0 : 2.0;
          if ( (m == 0) || (m == (CGsumcount - 1)) )  coeff = 1.0;
          const PetscScalar day = startday + m * CGsumstepdays;
          const PetscScalar temp = getTemperatureFromYearlyCycle(
                                       summer_warming,mean_annual,day, i,j);
          pdd_sum += coeff * CalovGreveIntegrand(temp);
        }
        pdd_sum = (CGsumstepdays / 3.0) * pdd_sum;
      }

      // now that we have the number of PDDs, compute mass balance from snow rate
      accum[i][j] = getSurfaceBalanceFromSnowAndPDD(snow_accum[i][j], dt, pdd_sum);

      if ((i == id) && (j == jd)) {
        ierr = verbPrintf(5,grid.com,
          "\nPDD at (i,j)=(id,jd)=(%d,%d):\n"
            "  mean_annual=%5.4f, summer_warming=%5.4f\n",
            i, j, mean_annual, summer_warming); CHKERRQ(ierr);
        ierr = verbPrintf(5,grid.com,
            "  CGsumcount = %d, CGsumstepdays = %5.2f, num_days = %d, intstartday = %d,\n"
            "  h = %6.2f, pdd_sum = %5.2f, snow_accum = %5.4f, accum = %5.4f\n",
            CGsumcount, CGsumstepdays, num_days, intstartday,
            h[i][j], pdd_sum, snow_accum[i][j]*secpera, accum[i][j]*secpera); CHKERRQ(ierr);
      }

    }
  }

  ierr = vTs.end_access(); CHKERRQ(ierr);
  ierr = vAccum.end_access(); CHKERRQ(ierr);
  ierr = vAccumSnow.end_access(); CHKERRQ(ierr);
  ierr = vh.end_access(); CHKERRQ(ierr);
  ierr = vLatitude.end_access(); CHKERRQ(ierr);
  return 0;
}


/*
Following are trials for the Calov-Greve versus monte carlo methods.  The 
conclusion is that Simpson's rule with at least 53 points per year---
weekly steps in the quadrature---probably suffices.  But the monthly steps
in \lo\cite{CalovGreve05}\elo are probably not sufficient.

DATA GENERATED BY CALLS
   $ mpiexec -n 2 pgrn -boot_from eis_green_smoothed.nc -Mx 83 -My 141 -Mz 201 \
       -ocean_kill -y 50 -o green_50yr_pdd
AND
   $ mpiexec -n 2 pgrn -boot_from eis_green_smoothed.nc -Mx 83 -My 141 -Mz 201 \
       -ocean_kill -y 50 -pdd_rand -o green_50yr_pdd_randN

initial:
S      0.00000:  2.82500  1.6708   0.2071   3042.000  270.5156

C-G pdd with 52 Simpson's points per year:
+0.57220
S     50.00000:  2.84077  1.6620   0.2385   3019.121  270.5353

C-G pdd with 366 Simpson's points per year:
+0.57177
S     50.00000:  2.83997  1.6580   0.2388   3019.117  270.5353

pdd_rand  trial 0:
+0.57156
S     50.00000:  2.83958  1.6604   0.2385   3019.114  270.5353

pdd_rand  trial 1:
+0.57174
S     50.00000:  2.83962  1.6632   0.2381   3019.112  270.5353

pdd_rand  trial 2:
+0.57158
S     50.00000:  2.83957  1.6608   0.2384   3019.112  270.5353
*/
