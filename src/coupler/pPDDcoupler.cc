// Copyright (C) 2007--2009 Ed Bueler, Nathan Shemonski, Constantine Khroulev
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


#include <ctime>  // for time(), used to initialize random number gen
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>       // for erfc() in CalovGreveIntegrand()
#include <petscda.h>
#include "../base/pism_const.hh"
#include "../base/grid.hh"
#include "../base/iceModelVec.hh"
#include "pccoupler.hh"


PISMPDDCoupler::PISMPDDCoupler() : PISMConstAtmosCoupler() {

  // normally use Calov-Greve expectation integral; if this is non-NULL
  //   then using actual random numbers
  pddRandGen = NULL;
  
  // see pPDDcoupler.hh for units and meaning
  pddFactorSnow = 0.003;  // m K^-1 day^-1; EISMINT-Greenland value; = 3 mm / (pos degree day)
  pddFactorIce  = 0.008;  // m K^-1 day^-1; EISMINT-Greenland value; = 8 mm / (pos degree day)
  pddRefreezeFrac = 0.6;  // [pure fraction]; EISMINT-Greenland value
  pddStdDev       = 5.0;  // K; std dev of daily temp variation; EISMINT-Greenland value
  pddSummerWarming = 15.0;// K
     // re SUMMER_WARMING:  (30.38 - 0.006277 * 1000.0 - 0.3262 * 75.0)
     //                    - (49.13 - 0.007992 * 1000.0 -0.7576 * 75.0)
     //                   =  15.32   K
     // is result of EISMINT-Greenland formulas for h=1000.0 m and lat=75.0 deg N
  pddSummerPeakDay = 243.0;  // Julian day; = August 1st
  
  usingMonthlyTemps = PETSC_FALSE; // -pdd_monthly_temps option required to use them
}


PISMPDDCoupler::~PISMPDDCoupler() {
  vannmeansurftemp.destroy();
  vsurfaccum.destroy();
  if (pddRandGen != NULL) {
    gsl_rng_free(pddRandGen);
    pddRandGen = NULL;
  }
  for (PetscInt j = 0; j < 12; ++j) {
    vmonthlysurftemp[j].destroy();
  }
}


//! Check if a positive degree day (PDD) model is desired.
/*!
This procedure determines if any PDD relates options are set.
If so, returns PETSC_TRUE in userWantsPDD, otherwise PETSC_FALSE.
 */
PetscErrorCode PISMPDDCoupler::userOptionsChoosePDD(PetscTruth &userWantsPDD) {
  PetscErrorCode  ierr;
  const int nNames = 10;
  const char pddOptNames[nNames][30] = {
           "-pdd",
           "-pdd_factor_snow", "-pdd_factor_ice", "-pdd_refreeze", 
           "-pdd_rand", "-pdd_rand_repeatable", "-pdd_std_dev",
           "-pdd_summer_warming", "-pdd_summer_peak_day",
           "-pdd_monthly_temps"};

  for (int k=0; k<nNames; k++) {
    PetscTruth check;
    ierr = PetscOptionsHasName(PETSC_NULL, pddOptNames[k], &check); CHKERRQ(ierr);
    if (check == PETSC_TRUE) {
      userWantsPDD = PETSC_TRUE;
      return 0;
    }
  }
  userWantsPDD = PETSC_FALSE;
  return 0;
}


PetscErrorCode PISMPDDCoupler::initFromOptions(IceGrid* g) {
  PetscErrorCode ierr;

ierr = verbPrintf(1,g->com,"entering PISMPDDCoupler::initFromOptions()\n"); CHKERRQ(ierr);

  // check options; we assume the PDD is desired; only overriding defaults here
  // check if truly random PDD is desired, and, if so, initialize it
  PetscTruth  pddRandSet, pddRepeatableSet;
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_rand", &pddRandSet); CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL, "-pdd_rand_repeatable", &pddRepeatableSet);
           CHKERRQ(ierr);
  if ( (pddRandSet == PETSC_TRUE) || (pddRepeatableSet == PETSC_TRUE) ) {
    // initialize the random number generator: use GSL's recommended default random
    // number generator, which seems to be "mt19937" and is DIEHARD
    pddRandGen = gsl_rng_alloc(gsl_rng_default);  // so pddRandGen != NULL now
    // seed with number of seconds since 1/1/1970 (?)
    gsl_rng_set(pddRandGen, pddRepeatableSet ? 0 : time(0));  
  }

  // check options for parameter values
  PetscTruth   pSet;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_factor_snow", &pddFactorSnow, &pSet);
             CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_factor_ice", &pddFactorIce, &pSet);
             CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_refreeze", &pddRefreezeFrac, &pSet);
             CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_std_dev", &pddStdDev, &pSet);
             CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_summer_warming", 
             &pddSummerWarming, &pSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_summer_peak_day", 
             &pddSummerPeakDay, &pSet); CHKERRQ(ierr);

  // read accumulation and surface temperature
  char  filename[PETSC_MAX_PATH_LEN];
  LocalInterpCtx* lic;

  ierr = PISMConstAtmosCoupler::initFromOptions(g); CHKERRQ(ierr); // sets grid and metadata
  ierr = findPISMInputFile((char*)filename, lic); CHKERRQ(ierr); // allocates lic

  // mean annual ice equivalent accumulation rate
  //   we READ by name 'artm' but change to name 'annavartm' below
  ierr = vannmeansurftemp.create(*g, "artm", true); CHKERRQ(ierr);
  ierr = vannmeansurftemp.set_attrs(
            "climate_state",
            "annual mean temperature at ice surface but below firn",
            "K", NULL); CHKERRQ(ierr);
  ierr = vannmeansurftemp.set(273.15); CHKERRQ(ierr);  // merely a default value

  // mean annual ice equivalent accumulation rate
  //   we READ by name 'acab' but change name to 'accum' below
  ierr = vsurfaccum.create(*g, "acab", false); CHKERRQ(ierr);  // global; no ghosts
  ierr = vsurfaccum.set_attrs(
            "climate_state", 
            "mean annual ice equivalent accumulation rate",
	    "m s-1", 
	    NULL);  // no CF standard_name
	    CHKERRQ(ierr);
  ierr = vsurfaccum.set_glaciological_units("m year-1");
  vsurfaccum.write_in_glaciological_units = true;
  ierr = vsurfaccum.set(0.0); CHKERRQ(ierr);  // merely a default value

  // now read two fields
  ierr = verbPrintf(2, g->com, 
     "initializing constant atmospheric climate and PDD: reading surface\n"
     "  accumulation 'acab' and surface temperature 'artm' from %s ... \n",
     filename); CHKERRQ(ierr); 

  ierr = vsurfaccum.regrid(filename, *lic, true); CHKERRQ(ierr); // it *is* critical
  ierr = vannmeansurftemp.regrid(filename, *lic, true); CHKERRQ(ierr); // it *is* critical

  // now that they are read, reset names to output name
  ierr = vsurfaccum.set_name("accum"); CHKERRQ(ierr);
  ierr = vannmeansurftemp.set_name("annavartm"); CHKERRQ(ierr);
  ierr = vsurfmassflux.set(0.0); CHKERRQ(ierr); // initialize to have some values

  // check on whether we should read monthly temperatures from file  
  char monthlyTempsFile[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-pdd_monthly_temps", 
             monthlyTempsFile, PETSC_MAX_PATH_LEN, &pSet); CHKERRQ(ierr);
  if (pSet == PETSC_TRUE) {
    usingMonthlyTemps = PETSC_TRUE;
    ierr = verbPrintf(2,grid->com,
       "using PDD based on monthly temperatures; reading from %s ...\n",
       monthlyTempsFile); CHKERRQ(ierr);
    ierr = readMonthlyTemps(monthlyTempsFile); CHKERRQ(ierr);
  } else {
    usingMonthlyTemps = PETSC_FALSE;
    ierr = verbPrintf(2,grid->com,
       "using PDD based on standard yearly surface temp cycle ...\n"); CHKERRQ(ierr);
  }
  
  delete lic;

ierr = verbPrintf(1,g->com,"ending PISMPDDCoupler::initFromOptions()\n"); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PISMPDDCoupler::writeCouplingFieldsToFile(const char *filename) {
  PetscErrorCode ierr;
  
  ierr = PISMConstAtmosCoupler::writeCouplingFieldsToFile(filename); CHKERRQ(ierr);
  
  ierr = vannmeansurftemp.write(filename, NC_FLOAT); CHKERRQ(ierr);

  ierr = vsurfaccum.write(filename, NC_FLOAT); CHKERRQ(ierr);

  if (usingMonthlyTemps == PETSC_TRUE) {
    for (PetscInt j = 0; j < 12; ++j) {
      if (vmonthlysurftemp[j].was_created()) {
        ierr = vmonthlysurftemp[j].write(filename, NC_FLOAT); CHKERRQ(ierr);
      }
    }
  }

  return 0;
}


PetscErrorCode PISMPDDCoupler::readMonthlyTemps(const char *filename) {
  PetscErrorCode ierr;
  
  // find file and set up info so regrid works
  bool file_exists = false;
  NCTool nc(grid);
  grid_info gi;
  ierr = nc.open_for_reading(filename, file_exists); CHKERRQ(ierr);
  if (!file_exists) {
    ierr = PetscPrintf(grid->com,
       "PISMPDDCoupler ERROR: Can't open file '%s' for reading monthly temps.\n",
       filename); CHKERRQ(ierr);
    PetscEnd();
  }
  ierr = nc.get_grid_info_2d(gi); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);
  LocalInterpCtx lic(gi, NULL, NULL, *grid); // 2D only; 

  // for each month, create an IceModelVec2 and assign attributes
  for (PetscInt j = 0; j < 12; ++j) {
    char monthlyTempName[20], mTstring[100];
    snprintf(monthlyTempName, 20, "temp_mon%d", j);
    ierr = verbPrintf(2, grid->com, 
       "  reading month %d surface temperature '%s' ...\n", j, monthlyTempName); CHKERRQ(ierr); 
    ierr = vmonthlysurftemp[j].create(*grid, monthlyTempName, false); // global; no ghosts
       CHKERRQ(ierr);
    snprintf(mTstring, 100, 
             "temperature at ice surface but below firn during month %d of {0,..,11}", j);
    ierr = vmonthlysurftemp[j].set_attrs("climate_state",mTstring,"K",NULL); CHKERRQ(ierr);
    ierr = vmonthlysurftemp[j].set(273.15); CHKERRQ(ierr);  // merely a default value
    ierr = vmonthlysurftemp[j].regrid(filename, lic, true); CHKERRQ(ierr); // it *is* critical
  }
  return 0;
}


//!  Return the amount of summer warming.
/*!
The amount of summer warming can also be a function of elevation, latitude, and mean 
annual surface temperature.

The date of peak summer warming is controlled by option <tt>-pdd_summer_peak_day</tt> 
with a given Julian day.  The magnitude of the yearly cycle is controlled by 
<tt>-pdd_summer_warming</tt> with a temperature change in degrees C.

This procedure is \e virtual and replacable.
 */
PetscScalar PISMPDDCoupler::getSummerWarming(
       const PetscScalar elevation, const PetscScalar latitude, const PetscScalar Tma) {
  // version here ignors elevation, latitude, and mean annual temperature (Tma)
  // and instead uses -pdd_summer_warming setable constant
  // see IceGRNModel for alternate implementation
  return pddSummerWarming;
}


//! Get the surface temperature at a point for a given time (in days) from yearly cycle.
/*!
A standard sinusoidal formula is used:
      \f[ T_s(d,i,j) = T_{ma} + S \cos((2\pi/365.24) * (d - P)) \f]
where \f$d\f$ is the day on a 365.24 day calendar, \f$T_{ma}\f$ is the annual
mean temperature in degrees C, \f$S\f$ is the amplitude of the sinusoid
(the ``summer warming'') in degrees C, and \f$P\f$ is the day of peak summer warming,
pddSummerPeakDay, which is taken to be 1 August by default.

This follows EISMINT-Greenland \lo\cite{RitzEISMINT}\elo.  See also IceGRNModel.
 */
PetscScalar PISMPDDCoupler::getTemperatureFromYearlyCycle(
       const PetscScalar summer_warming, const PetscScalar Tma, 
       const PetscScalar day) {
  const PetscScalar  rad_per_day = 2.0 * PETSC_PI / 365.24;
  return Tma + summer_warming * cos(rad_per_day * (day - pddSummerPeakDay));
}


/*!
Returns indices in {0,..,11} for the current and next months.  For the purpose
of indexing the monthly surface temperature data.
 */
PetscErrorCode PISMPDDCoupler::getMonthIndicesFromDay(const PetscScalar day, 
       PetscInt &curr, PetscInt &next) {
  PetscScalar month = 12.0 * day / 365.24;
  month = month - static_cast<PetscScalar> (((int) floor(month)) % 12);
  curr = (int) floor(month);
  curr = curr % 12;
  next = curr+1;
  if (next == 12)   next = 0;
  return 0;
}


/*!
Linearly interpolates between stored monthly temps.
 */
PetscScalar PISMPDDCoupler::getTemperatureFromMonthlyData(
       PetscScalar **currMonthSurfTemps, PetscScalar **nextMonthSurfTemps,
       const PetscInt i, const PetscInt j, const PetscScalar day) {
  PetscScalar month = 12.0 * day / 365.24;
  month = month - static_cast<PetscScalar> (((int) floor(month)) % 12);
  PetscInt  curr = (int) floor(month);
  curr = curr % 12;
  return (currMonthSurfTemps[i][j] - 273.15) 
       + (month - (PetscScalar)curr) * (nextMonthSurfTemps[i][j] - 273.15);
}


//! Compute the surface balance at a location given number of positive degree days and snowfall.
/*!  
The net surface mass balance, as ice equivalent thickness per time, is computed 
from the number of positive degree days and the yearly snowfall.

The number, or expected number, of positive degree days must be determined before 
calling this procedure.

We assume a constant rate of melting per positive degree day for snow.  The rate is set
by the option <tt>-pdd_factor_snow</tt>.  A fraction of the melted snow refreezes; 
this fraction is controlled by <tt>-pdd_refreeze</tt>.  If the number of positive 
degree days exceeds those needed to melt all of the snow then the excess number are 
used to melt both the ice that came from refreeze and perhaps ice which is already
present.

In either case, ice also melts at a constant rate per positive degree day, 
and this rate can be controlled by option <tt>-pdd_factor_ice</tt>.

If the rate of snowfall is negative then the rate is interpreted as an ice-equivalent
(direct) ablation rate and the PDD contribution is added (i.e. there is additional 
ablation), by melting ice.  Snowfall rates are generally positive nearly everywhere
on ice sheets, however.

The default values for the factors come from EISMINT-Greenland, \lo\cite{RitzEISMINT}\elo.

Arguments are snow fall rate snowrate in m * s^-1, dt in s, pddsum in degree (K) * day.
 */
PetscScalar PISMPDDCoupler::getSurfaceBalanceFromSnowAndPDD(
             const PetscScalar snowrate, const PetscScalar dt_secs, const PetscScalar pddsum) {

  if (snowrate < 0.0) {
    return snowrate - (pddsum * pddFactorIce / dt_secs);  // neg snowrate interpreted as ablation
  } else {
    const PetscScalar snow        = snowrate * dt_secs,  // units of m of ice-equivalent
                      snow_melted = pddsum * pddFactorSnow;  // m of ice-equivalent
    if (snow_melted <= snow) {
      return ((snow - snow_melted) + (snow_melted * pddRefreezeFrac)) / dt_secs;
    } else { // it is snowing, but all the snow melts and refreezes; this ice is
             // then removed, plus possibly more of the underlying ice
      const PetscScalar ice_created   = snow * pddRefreezeFrac,
                        excess_pddsum = pddsum - (snow / pddFactorSnow), // positive!
                        ice_melted    = excess_pddsum * pddFactorIce;
      return (ice_created - ice_melted) / dt_secs;
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

This integral is used for the expected number of positive degree days, unless the
user selects a random PDD implementation with <tt>-pdd_rand</tt> or 
<tt>-pdd_rand_repeatable</tt>.  The user can choose \f$\sigma\f$ by option
<tt>-pdd_std_dev</tt>.  Note that the integral is over a time interval of length
\c dt instead of a whole year as stated in \lo\cite{CalovGreve05}\elo.
 */
double PISMPDDCoupler::CalovGreveIntegrand(const double Tac) {

  return (pddStdDev / sqrt(2.0 * PETSC_PI))
         * exp(- Tac * Tac / (2.0 * pddStdDev * pddStdDev)) 
         + (Tac / 2.0) * gsl_sf_erfc(- Tac / (sqrt(2.0) * pddStdDev));
}


//! Compute the net surface mass balance from the PDD model and a stored map of snow rate. 
/*!
PISM implements two positive degree day models for computing surface mass balance 
from a stored map of snow fall rate.  There is a deterministic default method and an 
monte carlo (stochastic) method.  Both methods conceptually include a yearly cycle 
of temperature and additional weather-related variability according to a normally distributed 
random temperature change for each day and grid point.  The standard deviation of
this temperature change can be controlled by option <tt>-pdd_std_dev</tt>.  
The deterministic method computes only the expected amount of melting, while the 
monte carlo method uses pseudo-random numbers to simulate the melting.

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

The spatial distribution of temperature is, by default, read from vsurftemp.
This temperature is interpreted as the annual mean temperature at each location.

Alternatively, if option <tt>-pdd_monthly_temps</tt> provides a NetCDF file 
with the 12 monthly temperature maps, the temperature on a given day at a given
location is found by linear interpolation (in time) of the monthly temps.
 */
PetscErrorCode PISMPDDCoupler::updateSurfMassFluxAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years, 
             void *iceInfoNeeded, // will be interpreted as type iceInfoNeededByAtmosphereCoupler*
             IceModelVec2* &pvsmf) {

  PetscErrorCode ierr;
  IceInfoNeededByAtmosphereCoupler* info = (IceInfoNeededByAtmosphereCoupler*) iceInfoNeeded;

  PetscScalar **smflux, **amstemp, **saccum, **h, **lat, **smonthtemp[12];
  ierr = info->surfelev->get_array(h);   CHKERRQ(ierr);
  ierr = info->lat->get_array(lat); CHKERRQ(ierr);

  ierr =    vsurfmassflux.get_array(smflux);  CHKERRQ(ierr);
  ierr = vannmeansurftemp.get_array(amstemp); CHKERRQ(ierr);
  ierr =       vsurfaccum.get_array(saccum);  CHKERRQ(ierr);

  if (usingMonthlyTemps == PETSC_TRUE) {
    for (PetscInt j = 0; j < 12; ++j) {
      if (vmonthlysurftemp[j].was_created()) {
        ierr = vmonthlysurftemp[j].get_array(smonthtemp[j]); CHKERRQ(ierr);
      } else {
        SETERRQ1(1,"vmonthlysurftemp[%d] not created",j);
      }
    }
  }

  const PetscScalar startday = 365.24 * (t_years - floor(t_years));

  // set up for Calov-Greve method; use Simpson's rule to do integral, so number
  // of evaluations of integrand always odd; at least 53 evals per year (i.e. approximately
  // weekly); see trials at end of this file;  at least 3 evals in any case:
  PetscInt            CGsumcount = (int) ceil(52 * (dt_years) + 1);
  if (CGsumcount < 3) CGsumcount = 3;
  if ((CGsumcount % 2) == 0)  CGsumcount++;  // guarantee it is odd

  const PetscScalar CGsumstep     = dt_years * secpera / CGsumcount, // seconds
                    CGsumstepdays = (365.24 / secpera) * CGsumstep;  // days

  // set up for monte carlo method
  const PetscInt    intstartday = (int) ceil(startday),
                    num_days    = (int) ceil(365.24 * dt_years);

  PetscScalar pdd_sum;  // units of day^-1 (deg C)-1
       
  // run through grid and compute PDDs and then surface balance at each point
  for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
      const PetscScalar mean_annual = amstemp[i][j] - 273.15,  // in deg C
                        summer_warming = getSummerWarming(h[i][j],lat[i][j],mean_annual);

      // use one of the two methods for computing the number of positive degree days
      // at the given i,j grid point for the duration of time step (=dt)
      if (pddRandGen != NULL) { // since random stuff set up, do monte carlo
        pdd_sum = 0.0;
        // compute # of pos deg day:
        for (PetscInt day = intstartday; day < intstartday + num_days; day++){ 
          PetscScalar mytemp;
          if (usingMonthlyTemps == PETSC_TRUE) {
            PetscInt currMonthInd, nextMonthInd;
            ierr = getMonthIndicesFromDay((PetscScalar) day, currMonthInd, nextMonthInd);
                       CHKERRQ(ierr);
            mytemp = getTemperatureFromMonthlyData(
                       smonthtemp[currMonthInd], smonthtemp[nextMonthInd], 
                       i, j, (PetscScalar) day);
          } else {
            mytemp = getTemperatureFromYearlyCycle(summer_warming, mean_annual,
                                                   (PetscScalar) day);
          }
          const double randadd = gsl_ran_gaussian(pddRandGen, pddStdDev);
          const PetscScalar temp = mytemp + (PetscScalar) randadd;
          if (temp > 0.0)   pdd_sum += temp;
          //if ((i == grid->Mx/2) && (j == grid->My/2)) {
          //  ierr = PetscPrintf(grid->com,
          //    "  day=%d: mytemp=%5.4f, randadd=%5.4f, temp=%5.4f, pdd_sum=%5.4f\n",
          //    day, mytemp, randadd, temp, pdd_sum); CHKERRQ(ierr);
          //}
        }
      } else { // default Calov-Greve method; apply Simpson's rule for integral
        pdd_sum = 0.0;
        for (PetscInt m = 0; m < CGsumcount; ++m) {
          // Simpson's is: (h/3) * sum([1 4 2 4 2 4 ... 4 1] .* [f(x0) f(x1) ... f(xN)])
          PetscScalar  coeff = ((m % 2) == 1) ? 4.0 : 2.0;
          if ( (m == 0) || (m == (CGsumcount - 1)) )  coeff = 1.0;
          const PetscScalar day = startday + m * CGsumstepdays;
          PetscScalar temp;
          if (usingMonthlyTemps == PETSC_TRUE) {
            PetscInt currMonthInd, nextMonthInd;
            ierr = getMonthIndicesFromDay(day, currMonthInd, nextMonthInd); CHKERRQ(ierr);
            temp = getTemperatureFromMonthlyData(
                       smonthtemp[currMonthInd], smonthtemp[nextMonthInd], 
                       i, j, day);
          } else {
            temp = getTemperatureFromYearlyCycle(summer_warming, mean_annual,day);
          }
          pdd_sum += coeff * CalovGreveIntegrand(temp);
        }
        pdd_sum = (CGsumstepdays / 3.0) * pdd_sum;
      }

      // now that we have the number of PDDs, compute mass balance from snow rate
      smflux[i][j] = getSurfaceBalanceFromSnowAndPDD(
                            saccum[i][j], dt_years * secpera, pdd_sum);

      //if ((i == grid->Mx/2) && (j == grid->My/2)) {
      //  ierr = PetscPrintf(grid->com,
      //    "\nPDD at (i,j)=(%d,%d):\n"
      //      "  mean_annual=%5.4f, summer_warming=%5.4f\n",
      //      i, j, mean_annual, summer_warming); CHKERRQ(ierr);
      //  ierr = PetscPrintf(grid->com,
      //      "  CGsumcount = %d, CGsumstepdays = %5.2f, num_days = %d, intstartday = %d,\n"
      //      "  h = %6.2f, pdd_sum = %5.2f, saccum = %5.4f, smflux = %5.4f\n",
      //      CGsumcount, CGsumstepdays, num_days, intstartday,
      //      h[i][j], pdd_sum, saccum[i][j]*secpera, smflux[i][j]*secpera); CHKERRQ(ierr);
      //}

    }
  }

  if (usingMonthlyTemps == PETSC_TRUE) {
    for (PetscInt j = 0; j < 12; ++j) {
      if (vmonthlysurftemp[j].was_created()) {
        ierr = vmonthlysurftemp[j].end_access(); CHKERRQ(ierr);
      } else {
        SETERRQ1(2,"vmonthlysurftemp[%d] not created",j);
      }
    }
  }

  ierr = info->surfelev->end_access(); CHKERRQ(ierr);
  ierr = info->lat->end_access(); CHKERRQ(ierr);

  ierr =    vsurfmassflux.end_access(); CHKERRQ(ierr);
  ierr = vannmeansurftemp.end_access(); CHKERRQ(ierr);
  ierr =       vsurfaccum.end_access(); CHKERRQ(ierr);
  
  // now that it is up to date, return pointer to it
  pvsmf = &vsurfmassflux;
  return 0;
}


PetscErrorCode PISMPDDCoupler::updateSurfTempAndProvide(
             const PetscScalar t_years, const PetscScalar dt_years,
             void *iceInfoNeeded, // will be interpreted as type iceInfoNeededByAtmosphereCoupler*
             IceModelVec2* &pvst) {

  PetscErrorCode ierr;

  IceInfoNeededByAtmosphereCoupler* info = (IceInfoNeededByAtmosphereCoupler*) iceInfoNeeded;

  // need surface temperature representative of period [t_years,t_years+dt_years]
  //   get this by using midpoint of this interval
  const PetscScalar mid_years = t_years + (dt_years/2.0),
                    mid_day = 365.24 * (mid_years - floor(mid_years));
  PetscScalar **stemp;

  ierr = vsurftemp.get_array(stemp); CHKERRQ(ierr);

  if (usingMonthlyTemps == PETSC_TRUE) {
    PetscScalar **smonthtemp[12];
  
    for (PetscInt j = 0; j < 12; ++j) {
      if (vmonthlysurftemp[j].was_created()) {
        ierr = vmonthlysurftemp[j].get_array(smonthtemp[j]); CHKERRQ(ierr);
      } else {
        SETERRQ1(1,"vmonthlysurftemp[%d] not created",j);
      }
    }
  
    PetscInt currMonthInd, nextMonthInd;
    ierr = getMonthIndicesFromDay(mid_day, currMonthInd, nextMonthInd); CHKERRQ(ierr);
    for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
      for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
        stemp[i][j] = getTemperatureFromMonthlyData(
                       smonthtemp[currMonthInd], smonthtemp[nextMonthInd], i, j, mid_day);
      }
    }

    for (PetscInt j = 0; j < 12; ++j) {
      if (vmonthlysurftemp[j].was_created()) {
        ierr = vmonthlysurftemp[j].end_access(); CHKERRQ(ierr);
      } else {
        SETERRQ1(2,"vmonthlysurftemp[%d] not created",j);
      }
    }
  } else {
    // no monthly temps; usual case where surface temperature comes from stored
    //   annual mean (vannmeansurftemp) with additional yearly cycle controlled
    //   by summer_warming
    PetscScalar **h, **lat, **amstemp;
    ierr = info->surfelev->get_array(h);   CHKERRQ(ierr);
    ierr = info->lat->get_array(lat); CHKERRQ(ierr);
    ierr = vannmeansurftemp.get_array(amstemp); CHKERRQ(ierr);
    for (PetscInt i = grid->xs; i<grid->xs+grid->xm; ++i) {
      for (PetscInt j = grid->ys; j<grid->ys+grid->ym; ++j) {
        const PetscScalar mean_annual = amstemp[i][j] - 273.15,  // in deg C
                          summer_warming = getSummerWarming(h[i][j],lat[i][j],mean_annual);
        stemp[i][j] = getTemperatureFromYearlyCycle(summer_warming, mean_annual, mid_day);
      }
    }
    ierr = info->surfelev->end_access(); CHKERRQ(ierr);
    ierr = info->lat->end_access(); CHKERRQ(ierr);
    ierr = vannmeansurftemp.end_access(); CHKERRQ(ierr);
  }

  ierr = vsurftemp.end_access(); CHKERRQ(ierr);

  pvst = &vsurftemp;
  return 0;
}

