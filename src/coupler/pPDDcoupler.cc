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
#include "pPDDcoupler.hh"


PISMPDDCoupler::PISMPDDCoupler() : PISMAtmosphereCoupler() {

  // normally use Calov-Greve expectation integral; if this is non-NULL
  //   then using actual random numbers
  pddRandGen = NULL;
  
  // see pPDDcoupler.hh for units and meaning
  pddFactorSnow = 0.003;  // EISMINT-Greenland value
  pddFactorIce = 0.008;  // EISMINT-Greenland value
  pddRefreezeFrac = 0.6;  // EISMINT-Greenland value
  pddStdDev = 5.0;  // EISMINT-Greenland value
  pddSummerWarming = 15.0;
     // re SUMMER_WARMING:  (30.38 - 0.006277 * 1000.0 - 0.3262 * 75.0)
     //                    - (49.13 - 0.007992 * 1000.0 -0.7576 * 75.0)
     //                   =  15.32   K
     // is result of EISMINT-Greenland formulas for h=1000.0 m and lat=75.0 deg N
  pddSummerPeakDay = 243.0;  // August 1st
}


PISMPDDCoupler::~PISMPDDCoupler() {
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
           "-pdd", "-pdd_rand", "-pdd_rand_repeatable",
           "-pdd_factor_snow", "-pdd_factor_ice", "-pdd_refreeze", 
           "-pdd_std_dev",
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
  char actualfilename[PETSC_MAX_PATH_LEN];
  char* filename=&(actualfilename[0]);
  LocalInterpCtx* lic;

  ierr = PISMAtmosphereCoupler::initFromOptions(g); CHKERRQ(ierr); // sets grid and metadata

  // mean annual ice equivalent accumulation rate
  ierr = vsurfaccum.create(*g, "acab", false); CHKERRQ(ierr);  // global; no ghosts
  ierr = vsurfaccum.set_attrs(
            "climate_state", 
            "mean annual ice equivalent accumulation rate",
	    "m s-1", 
	    NULL);  // no CF standard_name
	    CHKERRQ(ierr);
  ierr = vsurfaccum.set_glaciological_units("m year-1", secpera);
  ierr = vsurfaccum.set(0.0); CHKERRQ(ierr);  // merely a default value

  // now read two fields
  ierr = findPISMInputFile(&filename, lic); CHKERRQ(ierr); // allocates lic

  ierr = verbPrintf(2, g->com, 
     "initializing constant atmospheric climate and PDD: reading surface\n"
     "  accumulation 'acab' and surface temperature 'artm' from %s ... \n",
     filename); CHKERRQ(ierr); 

  ierr = vsurfaccum.regrid(filename, *lic, true); CHKERRQ(ierr); // it *is* critical
  ierr = vsurftemp.regrid(filename, *lic, true); CHKERRQ(ierr); // it *is* critical

  // now that they are read, reset vsurfaccum name to output name
  ierr = vsurfaccum.set_name("accum"); CHKERRQ(ierr);
  ierr = vsurfmassflux.set(0.0); CHKERRQ(ierr); // initialize to have some values

  // check on whether we should read monthly temperatures from file  
  char monthlyTempsFile[PETSC_MAX_PATH_LEN];
  ierr = PetscOptionsGetString(PETSC_NULL, "-pdd_monthly_temps", 
             monthlyTempsFile, PETSC_MAX_PATH_LEN, &pSet); CHKERRQ(ierr);
  if (pSet == PETSC_TRUE) {
    ierr = verbPrintf(2,grid->com,
       "using PDD based on monthly temperatures read from %s ...\n",
       monthlyTempsFile); CHKERRQ(ierr);
    ierr = readMonthlyTemps(monthlyTempsFile); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid->com,
       "using PDD based on standard yearly surface temp cycle ...\n"); CHKERRQ(ierr);
  }
  
  delete lic;
  return 0;
}


PetscErrorCode PISMPDDCoupler::writeCouplingFieldsToFile(const char *filename) {
  PetscErrorCode ierr;
  
  ierr = PISMAtmosphereCoupler::writeCouplingFieldsToFile(filename); CHKERRQ(ierr);
  
  ierr = vsurfaccum.write(filename, NC_FLOAT); CHKERRQ(ierr);

  for (PetscInt j = 0; j < 12; ++j) {
    if (vmonthlysurftemp[j].was_created()) {
      ierr = vmonthlysurftemp[j].write(filename, NC_FLOAT); CHKERRQ(ierr);
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
    char monthlyTempName[20], mTstring[80];
    snprintf(monthlyTempName, 20, "temp_mon%d", j);
    ierr = verbPrintf(2, grid->com, 
       "  reading month %d surface temperature '%s' ...\n", j, monthlyTempName); CHKERRQ(ierr); 
    ierr = vmonthlysurftemp[j].create(*grid, monthlyTempName, false); // global; no ghosts
       CHKERRQ(ierr);
    snprintf(mTstring, 80, "temperature at ice surface but below firn during month %d", j);
    ierr = vmonthlysurftemp[j].set_attrs(
            "climate_state",
            mTstring,
            "K", 
            NULL); CHKERRQ(ierr);
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
There are two versions:

  - if "-pdd_monthly_temps" then the monthly data read by 
    readMonthlyTempDataPDD() is linearly interpolated, while

  - if not "-pdd_monthly_temps" then a standard sinusoidal formula is used:
      \f[ T_s(d,i,j) = T_{ma} + S \cos((2\pi/365.24) * (d - P)) \f]
    where \f$d\f$ is the day on a 365.24 day calendar,
    \f$T_{ma}\f$ is the annual mean temperature in degrees C,
    \f$S\f$ is the amplitude of the sinusoid (the ``summer warming'') in degrees C,
    and \f$P\f$ is the day of peak summer warming, pddSummerPeakDay, which is rather
    arbitrarily taken to be August 1st by default.

EISMINT-Greenland assumes the second version \lo\cite{RitzEISMINT}\elo.  
See <tt>iceGRNModel.cc</tt>.

This procedure is \e virtual and replacable.
 */
PetscScalar PISMPDDCoupler::getTemperatureFromYearlyCycle(
       const PetscScalar summer_warming, const PetscScalar Tma, const PetscScalar day,
       const PetscInt i, const PetscInt j) {
  PetscTruth pSet;
  PetscOptionsHasName(PETSC_NULL, "-pdd_monthly_temps", &pSet);
  if (pSet == PETSC_FALSE) {
    // this is the default case for EISMINT-Greenland
    const PetscScalar  rad_per_day = 2.0 * PETSC_PI / 365.24;
    return Tma + summer_warming * cos(rad_per_day * (day - pddSummerPeakDay));
  }

  // this case is when "-pdd_monthly_temps foo.nc"
  PetscScalar month = 12.0 * day / 365.24;
  month = month - static_cast<PetscScalar> (((int) floor(month)) % 12);
  PetscInt  curr = (int) floor(month);
  curr = curr % 12;
  PetscInt  next = curr+1;
  if (next == 12)   next = 0;

  verbPrintf(1,grid->com,"NOT IMPLEMENTED");
  PetscEnd();

#if 0  
  PetscScalar **currTs, **nextTs;
  DAVecGetArray(grid.da2, vmonthlyTs[curr], &currTs);
  DAVecGetArray(grid.da2, vmonthlyTs[next], &nextTs);
  const PetscScalar myTs = (currTs[i][j] - 273.15) 
       + (month - (PetscScalar)curr) * (nextTs[i][j] - 273.15);
  DAVecRestoreArray(grid.da2, vmonthlyTs[curr], &currTs);
  DAVecRestoreArray(grid.da2, vmonthlyTs[next], &nextTs);
  return myTs;
#endif
  return 0.0;
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

