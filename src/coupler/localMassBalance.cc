// Copyright (C) 2009 Ed Bueler
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

#include <petsc.h>
#include <ctime>  // for time(), used to initialize random number gen
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>       // for erfc() in CalovGreveIntegrand()
#include "../base/pism_const.hh"
#include "../base/NCVariable.hh"
#include "localMassBalance.hh"


LocalMassBalance::LocalMassBalance(NCConfigVariable* myconfig) {
  config = myconfig;
}


PetscErrorCode LocalMassBalance::init() {
  PetscPrintf(PETSC_COMM_WORLD,"LocalMassBalance is a virtual class.  ENDING ...\n");
  PetscEnd();
  return 0;
}


PetscErrorCode LocalMassBalance::getNForTemperatureSeries(
                   PetscScalar /* t */, PetscScalar /* dt */, PetscInt & /* N */) {
  PetscPrintf(PETSC_COMM_WORLD,"LocalMassBalance is a virtual class.  ENDING ...\n");
  PetscEnd();
  return 0;
}


PetscScalar LocalMassBalance::getMassFluxFromTemperatureTimeSeries(
             PetscScalar /* t */, PetscScalar /* dt_series */, PetscScalar * /* T */, PetscInt /* N */,
             PetscScalar /* precip */) {
  PetscPrintf(PETSC_COMM_WORLD,"LocalMassBalance is a virtual class.  ENDING ...\n");
  PetscEnd();
  return 0.0;
}


PDDMassBalance::PDDMassBalance(NCConfigVariable* myconfig) : LocalMassBalance(myconfig) {
  // FIXME: switch over scheme and defaults to Fausto choice; make EISMINT-Greenland a special case,
  //   but not needing a derived class (I think)
  pddFactorSnow   = config->get("pdd_factor_snow");
  pddFactorIce    = config->get("pdd_factor_ice");
  pddRefreezeFrac = config->get("pdd_refreeze");
  pddStdDev       = config->get("pdd_std_dev");
}


PetscErrorCode PDDMassBalance::init() {
  // check options for parameter values
  PetscErrorCode ierr;
  PetscTruth     pSet;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_factor_snow", &pddFactorSnow, &pSet);
             CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_factor_ice", &pddFactorIce, &pSet);
             CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_refreeze", &pddRefreezeFrac, &pSet);
             CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-pdd_std_dev", &pddStdDev, &pSet);
             CHKERRQ(ierr);
  return 0;
}


PetscErrorCode PDDMassBalance::setDegreeDayFactorsFromSpecialInfo(
                                  PetscScalar latitude, PetscScalar T_mj) {
  const PetscScalar
    beta_ice_w = config->get("pdd_fausto_beta_ice_w"),
    beta_snow_w = config->get("pdd_fausto_beta_snow_w");
  if (latitude < config->get("pdd_fausto_latitude_beta_w")) { // case lat < 72 deg N
    pddFactorIce  = beta_ice_w;
    pddFactorSnow = beta_snow_w;
  } else { // case > 72 deg N
    const PetscScalar  
      T_c = config->get("pdd_fausto_T_c"),
      T_w = config->get("pdd_fausto_T_w"),
      beta_ice_c = config->get("pdd_fausto_beta_ice_c"),
      beta_snow_c = config->get("pdd_fausto_beta_snow_c");
    if (T_mj >= T_w) {
      pddFactorIce  = beta_ice_w;
      pddFactorSnow = beta_snow_w;
    } else if (T_mj <= T_c) {
      pddFactorIce  = beta_ice_c;
      pddFactorSnow = beta_snow_c;
    } else { // middle case   T_c < T_mj < T_w
      const PetscScalar
         lam_i = pow( (T_w - T_mj) / (T_w - T_c) , 3.0),
         lam_s = (T_mj - T_c) / (T_w - T_c);
      pddFactorIce  = beta_ice_w + (beta_ice_c - beta_ice_w) * lam_i;
      pddFactorSnow = beta_snow_w + (beta_snow_c - beta_snow_w) * lam_s;
    }
  }

  // degree-day factors in \ref Faustoetal2009 are water-equivalent
  //   thickness per degree day; ice-equivalent thickness melted per degree
  //   day is slightly larger; for example, iwfactor = 1000/910
  const PetscScalar iwfactor = config->get("fresh_water_density") / config->get("ice_density");
  pddFactorSnow *= iwfactor;
  pddFactorIce  *= iwfactor;
  return 0;
}


/*!
Because Calov-Greve method uses Simpson's rule to do integral, we choose the number 
of times to be odd.  Numerical integration accuracy, assuming a smooth-in-time yearly cycle,
suggests that 53 evals per year (i.e. approximately weekly) should be sufficiently accurate.
 */
PetscErrorCode PDDMassBalance::getNForTemperatureSeries(
                   PetscScalar /* t */, PetscScalar dt, PetscInt &N) {
  PetscScalar dt_years = dt / secpera;
  N = (int) ceil(52 * (dt_years) + 1);
  if (N < 3) N = 3;
  if ((N % 2) == 0)  N++;  // guarantee it is odd
  return 0;
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

The scheme here came from EISMINT-Greenland \ref RitzEISMINT .

Arguments are snow fall rate precip in (ice-equivalent) m s-1, t and dt in s, T[] in
K.  There are N temperature values T[0],...,T[N-1] at times
t,t+dt_series,...,t+(N-1)*dt_series.  This time series covers the interval [t,t+dt]
where dt = (N-1)*dt_series.
 */
PetscScalar PDDMassBalance::getMassFluxFromTemperatureTimeSeries(
             PetscScalar t, PetscScalar dt_series, PetscScalar *T, PetscInt N,
             PetscScalar precip) {
#if 0
  static bool flag = true;
  if (flag) {
     PetscPrintf(PETSC_COMM_WORLD,
     "\nfirst call to PDDMassBalance::getMassFluxFromTemperatureTimeSeries():\n"
       "  pddFactorSnow,pddFactorIce,pddRefreezeFrac,pddStdDev=%f,%f,%f,%f\n\n",
     pddFactorSnow,pddFactorIce,pddRefreezeFrac,pddStdDev); 
     flag = false;
  }
#endif
  const PetscScalar
      pddsum = getPDDSumFromTemperatureTimeSeries(t,dt_series,T,N), // units: K day
      dt     = (N-1) * dt_series; 
  if (precip < 0.0) {
    // neg precip interpreted as ablation, so positive degree-days are ignored
    return precip;
  } else {
    // positive precip: it snowed (precip = snow; never rain)
    const PetscScalar snow        = precip * dt,   // units: m (ice-equivalent)
                      snow_melted = pddsum * pddFactorSnow;  // units: m (ice-equivalent)
    if (snow_melted <= snow) {
      return ((snow - snow_melted) + (snow_melted * pddRefreezeFrac)) / dt;
    } else { // it is snowing, but all the snow melts and refreezes; this ice is
             // then removed, plus possibly more of the underlying ice
      const PetscScalar ice_deposited = snow * pddRefreezeFrac,
                        excess_pddsum = pddsum - (snow / pddFactorSnow), // positive!; units K day
                        ice_melted    = excess_pddsum * pddFactorIce; // units: K day
      return (ice_deposited - ice_melted) / dt;
    }
  }
}


//! Compute the integrand in integral (6) in \ref CalovGreve05 .
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
\c dt instead of a whole year as stated in \ref CalovGreve05 .

The argument \c Tac is the temperature in K.  The value \f$T_{ac}(t)\f$
in the above integral must be in degrees C, so the shift is done within this 
procedure.
 */
PetscScalar PDDMassBalance::CalovGreveIntegrand(
               PetscScalar sigma, PetscScalar Tac) {
  const PetscScalar TacC = Tac - 273.15;
  return (sigma / sqrt(2.0 * pi)) * exp(- TacC * TacC / (2.0 * sigma * sigma)) 
           + (TacC / 2.0) * gsl_sf_erfc(- TacC / (sqrt(2.0) * sigma));
}


PetscScalar PDDMassBalance::getPDDSumFromTemperatureTimeSeries(
             PetscScalar /* t */, PetscScalar dt_series, PetscScalar *T, PetscInt N) {
  PetscScalar  pdd_sum = 0.0;  // return value has units  K day
  const PetscScalar sperd = 8.64e4, // exact seconds per day
                    h_days = dt_series / sperd;
  const PetscInt Nsimp = ((N % 2) == 1) ? N : N-1; // odd N case is pure simpson's
  // Simpson's rule is:
  //   integral \approx (h/3) * sum( [1 4 2 4 2 4 ... 4 1] .* [f(t_0) f(t_1) ... f(t_N-1)] )
  for (PetscInt m = 0; m < Nsimp; ++m) {
    PetscScalar  coeff = ((m % 2) == 1) ? 4.0 : 2.0;
    if ( (m == 0) || (m == (Nsimp-1)) )  coeff = 1.0;
    pdd_sum += coeff * CalovGreveIntegrand(pddStdDev,T[m]);  // pass in temp in K
  }
  pdd_sum = (h_days / 3.0) * pdd_sum;
  if (Nsimp < N) { // add one more subinterval by trapezoid
    pdd_sum += (h_days / 2.0) * ( CalovGreveIntegrand(pddStdDev,T[N-2])
                                  + CalovGreveIntegrand(pddStdDev,T[N-1]) );
  }
#if 0
  static bool flag = true;
  if (flag) {
     PetscPrintf(PETSC_COMM_WORLD,
     "\nfirst call to PDDMassBalance::getPDDSumFromTemperatureTimeSeries():\n"
       "  pdd_sum=%f\n\n",pdd_sum); 
     flag = false;
  }
#endif
  return pdd_sum;
}


/*!
Initializes the random number generator (RNG).  The RNG is GSL's recommended default,
which seems to be "mt19937" and is DIEHARD (whatever that means ...). Seed with
wall clock time in seconds in non-repeatable case, and with 0 in repeatable case.
 */
PDDrandMassBalance::PDDrandMassBalance(NCConfigVariable* myconfig, bool repeatable) : PDDMassBalance(myconfig) {
  pddRandGen = gsl_rng_alloc(gsl_rng_default);  // so pddRandGen != NULL now
  gsl_rng_set(pddRandGen, repeatable ? 0 : time(0)); 
#if 0
  PetscTruth     pSet;
  PetscOptionsGetScalar(PETSC_NULL, "-pdd_std_dev", &pddStdDev, &pSet);
  PetscPrintf(PETSC_COMM_WORLD,"\nPDDrandMassBalance constructor; pddStdDev = %10.4f\n",
              pddStdDev);
  for (int k=0; k<20; k++) {
    PetscPrintf(PETSC_COMM_WORLD,"  %9.4f\n",gsl_ran_gaussian(pddRandGen, pddStdDev));
  }
  PetscPrintf(PETSC_COMM_WORLD,"\n\n");
#endif
}


PDDrandMassBalance::~PDDrandMassBalance() {
  if (pddRandGen != NULL) {
    gsl_rng_free(pddRandGen);
    pddRandGen = NULL;
  }
}


/*!
We need to compute simulated random temperature each actual \e day, or at least as
close as we can reasonably get.  Output Nseries is number of days plus one.

Implementation of getPDDSumFromTemperatureTimeSeries() does require Nseries >= 2,
however.

(Alternatively we could adjust the standard deviation
to account for non-day sub-intervals.)
 */
PetscErrorCode PDDrandMassBalance::getNForTemperatureSeries(
                   PetscScalar /* t */, PetscScalar dt, PetscInt &Nseries) {
  const PetscScalar sperd = 8.64e4;
  Nseries = (int) ceil(dt / sperd);
  if (Nseries < 2) Nseries = 2;
  return 0;
}


PetscScalar PDDrandMassBalance::getPDDSumFromTemperatureTimeSeries(
             PetscScalar /* t */, PetscScalar dt_series, PetscScalar *T, PetscInt N) {
  PetscScalar       pdd_sum = 0.0;  // return value has units  K day
  const PetscScalar sperd = 8.64e4, // exact seconds per day
                    h_days = dt_series / sperd;
  // there are N-1 intervals [t,t+dt],...,[t+(N-2)dt,t+(N-1)dt]
  for (PetscInt m = 0; m < N-1; ++m) {
    PetscScalar temp = 0.5 * (T[m] + T[m+1]); // av temp in [t+m*dt,t+(m+1)*dt]
    temp += gsl_ran_gaussian(pddRandGen, pddStdDev); // add random: N(0,sigma)
    if (temp > 273.15)   pdd_sum += h_days * (temp - 273.15);
  }
#if 0
  static bool flag = true;
  if (flag) {
     PetscPrintf(PETSC_COMM_WORLD,
     "\nfirst call to PDDrandMassBalance::getPDDSumFromTemperatureTimeSeries():\n"
       "  pdd_sum=%f\n\n",pdd_sum); 
     flag = false;
  }
#endif
  return pdd_sum;
}

