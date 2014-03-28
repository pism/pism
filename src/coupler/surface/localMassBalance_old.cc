// Copyright (C) 2009, 2010, 2011, 2014 Ed Bueler and Constantine Khroulev and Andy Aschwanden
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
#include "pism_const.hh"
#include "NCVariable.hh"
#include "localMassBalance_old.hh"
#include "IceGrid.hh"
#include "PISMConfig.hh"
#include <cmath>

PDDMassBalance_Old::PDDMassBalance_Old(const PISMConfig& myconfig) : LocalMassBalance_Old(myconfig) {
  precip_as_snow = config.get_flag("interpret_precip_as_snow");
  Tmin = config.get("air_temp_all_precip_as_snow");
  Tmax = config.get("air_temp_all_precip_as_rain");

  secpera = 3.15569259747e7;
}


/*!
Because Calov-Greve method uses Simpson's rule to do integral, we choose the 
returned number of times N to be odd.
 */
PetscErrorCode PDDMassBalance_Old::getNForTemperatureSeries(PetscScalar /* t */,
                                                        PetscScalar dt, PetscInt &N) {
  PetscInt    NperYear = static_cast<PetscInt>( 
                           config.get("pdd_max_evals_per_year") );
  PetscScalar dt_years = dt / secpera;
  N = (int) ceil((NperYear - 1) * (dt_years) + 1);
  if (N < 3) N = 3;
  if ((N % 2) == 0)  N++;  // guarantee N is odd
  return 0;
}


//! Compute the integrand in integral (6) in [\ref CalovGreve05].
/*!
The integral is
   \f[\mathrm{PDD} = \int_{t_0}^{t_0+\mathtt{dt}} dt\,
         \bigg[\frac{\sigma}{\sqrt{2\pi}}\,\exp\left(-\frac{T_{ac}(t)^2}{2\sigma^2}\right)
               + \frac{T_{ac}(t)}{2}\,\mathrm{erfc}
               \left(-\frac{T_{ac}(t)}{\sqrt{2}\,\sigma}\right)\bigg] \f]
This procedure computes the quantity in square brackets.  The value \f$T_{ac}(t)\f$
in the above integral is in degrees C.  Here we think of the argument \c TacC
as temperature in Celsius, but really it is the temperature above a threshold
at which it is "positive".

This integral is used for the expected number of positive degree days, unless the
user selects a random PDD implementation with <tt>-pdd_rand</tt> or 
<tt>-pdd_rand_repeatable</tt>.  The user can choose \f$\sigma\f$ by option
<tt>-pdd_std_dev</tt>.  Note that the integral is over a time interval of length
\c dt instead of a whole year as stated in \ref CalovGreve05 .
 */
PetscScalar PDDMassBalance_Old::CalovGreveIntegrand(PetscScalar sigma, PetscScalar TacC) {
  const PetscScalar Z    = TacC / (sqrt(2.0) * sigma);
  return (sigma / sqrt(2.0 * M_PI)) * exp(-Z*Z) + (TacC / 2.0) * gsl_sf_erfc(-Z);
}


//! Compute the expected number of positive degree days from the input temperature time-series.
PetscScalar PDDMassBalance_Old::getPDDSumFromTemperatureTimeSeries(
               PetscScalar pddStdDev, PetscScalar pddThresholdTemp,
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
    pdd_sum += coeff * CalovGreveIntegrand(pddStdDev,T[m]-pddThresholdTemp);  // pass in temp in K
  }
  pdd_sum = (h_days / 3.0) * pdd_sum;
  if (Nsimp < N) { // add one more subinterval by trapezoid
    pdd_sum += (h_days / 2.0) * ( CalovGreveIntegrand(pddStdDev,T[N-2]-pddThresholdTemp)
                                  + CalovGreveIntegrand(pddStdDev,T[N-1]-pddThresholdTemp) );
  }
  return pdd_sum;
}


//! /brief Report the amount of snow fallen in the given time period, according
//!        to the temperature time-series; remove the rain.
/*! 
Use the temperature time-series to determine whether the
precipitation is snow or rain.  Rain is removed entirely from the surface mass
balance, and will not be included in the computed runoff, which is meltwater runoff.
There is an allowed linear transition for Tmin below which all
precipitation is interpreted as snow, and Tmax above which all precipitation is
rain (see, e.g. \ref Hock2005b).
 */
PetscScalar PDDMassBalance_Old::getSnowFromPrecipAndTemperatureTimeSeries(
                 PetscScalar precip_rate,
                 PetscScalar /*t*/, PetscScalar dt_series, PetscScalar *T, PetscInt N) {

  PetscScalar snow = 0.0;
  if (precip_as_snow) {
    // positive precip_rate: it snowed (precip = snow; never rain)
    snow = precip_rate * (N-1) * dt_series;   // units: m (ice-equivalent)
  } else {
    // Following \ref Hock2005b we employ a linear transition from Tmin to Tmax
    for (PetscInt i=0; i<N-1; i++) { // go over all N-1 subintervals in interval[t,t+dt_series]
      const PetscScalar Tav = (T[i] + T[i+1]) / 2.0; // use midpt of temp series for subinterval
      PetscScalar acc_rate = 0.0;                    // accumulation rate during a sub-interval

      if (Tav <= Tmin) { // T <= Tmin, all precip is snow
        acc_rate = precip_rate;
      } else if (Tav < Tmax) { // linear transition from Tmin to Tmax
        acc_rate = ((Tmax-Tav)/(Tmax-Tmin)) * precip_rate;
      } else { // T >= Tmax, all precip is rain -- ignore it
        acc_rate = 0;
      }
      snow += acc_rate * dt_series;  // units: m (ice-equivalent)
    }
  }
  return snow;
}


//! \brief Compute the surface mass balance at a location from the number of positive
//! degree days and the accumulation amount in a time interval.
/*!
Several mass fluxes, including the final surface mass balance, are computed 
as ice-equivalent meters per second from the number of positive degree days
and an accumulation amount in the time interval.

This is a PDD scheme.  The input parameter \c ddf.snow is a rate of melting per
positive degree day for snow.  If the number of PDDs, which describes an energy
available for melt, exceeds those needed to melt all of the accumulation (snow, using
\c ddf.snow) then:
  \li a fraction of the melted snow refreezes, conceptualized as superimposed
      ice, and this is controlled by parameter \c ddf.refreezeFrac, and
  \li the excess number of PDDs is used to melt both the ice that came from
      refreeze and then any ice which is already present.
    
In either case, \e ice melts at a constant rate per positive degree day,
controlled by parameter \c ddf.ice.

In the weird case where the accumulation amount is negative, it is interpreted
as an (ice-equivalent) direct ablation rate.  Accumulation amounts are generally
positive everywhere on ice sheets, however.

The scheme here came from EISMINT-Greenland [\ref RitzEISMINT], but is influenced
by R. Hock (personal communication).

The input argument \c dt is in seconds.  The input argument \c pddsum is in K day.
These strange units work because the degree day factors in \c ddf have their
usual glaciological units, namely m K-1 day-1, and we only multiply degree 
day factors times \e amounts of snow, not rates.

There are four outputs
  \li \c accumulation_rate,
  \li \c melt_rate,
  \li \c runoff_rate,
  \li \c smb_rate.

All are rates and all are in ice-equivalent m s-1.  Note "accumulation"
is the part of precipitation which is not rain, so, if needed, the modeled
rain rate is <tt>rain_rate = precip_rate - accum_rate</tt>.  Again, "runoff" is 
meltwater runoff and does not include rain.

In normal areas where \c accumulation > 0, the output quantities satisfy
  \li <tt>smb_rate = accumulation_rate - runoff_rate</tt>

in all cases.
 */
PetscErrorCode PDDMassBalance_Old::getMassFluxesFromPDDs(const DegreeDayFactors_Old &ddf,
                                                     PetscScalar dt, PetscScalar pddsum,
                                                     PetscScalar accumulation,
                                                     PetscScalar &accumulation_rate,
                                                     PetscScalar &melt_rate,
                                                     PetscScalar &runoff_rate,
                                                     PetscScalar &smb_rate) {

#if (PISM_DEBUG==1)
  if (dt <= 0) {
    SETERRQ1(PETSC_COMM_SELF, 1, "dt = %f is not alloved", dt);
  }
#endif

  if (accumulation < 0.0) {           // weird, but allowed, case
    accumulation_rate = 0.0;
    melt_rate         = accumulation / dt;
    runoff_rate       = melt_rate;
    smb_rate          = melt_rate;
    return 0;
  }

  accumulation_rate = accumulation / dt;
  
  // no melt case: we're done
  if (pddsum <= 0.0) {
    melt_rate   = 0.0;
    runoff_rate = 0.0;
    smb_rate    = accumulation_rate; // = accumulation_rate - runoff_rate
    return 0;
  }

  const PetscScalar snow_max_melted = pddsum * ddf.snow;  // units: m (ice-equivalent)
  if (snow_max_melted <= accumulation) {
    // some of the snow melted and some is left; in any case, all of the energy
    //   available for melt, namely all of the positive degree days (PDDs) were
    //   used-up in melting snow; but because the snow does not completely melt,
    //   we generate 0 modeled runoff, and the surface balance is the accumulation
    melt_rate   = snow_max_melted / dt;
    runoff_rate = 0.0;
    smb_rate    = accumulation_rate; // = accumulation_rate - runoff_rate
    return 0;
  }

  // at this point: all of the snow melted; some of melt mass is kept as 
  //   refrozen ice; excess PDDs remove some of this ice and perhaps more of the
  //   underlying ice
  const PetscScalar  ice_created_by_refreeze = accumulation * ddf.refreezeFrac;  // units: m (ice-equivalent)

  const PetscScalar
      excess_pddsum = pddsum - (accumulation / ddf.snow), // units: K day
      ice_melted    = excess_pddsum * ddf.ice; // units: m (ice-equivalent)

  // all of the snow, plus some of the ice, was melted
  melt_rate   = (accumulation + ice_melted) / dt;
  // all of the snow which melted but which did not refreeze, plus all of the
  //   ice which melted (ice includes refrozen snow!), combine to give meltwater
  //   runoff rate 
  runoff_rate = (accumulation - ice_created_by_refreeze + ice_melted) / dt;
  smb_rate    = accumulation_rate - runoff_rate;
  return 0;
}


/*!
Initializes the random number generator (RNG).  The RNG is GSL's recommended default,
which seems to be "mt19937" and is DIEHARD (whatever that means ...). Seed with
wall clock time in seconds in non-repeatable case, and with 0 in repeatable case.
 */
PDDrandMassBalance_Old::PDDrandMassBalance_Old(const PISMConfig& myconfig, bool repeatable)
    : PDDMassBalance_Old(myconfig) {
  pddRandGen = gsl_rng_alloc(gsl_rng_default);  // so pddRandGen != NULL now
  gsl_rng_set(pddRandGen, repeatable ? 0 : time(0));
}


PDDrandMassBalance_Old::~PDDrandMassBalance_Old() {
  if (pddRandGen != NULL) {
    gsl_rng_free(pddRandGen);
    pddRandGen = NULL;
  }
}


/*!
We need to compute simulated random temperature each actual \e day, or at least as
close as we can reasonably get.  Output \c N is number of days or number of days
plus one.

Thus this method ignores <tt>config.get("pdd_max_temperature_evals_per_year")</tt>,
which is used in the base class PDDMassBalance_Old.

Implementation of getPDDSumFromTemperatureTimeSeries() requires returned
N >= 2, so we guarantee that.
 */
PetscErrorCode PDDrandMassBalance_Old::getNForTemperatureSeries(
                   PetscScalar /* t */, PetscScalar dt, PetscInt &N) {
  const PetscScalar sperd = 8.64e4;
  N = (int) ceil(dt / sperd);
  if (N < 2) N = 2;
  return 0;
}


PetscScalar PDDrandMassBalance_Old::getPDDSumFromTemperatureTimeSeries(
             PetscScalar pddStdDev, PetscScalar pddThresholdTemp,
             PetscScalar /* t */, PetscScalar dt_series, PetscScalar *T, PetscInt N) {
  PetscScalar       pdd_sum = 0.0;  // return value has units  K day
  const PetscScalar sperd = 8.64e4, // exact seconds per day
                    h_days = dt_series / sperd;
  // there are N-1 intervals [t,t+dt],...,[t+(N-2)dt,t+(N-1)dt]
  for (PetscInt m = 0; m < N-1; ++m) {
    PetscScalar temp = 0.5*(T[m] + T[m+1]); // av temp in [t+m*dt,t+(m+1)*dt]
    temp += gsl_ran_gaussian(pddRandGen, pddStdDev); // add random: N(0,sigma)
    if (temp > pddThresholdTemp)
      pdd_sum += h_days * (temp - pddThresholdTemp);
  }
  return pdd_sum;
}


FaustoGrevePDDObject_Old::FaustoGrevePDDObject_Old(IceGrid &g, const PISMConfig &myconfig)
  : grid(g), config(myconfig) {

  beta_ice_w = config.get("pdd_fausto_beta_ice_w");
  beta_snow_w = config.get("pdd_fausto_beta_snow_w");

  T_c = config.get("pdd_fausto_T_c");
  T_w = config.get("pdd_fausto_T_w");
  beta_ice_c = config.get("pdd_fausto_beta_ice_c");
  beta_snow_c = config.get("pdd_fausto_beta_snow_c");

  fresh_water_density = config.get("fresh_water_density");
  ice_density = config.get("ice_density");
  pdd_fausto_latitude_beta_w = config.get("pdd_fausto_latitude_beta_w");

  temp_mj.create(grid, "temp_mj_faustogreve", WITHOUT_GHOSTS);
  temp_mj.set_attrs("internal",
                    "mean July air temp from Fausto et al (2009) parameterization",
                    "K", "");
}


PetscErrorCode FaustoGrevePDDObject_Old::setDegreeDayFactors(PetscInt i, PetscInt j,
                                                         PetscScalar /* usurf */,
                                                         PetscScalar lat, PetscScalar /* lon */,
                                                         DegreeDayFactors_Old &ddf) {

  PetscErrorCode ierr = temp_mj.begin_access(); CHKERRQ(ierr);
  const PetscScalar T_mj = temp_mj(i,j);
  ierr = temp_mj.end_access(); CHKERRQ(ierr);

  if (lat < pdd_fausto_latitude_beta_w) { // case lat < 72 deg N
    ddf.ice  = beta_ice_w;
    ddf.snow = beta_snow_w;
  } else { // case > 72 deg N
    if (T_mj >= T_w) {
      ddf.ice  = beta_ice_w;
      ddf.snow = beta_snow_w;
    } else if (T_mj <= T_c) {
      ddf.ice  = beta_ice_c;
      ddf.snow = beta_snow_c;
    } else { // middle case   T_c < T_mj < T_w
      const PetscScalar
         lam_i = pow( (T_w - T_mj) / (T_w - T_c) , 3.0),
         lam_s = (T_mj - T_c) / (T_w - T_c);
      ddf.ice  = beta_ice_w + (beta_ice_c - beta_ice_w) * lam_i;
      ddf.snow = beta_snow_w + (beta_snow_c - beta_snow_w) * lam_s;
    }
  }

  // degree-day factors in \ref Faustoetal2009 are water-equivalent
  //   thickness per degree day; ice-equivalent thickness melted per degree
  //   day is slightly larger; for example, iwfactor = 1000/910
  const PetscScalar iwfactor = fresh_water_density / ice_density;
  ddf.snow *= iwfactor;
  ddf.ice  *= iwfactor;
  return 0;
}


//! Updates mean July near-surface air temperature.
/*!
Unfortunately this duplicates code in PA_SeaRISE_Greenland::update();
 */
PetscErrorCode FaustoGrevePDDObject_Old::update_temp_mj(
    IceModelVec2S *surfelev, IceModelVec2S *lat, IceModelVec2S *lon) {
  PetscErrorCode ierr;

  const PetscReal 
    d_mj     = config.get("snow_temp_fausto_d_mj"),      // K
    gamma_mj = config.get("snow_temp_fausto_gamma_mj"),  // K m-1
    c_mj     = config.get("snow_temp_fausto_c_mj"),      // K (degN)-1
    kappa_mj = config.get("snow_temp_fausto_kappa_mj");  // K (degW)-1
  
  PetscScalar **lat_degN, **lon_degE, **h;
  ierr = surfelev->get_array(h);   CHKERRQ(ierr);
  ierr = lat->get_array(lat_degN); CHKERRQ(ierr);
  ierr = lon->get_array(lon_degE); CHKERRQ(ierr);
  ierr = temp_mj.begin_access();  CHKERRQ(ierr);

  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      temp_mj(i,j) = d_mj + gamma_mj * h[i][j] + c_mj * lat_degN[i][j] + kappa_mj * (-lon_degE[i][j]);
    }
  }
  
  ierr = surfelev->end_access();   CHKERRQ(ierr);
  ierr = lat->end_access(); CHKERRQ(ierr);
  ierr = lon->end_access(); CHKERRQ(ierr);
  ierr = temp_mj.end_access();  CHKERRQ(ierr);

  return 0;
}

