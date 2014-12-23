// Copyright (C) 2009, 2010, 2011, 2013, 2014 Ed Bueler and Constantine Khroulev and Andy Aschwanden
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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
#include <cmath>                // for erfc() in CalovGreveIntegrand()
#include <assert.h>
#include "pism_const.hh"
#include "PISMConfig.hh"
#include "localMassBalance.hh"
#include "IceGrid.hh"
#include <algorithm>

namespace pism {

LocalMassBalance::LocalMassBalance(const Config &myconfig)
  : config(myconfig), m_unit_system(config.get_unit_system()),
    seconds_per_day(86400) {
  // empty
}

PDDMassBalance::PDDMassBalance(const Config& myconfig)
  : LocalMassBalance(myconfig) {
  precip_as_snow     = config.get_flag("interpret_precip_as_snow");
  Tmin               = config.get("air_temp_all_precip_as_snow");
  Tmax               = config.get("air_temp_all_precip_as_rain");
  pdd_threshold_temp = config.get("pdd_positive_threshold_temp");
  refreeze_ice_melt  = config.get_flag("pdd_refreeze_ice_melt");
}


/*! \brief Compute the number of points for temperature and
    precipitation time-series.
 */
unsigned int PDDMassBalance::get_timeseries_length(double dt) {
  const unsigned int    NperYear = static_cast<unsigned int>(config.get("pdd_max_evals_per_year"));
  const double dt_years = m_unit_system.convert(dt, "seconds", "years");

  return std::max(1U, static_cast<unsigned int>(ceil(NperYear * dt_years)));
}


//! Compute the integrand in integral (6) in [\ref CalovGreve05].
/*!
The integral is
   \f[\mathrm{PDD} = \int_{t_0}^{t_0+\mathtt{dt}} dt\,
         \bigg[\frac{\sigma}{\sqrt{2\pi}}\,\exp\left(-\frac{T_{ac}(t)^2}{2\sigma^2}\right)
               + \frac{T_{ac}(t)}{2}\,\mathrm{erfc}
               \left(-\frac{T_{ac}(t)}{\sqrt{2}\,\sigma}\right)\bigg] \f]
This procedure computes the quantity in square brackets.  The value \f$T_{ac}(t)\f$
in the above integral is in degrees C.  Here we think of the argument `TacC`
as temperature in Celsius, but really it is the temperature above a threshold
at which it is "positive".

This integral is used for the expected number of positive degree days, unless the
user selects a random PDD implementation with `-pdd_rand` or
`-pdd_rand_repeatable`.  The user can choose \f$\sigma\f$ by option
`-pdd_std_dev`.  Note that the integral is over a time interval of length
`dt` instead of a whole year as stated in \ref CalovGreve05 . If `sigma` is zero, return the positive part of `TacC`.
 */
double PDDMassBalance::CalovGreveIntegrand(double sigma, double TacC) {

  if (sigma == 0) {
    return std::max(TacC, 0.0);
  } else {
    const double Z = TacC / (sqrt(2.0) * sigma);
    return (sigma / sqrt(2.0 * M_PI)) * exp(-Z*Z) + (TacC / 2.0) * erfc(-Z);
  }
}


//! Compute the expected number of positive degree days from the input temperature time-series.
/**
 * Use the rectangle method for simplicity.
 *
 * @param S standard deviation for air temperature excursions
 * @param dt_series length of the step for the time-series
 * @param T air temperature (array of length N)
 * @param N length of the T array
 * @param[out] PDDs pointer to a pre-allocated array with N-1 elements
 */
void PDDMassBalance::get_PDDs(double *S, double dt_series,
                              double *T, unsigned int N, double *PDDs) {
  const double h_days = dt_series / seconds_per_day;

  for (unsigned int k = 0; k < N; ++k) {
    PDDs[k] = h_days * CalovGreveIntegrand(S[k], T[k] - pdd_threshold_temp);
  }
}


//! \brief Extract snow accumulation from mixed (snow and rain)
//! precipitation using the temperature time-series.
/** Uses the temperature time-series to determine whether the
 * precipitation is snow or rain. Rain is removed entirely from the
 * surface mass balance, and will not be included in the computed
 * runoff, which is meltwater runoff. There is an allowed linear
 * transition for Tmin below which all precipitation is interpreted as
 * snow, and Tmax above which all precipitation is rain (see, e.g.
 * [\ref Hock2005b]).
 *
 * Sets P[i] to the *solid* (snow) accumulation *rate*.
 *
 * @param[in,out] P precipitation rate (array of length N)
 * @param[in] T air temperature (array of length N)
 * @param[in] N array length
 */
void PDDMassBalance::get_snow_accumulation(double *P, double *T,
                                           unsigned int N) {

  // Following \ref Hock2005b we employ a linear transition from Tmin to Tmax
  for (unsigned int i = 0; i < N; i++) {
    // do not allow negative precipitation
    if (P[i] < 0.0) {
      P[i] = 0.0;
      continue;
    }

    if (precip_as_snow || T[i] <= Tmin) { // T <= Tmin, all precip is snow
      // no change
    } else if (T[i] < Tmax) { // linear transition from Tmin to Tmax
      P[i] *= (Tmax - T[i]) / (Tmax - Tmin);
    } else { // T >= Tmax, all precip is rain -- ignore it
      P[i] = 0.0;
    }
  }

}


//! \brief Compute the surface mass balance at a location from the number of positive
//! degree days and the accumulation amount in a time interval.
/*!
 * This is a PDD scheme. The input parameter `ddf.snow` is a rate of
 * melting per positive degree day for snow.
 *
 * - a fraction of the melted snow and ice refreezes, conceptualized
 *   as superimposed ice, and this is controlled by parameter \c
 *   ddf.refreezeFrac
 *
 * - the excess number of PDDs is used to melt both the ice that came
 *   from refreeze and then any ice which is already present.
 *
 * Ice melts at a constant rate per positive degree day, controlled by
 * parameter `ddf.ice`.
 *
 * The scheme here came from EISMINT-Greenland [\ref RitzEISMINT], but
 * is influenced by R. Hock (personal communication).
 */
void PDDMassBalance::step(const DegreeDayFactors &ddf,
                          double PDDs,
                          double accumulation,
                          double &snow_depth,
                          double &cumulative_melt,
                          double &cumulative_runoff,
                          double &cumulative_smb) {

  double
    max_snow_melted = PDDs * ddf.snow,
    snow_melted = 0.0, excess_pdds = 0.0;

  snow_depth += accumulation;

  if (PDDs <= 0.0) {       // The "no melt" case.
    snow_melted = 0.0;
    excess_pdds = 0.0;
  } else if (max_snow_melted <= snow_depth) {
    // Some of the snow melted and some is left; in any case, all of
    // the energy available for melt, namely all of the positive
    // degree days (PDDs) were used up in melting snow.

    snow_melted = max_snow_melted;
    excess_pdds = 0.0;
  } else {
    // All (snow_depth meters) of snow melted. Excess_pddsum is the
    // positive degree days available to melt ice.
    snow_melted = snow_depth;
    excess_pdds = PDDs - (snow_melted / ddf.snow); // units: K day
  }

  double
    ice_melted              = excess_pdds * ddf.ice,
    melt                    = snow_melted + ice_melted,
    ice_created_by_refreeze = 0.0, runoff = 0.0;

  if (refreeze_ice_melt) {
    ice_created_by_refreeze = melt * ddf.refreezeFrac;
  } else {
    ice_created_by_refreeze = snow_melted * ddf.refreezeFrac;
  }

  runoff = melt - ice_created_by_refreeze;

  snow_depth -= snow_melted;
  if (snow_depth < 0.0) {
    snow_depth = 0.0;
  }

  cumulative_melt   += melt;
  cumulative_runoff += runoff;
  cumulative_smb    += accumulation - runoff;
}


/*!
Initializes the random number generator (RNG).  The RNG is GSL's recommended default,
which seems to be "mt19937" and is DIEHARD (whatever that means ...). Seed with
wall clock time in seconds in non-repeatable case, and with 0 in repeatable case.
 */
PDDrandMassBalance::PDDrandMassBalance(const Config& myconfig,
                                       bool repeatable)
  : PDDMassBalance(myconfig) {
  pddRandGen = gsl_rng_alloc(gsl_rng_default);  // so pddRandGen != NULL now
  gsl_rng_set(pddRandGen, repeatable ? 0 : time(0));
}


PDDrandMassBalance::~PDDrandMassBalance() {
  if (pddRandGen != NULL) {
    gsl_rng_free(pddRandGen);
    pddRandGen = NULL;
  }
}


/*! We need to compute simulated random temperature each actual \e
  day, or at least as close as we can reasonably get. Output `N` is
  number of days or number of days plus one.

  Thus this method ignores
  `config.get("pdd_max_evals_per_year")`, which is
  used in the base class PDDMassBalance.

  Implementation of get_PDDs() requires returned N >= 2, so we
  guarantee that.
 */
unsigned int PDDrandMassBalance::get_timeseries_length(double dt) {
  return std::max(static_cast<size_t>(ceil(dt / seconds_per_day)), (size_t)2);
}

/** 
 * Computes
 * \f[
 * \text{PDD} = \sum_{i=0}^{N-1} h_{\text{days}} \cdot \text{max}(T_i-T_{\text{threshold}}, 0).
 * \f]
 * 
 * @param S \f$\sigma\f$ (standard deviation for daily temperature excursions)
 * @param dt_series time-series step, in seconds
 * @param T air temperature
 * @param N number of points in the temperature time-series, each corresponds to a sub-interval
 * @param PDDs pointer to a pre-allocated array of length N
 */
void PDDrandMassBalance::get_PDDs(double *S, double dt_series,
                                  double *T, unsigned int N, double *PDDs) {
  const double h_days = dt_series / seconds_per_day;

  for (unsigned int k = 0; k < N; ++k) {
    // average temperature in k-th interval
    double T_k = T[k] + gsl_ran_gaussian(pddRandGen, S[k]); // add random: N(0,sigma)

    if (T_k > pdd_threshold_temp) {
      PDDs[k] = h_days * (T_k - pdd_threshold_temp);
    }
  }
}


FaustoGrevePDDObject::FaustoGrevePDDObject(const IceGrid &g)
  : grid(g), config(g.config) {

  beta_ice_w  = config.get("pdd_fausto_beta_ice_w");
  beta_snow_w = config.get("pdd_fausto_beta_snow_w");

  T_c         = config.get("pdd_fausto_T_c");
  T_w         = config.get("pdd_fausto_T_w");
  beta_ice_c  = config.get("pdd_fausto_beta_ice_c");
  beta_snow_c = config.get("pdd_fausto_beta_snow_c");

  fresh_water_density        = config.get("fresh_water_density");
  ice_density                = config.get("ice_density");
  pdd_fausto_latitude_beta_w = config.get("pdd_fausto_latitude_beta_w");

  temp_mj.create(grid, "temp_mj_faustogreve", WITHOUT_GHOSTS);
  temp_mj.set_attrs("internal",
                    "mean July air temp from Fausto et al (2009) parameterization",
                    "K", "");
}


PetscErrorCode FaustoGrevePDDObject::setDegreeDayFactors(int i, int j,
                                                         double /* usurf */,
                                                         double lat, double /* lon */,
                                                         DegreeDayFactors &ddf) {

  IceModelVec::AccessList list(temp_mj);
  const double T_mj = temp_mj(i,j);

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
      const double
         lam_i = pow((T_w - T_mj) / (T_w - T_c) , 3.0),
         lam_s = (T_mj - T_c) / (T_w - T_c);
      ddf.ice  = beta_ice_w + (beta_ice_c - beta_ice_w) * lam_i;
      ddf.snow = beta_snow_w + (beta_snow_c - beta_snow_w) * lam_s;
    }
  }

  // degree-day factors in \ref Faustoetal2009 are water-equivalent
  //   thickness per degree day; ice-equivalent thickness melted per degree
  //   day is slightly larger; for example, iwfactor = 1000/910
  const double iwfactor = fresh_water_density / ice_density;
  ddf.snow *= iwfactor;
  ddf.ice  *= iwfactor;
  return 0;
}


//! Updates mean July near-surface air temperature.
/*!
Unfortunately this duplicates code in PA_SeaRISE_Greenland::update();
 */
PetscErrorCode FaustoGrevePDDObject::update_temp_mj(IceModelVec2S *surfelev,
                                                    IceModelVec2S *lat, IceModelVec2S *lon) {
  const double
    d_mj     = config.get("snow_temp_fausto_d_mj"),      // K
    gamma_mj = config.get("snow_temp_fausto_gamma_mj"),  // K m-1
    c_mj     = config.get("snow_temp_fausto_c_mj"),      // K (degN)-1
    kappa_mj = config.get("snow_temp_fausto_kappa_mj");  // K (degW)-1

  IceModelVec2S &h = *surfelev, &lat_degN = *lat, &lon_degE = *lon;

  IceModelVec::AccessList list;
  list.add(h);
  list.add(lat_degN);
  list.add(lon_degE);
  list.add(temp_mj);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    temp_mj(i,j) = d_mj + gamma_mj * h(i,j) + c_mj * lat_degN(i,j) + kappa_mj * (-lon_degE(i,j));
  }

  return 0;
}

} // end of namespace pism
