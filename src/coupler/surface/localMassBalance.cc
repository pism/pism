// Copyright (C) 2009, 2010, 2011, 2013, 2014, 2015, 2016, 2017, 2018 Ed Bueler and Constantine Khroulev and Andy Aschwanden
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

#include <cassert>
#include <ctime>  // for time(), used to initialize random number gen
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>       // M_PI
#include <cmath>                // for erfc() in CalovGreveIntegrand()
#include <algorithm>

#include "pism/util/pism_utilities.hh"
#include "pism/util/ConfigInterface.hh"
#include "localMassBalance.hh"
#include "pism/util/IceGrid.hh"

namespace pism {
namespace surface {

LocalMassBalance::Changes::Changes() {
  firn_depth    = 0.0;
  snow_depth    = 0.0;
  melt          = 0.0;
  runoff        = 0.0;
  smb           = 0.0;
}

LocalMassBalance::LocalMassBalance(Config::ConstPtr myconfig, units::System::Ptr system)
  : m_config(myconfig), m_unit_system(system),
    m_seconds_per_day(86400) {
  // empty
}

LocalMassBalance::~LocalMassBalance() {
  // empty
}

std::string LocalMassBalance::method() const {
  return m_method;
}

PDDMassBalance::PDDMassBalance(Config::ConstPtr config, units::System::Ptr system)
  : LocalMassBalance(config, system) {
  precip_as_snow     = m_config->get_boolean("surface.pdd.interpret_precip_as_snow");
  Tmin               = m_config->get_double("surface.pdd.air_temp_all_precip_as_snow");
  Tmax               = m_config->get_double("surface.pdd.air_temp_all_precip_as_rain");
  pdd_threshold_temp = m_config->get_double("surface.pdd.positive_threshold_temp");
  refreeze_ice_melt  = m_config->get_boolean("surface.pdd.refreeze_ice_melt");

  m_method = "an expectation integral";
}


/*! \brief Compute the number of points for temperature and
    precipitation time-series.
 */
unsigned int PDDMassBalance::get_timeseries_length(double dt) {
  const unsigned int    NperYear = static_cast<unsigned int>(m_config->get_double("surface.pdd.max_evals_per_year"));
  const double dt_years = units::convert(m_unit_system, dt, "seconds", "years");

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

This integral is used for the expected number of positive degree days. The user can choose
\f$\sigma\f$ by option `-pdd_std_dev`. Note that the integral is over a time interval of
length `dt` instead of a whole year as stated in \ref CalovGreve05 . If `sigma` is zero,
return the positive part of `TacC`.
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
void PDDMassBalance::get_PDDs(double dt_series,
                              const std::vector<double> &S,
                              const std::vector<double> &T,
                              std::vector<double> &PDDs) {
  assert(S.size() == T.size() and T.size() == PDDs.size());
  assert(dt_series > 0.0);

  const double h_days = dt_series / m_seconds_per_day;
  const size_t N = S.size();

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
void PDDMassBalance::get_snow_accumulation(const std::vector<double> &T,
                                           std::vector<double> &P) {

  assert(T.size() == P.size());
  const size_t N = T.size();

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
 * `accumulation` has units "meter / second".
 *
 * - a fraction of the melted snow and ice refreezes, conceptualized
 *   as superimposed ice, and this is controlled by parameter \c
 *   ddf.refreeze_fraction
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
PDDMassBalance::Changes PDDMassBalance::step(const DegreeDayFactors &ddf,
                                             double PDDs,
                                             double thickness,
                                             double old_firn_depth,
                                             double old_snow_depth,
                                             double accumulation) {
  double
    firn_depth      = old_firn_depth,
    snow_depth      = old_snow_depth,
    max_snow_melted = PDDs * ddf.snow,
    firn_melted     = 0.0,
    snow_melted     = 0.0,
    excess_pdds     = 0.0;

  assert(thickness >= 0);

  // snow depth cannot exceed total thickness
  snow_depth = std::min(snow_depth, thickness);

  assert(snow_depth >= 0);

  // firn depth cannot exceed thickness - snow_depth
  firn_depth = std::min(firn_depth, thickness - snow_depth);

  assert(firn_depth >= 0);

  double ice_thickness = thickness - snow_depth - firn_depth;

  assert(ice_thickness >= 0);

  snow_depth += accumulation;

  if (PDDs <= 0.0) {            // The "no melt" case.
    snow_melted = 0.0;
    firn_melted = 0.0,
    excess_pdds = 0.0;
  } else if (max_snow_melted <= snow_depth) {
    // Some of the snow melted and some is left; in any case, all of
    // the energy available for melt, namely all of the positive
    // degree days (PDDs) were used up in melting snow.
    snow_melted = max_snow_melted;
    firn_melted = 0.0;
    excess_pdds = 0.0;
  } else if (max_snow_melted <= firn_depth + snow_depth) {
    // All of the snow is melted but some firn is left; in any case, all of
    // the energy available for melt, namely all of the positive
    // degree days (PDDs) were used up in melting snow and firn.
    snow_melted = snow_depth;
    firn_melted = max_snow_melted - snow_melted;
    excess_pdds = 0.0;
  } else {
    // All (firn_depth and snow_depth meters) of snow and firn melted. Excess_pdds is the
    // positive degree days available to melt ice.
    firn_melted = firn_depth;
    snow_melted = snow_depth;
    excess_pdds = PDDs - ((firn_melted + snow_melted) / ddf.snow); // units: K day
  }

  double
    ice_melted              = std::min(excess_pdds * ddf.ice, ice_thickness),
    melt                    = snow_melted + firn_melted + ice_melted,
    ice_created_by_refreeze = 0.0;

  if (refreeze_ice_melt) {
    ice_created_by_refreeze = melt * ddf.refreeze_fraction;
  } else {
    // Should this only be snow melted?
    ice_created_by_refreeze = (firn_melted + snow_melted) * ddf.refreeze_fraction;
  }

  const double runoff = melt - ice_created_by_refreeze;

  snow_depth = std::max(snow_depth - snow_melted, 0.0);
  firn_depth = std::max(firn_depth - firn_melted, 0.0);

  // FIXME: need to add snow that hasn't melted, is this correct?
  // firn_depth += (snow_depth - snow_melted);
  // Turn firn into ice at X times accumulation
  // firn_depth -= accumulation *  m_config->get_double("surface.pdd.firn_compaction_to_accumulation_ratio");

  const double smb = accumulation - runoff;

  Changes result;
  // Ensure that we never generate negative ice thicknesses. As far as I can tell the code
  // above guarantees that thickness + smb >= *in exact arithmetic*. The check below
  // should make sure that we don't get bitten by rounding errors.
  result.smb        = thickness + smb >= 0 ? smb : -thickness;
  result.firn_depth = firn_depth - old_firn_depth;
  result.snow_depth = snow_depth - old_snow_depth;
  result.melt       = melt;
  result.runoff     = runoff;

  assert(thickness + result.smb >= 0);

  return result;
}


/*!
Initializes the random number generator (RNG).  The RNG is GSL's recommended default,
which seems to be "mt19937" and is DIEHARD (whatever that means ...). Seed with
wall clock time in seconds in non-repeatable case, and with 0 in repeatable case.
 */
PDDrandMassBalance::PDDrandMassBalance(Config::ConstPtr config, units::System::Ptr system,
                                       Kind kind)
  : PDDMassBalance(config, system) {
  pddRandGen = gsl_rng_alloc(gsl_rng_default);  // so pddRandGen != NULL now
  gsl_rng_set(pddRandGen, kind == REPEATABLE ? 0 : time(0));

  m_method = (kind == NOT_REPEATABLE
              ? "simulation of a random process"
              : "repeatable simulation of a random process");
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
  `config.get_double("surface.pdd.max_evals_per_year")`, which is
  used in the base class PDDMassBalance.

  Implementation of get_PDDs() requires returned N >= 2, so we
  guarantee that.
 */
unsigned int PDDrandMassBalance::get_timeseries_length(double dt) {
  return std::max(static_cast<size_t>(ceil(dt / m_seconds_per_day)), (size_t)2);
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
void PDDrandMassBalance::get_PDDs(double dt_series,
                                  const std::vector<double> &S,
                                  const std::vector<double> &T,
                                  std::vector<double> &PDDs) {
  assert(S.size() == T.size() and T.size() == PDDs.size());
  assert(dt_series > 0.0);

  const double h_days = dt_series / m_seconds_per_day;
  const size_t N = S.size();

  for (unsigned int k = 0; k < N; ++k) {
    // average temperature in k-th interval
    double T_k = T[k] + gsl_ran_gaussian(pddRandGen, S[k]); // add random: N(0,sigma)

    if (T_k > pdd_threshold_temp) {
      PDDs[k] = h_days * (T_k - pdd_threshold_temp);
    }
  }
}


FaustoGrevePDDObject::FaustoGrevePDDObject(IceGrid::ConstPtr g)
  : m_grid(g), m_config(g->ctx()->config()) {

  m_beta_ice_w  = m_config->get_double("surface.pdd.fausto.beta_ice_w");
  m_beta_snow_w = m_config->get_double("surface.pdd.fausto.beta_snow_w");

  m_T_c         = m_config->get_double("surface.pdd.fausto.T_c");
  m_T_w         = m_config->get_double("surface.pdd.fausto.T_w");
  m_beta_ice_c  = m_config->get_double("surface.pdd.fausto.beta_ice_c");
  m_beta_snow_c = m_config->get_double("surface.pdd.fausto.beta_snow_c");

  m_fresh_water_density        = m_config->get_double("constants.fresh_water.density");
  m_ice_density                = m_config->get_double("constants.ice.density");
  m_pdd_fausto_latitude_beta_w = m_config->get_double("surface.pdd.fausto.latitude_beta_w");
  m_refreeze_fraction = m_config->get_double("surface.pdd.refreeze");


  m_temp_mj.create(m_grid, "temp_mj_faustogreve", WITHOUT_GHOSTS);
  m_temp_mj.set_attrs("internal",
                    "mean July air temp from Fausto et al (2009) parameterization",
                    "K", "");
}

FaustoGrevePDDObject::~FaustoGrevePDDObject() {
  // empty
}

LocalMassBalance::DegreeDayFactors FaustoGrevePDDObject::degree_day_factors(int i, int j,
                                                                            double latitude) {

  LocalMassBalance::DegreeDayFactors ddf;
  ddf.refreeze_fraction = m_refreeze_fraction;

  IceModelVec::AccessList list(m_temp_mj);
  const double T_mj = m_temp_mj(i,j);

  if (latitude < m_pdd_fausto_latitude_beta_w) { // case latitude < 72 deg N
    ddf.ice  = m_beta_ice_w;
    ddf.snow = m_beta_snow_w;
  } else { // case > 72 deg N
    if (T_mj >= m_T_w) {
      ddf.ice  = m_beta_ice_w;
      ddf.snow = m_beta_snow_w;
    } else if (T_mj <= m_T_c) {
      ddf.ice  = m_beta_ice_c;
      ddf.snow = m_beta_snow_c;
    } else { // middle case   T_c < T_mj < T_w
      const double
        lam_i = pow((m_T_w - T_mj) / (m_T_w - m_T_c) , 3.0),
        lam_s = (T_mj - m_T_c) / (m_T_w - m_T_c);
      ddf.ice  = m_beta_ice_w + (m_beta_ice_c - m_beta_ice_w) * lam_i;
      ddf.snow = m_beta_snow_w + (m_beta_snow_c - m_beta_snow_w) * lam_s;
    }
  }

  // degree-day factors in \ref Faustoetal2009 are water-equivalent
  //   thickness per degree day; ice-equivalent thickness melted per degree
  //   day is slightly larger; for example, iwfactor = 1000/910
  const double iwfactor = m_fresh_water_density / m_ice_density;
  ddf.snow *= iwfactor;
  ddf.ice  *= iwfactor;

  return ddf;
}


//! Updates mean July near-surface air temperature.
/*!
Unfortunately this duplicates code in SeaRISEGreenland::update();
 */
void FaustoGrevePDDObject::update_temp_mj(const IceModelVec2S &surfelev,
                                          const IceModelVec2S &lat,
                                          const IceModelVec2S &lon) {
  const double
    d_mj     = m_config->get_double("atmosphere.fausto_air_temp.d_mj"),      // K
    gamma_mj = m_config->get_double("atmosphere.fausto_air_temp.gamma_mj"),  // K m-1
    c_mj     = m_config->get_double("atmosphere.fausto_air_temp.c_mj"),      // K (degN)-1
    kappa_mj = m_config->get_double("atmosphere.fausto_air_temp.kappa_mj");  // K (degW)-1

  const IceModelVec2S
    &h        = surfelev,
    &lat_degN = lat,
    &lon_degE = lon;

  IceModelVec::AccessList list{&h, &lat_degN, &lon_degE, &m_temp_mj};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    m_temp_mj(i,j) = d_mj + gamma_mj * h(i,j) + c_mj * lat_degN(i,j) + kappa_mj * (-lon_degE(i,j));
  }
}

} // end of namespace surface
} // end of namespace pism
