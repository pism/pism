// Copyright (C) 2009, 2010, 2011, 2013, 2014, 2015, 2016, 2017, 2022 Ed Bueler and Constantine Khroulev and Andy Aschwanden
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

#include <algorithm>
#include <cassert>
#include <cmath>          // for erfc() in CalovGreveIntegrand()
#include <ctime>          // for time(), used to initialize random number gen
#include <gsl/gsl_math.h> // M_PI
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iostream>

#include "localITM.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace surface {

LocalMassBalanceITM::Changes::Changes() {
  snow_depth = 0.0;
  melt       = 0.0;
  runoff     = 0.0;
  smb        = 0.0;
}

LocalMassBalanceITM::Melt::Melt() {
  T_melt   = 0.0;
  I_melt   = 0.0;
  c_melt   = 0.0;
  ITM_melt = 0.0;
}

LocalMassBalanceITM::LocalMassBalanceITM(Config::ConstPtr myconfig, units::System::Ptr system)
    : m_config(myconfig), m_unit_system(system), m_seconds_per_day(86400) {
  // empty
}

LocalMassBalanceITM::~LocalMassBalanceITM() {
  // empty
}

std::string LocalMassBalanceITM::method() const {
  return m_method;
}

ITMMassBalance::ITMMassBalance(Config::ConstPtr config, units::System::Ptr system)
    : LocalMassBalanceITM(config, system) {
  precip_as_snow     = m_config->get_flag("surface.itm.interpret_precip_as_snow");
  Tmin               = m_config->get_number("surface.itm.air_temp_all_precip_as_snow");
  Tmax               = m_config->get_number("surface.itm.air_temp_all_precip_as_rain");
  refreeze_ice_melt  = m_config->get_flag("surface.itm.refreeze_ice_melt");
  pdd_threshold_temp = m_config->get_number("surface.itm.positive_threshold_temp");


  m_method = "insolation temperature melt";
}


/*! \brief Compute the number of points for temperature and
    precipitation time-series.
 */
unsigned int ITMMassBalance::get_timeseries_length(double dt) {
  const unsigned int NperYear = static_cast<unsigned int>(m_config->get_number("surface.itm.max_evals_per_year"));
  const double dt_years       = units::convert(m_unit_system, dt, "seconds", "years");

  return std::max(1U, static_cast<unsigned int>(ceil(NperYear * dt_years)));
}


double ITMMassBalance::CalovGreveIntegrand(double sigma, double TacC) {

  if (sigma == 0) {
    return std::max(TacC, 0.0);
  } else {
    const double Z = TacC / (sqrt(2.0) * sigma);
    return (sigma / sqrt(2.0 * M_PI)) * exp(-Z * Z) + (TacC / 2.0) * erfc(-Z);
  }
}


double ITMMassBalance::get_albedo_melt(double melt, int mask_value, double dtseries) {
  const double ice_density  = m_config->get_number("constants.ice.density");
  double albedo             = m_config->get_number("surface.itm.albedo_snow");
  const double albedo_land  = m_config->get_number("surface.itm.albedo_land");  //0.2
  const double albedo_ocean = m_config->get_number("surface.itm.albedo_ocean"); // 0.1;
  // melt has a unit of meters ice equivalent
  // dtseries has a unit of seconds
  double intersection = m_config->get_number("surface.itm.albedo_snow");  //0.82;
  double slope        = m_config->get_number("surface.itm.albedo_slope"); //-790;
  double albedo_ice   = m_config->get_number("surface.itm.albedo_ice");

  if (mask_value == 4) { // mask value for ice free ocean
    albedo = albedo_ocean;
  }
  if (mask_value == 0) { // mask value for bedrock
    albedo = albedo_land;
  } else {
    albedo = intersection + slope * melt * ice_density / (dtseries);
    if (albedo < albedo_ice) {
      albedo = albedo_ice;
    }
  }

  return albedo;
}


double ITMMassBalance::get_tau_a(double surface_elevation) {
  const double tau_a_slope     = m_config->get_number("surface.itm.tau_a_slope");
  const double tau_a_intercept = m_config->get_number("surface.itm.tau_a_intercept");
  return tau_a_intercept + tau_a_slope * surface_elevation; // transmissivity of the atmosphere, linear fit
}


double ITMMassBalance::get_h_phi(const double &phi, const double &lat, const double &delta) {
  // calculate the hour angle at which the sun reaches phi (for melting period during the day)
  double input_h_phi         = (sin(phi) - sin(lat) * sin(delta)) / (cos(lat) * cos(delta));
  double input_h_phi_clipped = std::max(-1., std::min(input_h_phi, 1.));
  return acos(input_h_phi_clipped);
}


double ITMMassBalance::get_q_insol(const double &solar_constant, const double &distance2, const double &h_phi,
                                   const double &lat, const double &delta) {
  if (h_phi == 0) {
    return 0.;
  } else {
    return solar_constant * distance2 * (h_phi * sin(lat) * sin(delta) + cos(lat) * cos(delta) * sin(h_phi)) / h_phi;
  }
}

double ITMMassBalance::get_TOA_insol(const double &solar_constant, const double &distance2, const double &h0,
                                     const double &lat, const double &delta) {
  if (h0 == 0) {
    return 0.;
  } else {
    return solar_constant * distance2 * (h0 * sin(lat) * sin(delta) + cos(lat) * cos(delta) * sin(h0)) / M_PI;
  }
}


//!
/* compute diurnal melt scheme  by equation (6) by Uta Krebs-Kanzow et al., The Cryosphere, 2018
 * @param dt_series length of the step for the time-series
 * @param T air temperature at time [k]
 * @param insolation at time [k]
 * @param surface_elevation
 * @param albedo which was should be figured by get_albedo (?)
 * @param[out] melt pointer to a pre-allocated array with N-1 elements
 * output in mm water equivalent
 */


ITMMassBalance::Melt ITMMassBalance::calculate_ETIM_melt(double dt_series, const double &S, const double &T,
                                                         const double &surface_elevation, const double &delta,
                                                         const double &distance2, const double &lat,
                                                         const double &albedo) {
  assert(dt_series > 0.0);

  Melt ETIM_melt;

  const double rho_w = m_config->get_number("constants.fresh_water.density");             // mass density of water
  const double L_m = m_config->get_number("constants.fresh_water.latent_heat_of_fusion"); // latent heat of ice melting
  const double tau_a      = get_tau_a(surface_elevation);
  const double itm_c      = m_config->get_number("surface.itm.itm_c");
  const double itm_lambda = m_config->get_number("surface.itm.itm_lambda");
  const double bm_temp    = m_config->get_number("surface.itm.background_melting_temp");
  // if background melting is true, use effective pdd temperatures and do not allow melting below background meltin temp.
  const double solar_constant = m_config->get_number("surface.itm.solar_constant");


  const double phi = m_config->get_number("surface.itm.phi") * M_PI / 180.;


  double h_phi = get_h_phi(phi, lat, delta);


  double h0 = get_h_phi(0, lat, delta);


  double quotient_delta_t = h_phi / M_PI;

  ETIM_melt.transmissivity = tau_a;

  double q_insol = get_q_insol(solar_constant, distance2, h_phi, lat, delta);

  ETIM_melt.TOA_insol = get_TOA_insol(solar_constant, distance2, h0, lat, delta);
  ETIM_melt.q_insol   = q_insol;

  assert(dt_series > 0.0);
  double Teff = 0;
  Teff        = CalovGreveIntegrand(S, T - pdd_threshold_temp);
  if (Teff < 1.e-4) {
    Teff = 0;
  }

  ETIM_melt.T_melt = quotient_delta_t * dt_series / (rho_w * L_m) * itm_lambda * (Teff);

  if (T < bm_temp) {
    ETIM_melt.ITM_melt = 0.;
  } else {
    ETIM_melt.ITM_melt =
        quotient_delta_t * dt_series / (rho_w * L_m) * (tau_a * (1. - albedo) * q_insol + itm_c + itm_lambda * (Teff));
  }


  ETIM_melt.I_melt = dt_series / (rho_w * L_m) * (tau_a * (1. - albedo) * q_insol) * quotient_delta_t;
  ETIM_melt.c_melt = dt_series / (rho_w * L_m) * itm_c * quotient_delta_t;


  return ETIM_melt;
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
void ITMMassBalance::get_snow_accumulationITM(const std::vector<double> &T, std::vector<double> &P) {


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


double ITMMassBalance::get_refreeze_fraction(const double &T) {
  double refreeze;
  double Tmin_refreeze = m_config->get_number("surface.itm.air_temp_all_refreeze");
  double Tmax_refreeze = m_config->get_number("surface.itm.air_temp_no_refreeze");
  if (T <= Tmin_refreeze) {
    refreeze = 1.;
  } else if ((Tmin_refreeze < T) and (T <= Tmax_refreeze)) {
    refreeze = 1. / (Tmin_refreeze - Tmax_refreeze) * T + Tmax_refreeze / (Tmax_refreeze - Tmin_refreeze);
  } else {
    refreeze = 0.;
  }
  return refreeze;
}


//! \brief Compute the surface mass balance at a location from the amount of
//! melted snow and the accumulation amount in a time interval.
/*!
 *
 * `accumulation` has units "meter / second".
 *
 * - a fraction of the melted snow and ice refreezes, conceptualized
 *   as superimposed ice, and this is controlled by parameter \c
 *   refreeze_fraction
 *
 * - the excess of 'ITM_melt' is used to melt both the ice that came
 *   from refreeze and then any ice which is already present.
 *
 * The scheme here came from EISMINT-Greenland [\ref RitzEISMINT], but
 * is influenced by R. Hock (personal communication).
 */
ITMMassBalance::Changes ITMMassBalance::step(const double &refreeze_fraction, double thickness, double ITM_melt,
                                             double old_firn_depth, double old_snow_depth, double accumulation) {


  Changes result;


  double firn_depth = old_firn_depth, snow_depth = old_snow_depth, max_snow_melted = ITM_melt, firn_melted = 0.0,
         snow_melted = 0.0, ice_melted = 0.0;

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


  if (ITM_melt <= 0.0) { // The "no melt" case.
    snow_melted = 0.0;
    firn_melted = 0.0, ice_melted = 0.0;
  } else if (max_snow_melted <= snow_depth) {
    // Some of the snow melted and some is left; in any case, all of
    // the energy available for melt was used up in melting snow.
    snow_melted = max_snow_melted;
    firn_melted = 0.0;
    ice_melted  = 0.0;
  } else if (max_snow_melted <= firn_depth + snow_depth) {
    // All of the snow is melted but some firn is left; in any case, all of
    // the energy available for melt was used up in melting snow.
    snow_melted = snow_depth;
    firn_melted = max_snow_melted - snow_melted;
    ice_melted  = 0.0;
  } else {
    // All (firn and snow_depth meters) of snow melted. Excess_melt is
    // available to melt ice.
    firn_melted = firn_depth;
    snow_melted = snow_depth;
    ice_melted  = ITM_melt - (firn_melted + snow_melted); //
  }

  double melt                 = ITM_melt, //snow_melted + firn_melted + ice_melted,
      ice_created_by_refreeze = 0.0;


  if (refreeze_ice_melt) {
    ice_created_by_refreeze = melt * refreeze_fraction;
  } else {
    // Should this only be snow melted?
    ice_created_by_refreeze = (firn_melted + snow_melted) * refreeze_fraction;
  }


  snow_depth = std::max(snow_depth - snow_melted, 0.0);

  firn_depth = std::max(firn_depth - firn_melted, 0.0);
  // FIXME: need to add snow that hasn't melted, is this correct?
  // firn_depth += (snow_depth - snow_melted);
  // Turn firn into ice at X times accumulation
  // firn_depth -= accumulation *  m_config->get_number("surface.pdd.firn_compaction_to_accumulation_ratio");


  const double runoff = melt - ice_created_by_refreeze;
  const double smb    = accumulation - runoff;

  result.firn_depth = firn_depth - old_firn_depth;
  result.snow_depth = snow_depth - old_snow_depth;
  result.melt       = melt;
  result.runoff     = runoff;
  result.smb        = thickness + smb >= 0 ? smb : -thickness;


  assert(thickness + result.smb >= 0);

  return result;
}


} // end of namespace surface
} // end of namespace pism
