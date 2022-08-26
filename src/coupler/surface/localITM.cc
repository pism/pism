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
#include <cmath>          // M_PI and erfc() in CalovGreveIntegrand()

#include "localITM.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace surface {

ITMMassBalance::Changes::Changes() {
  snow_depth = 0.0;
  melt       = 0.0;
  runoff     = 0.0;
  smb        = 0.0;
}

ITMMassBalance::Melt::Melt() {
  T_melt   = 0.0;
  I_melt   = 0.0;
  c_melt   = 0.0;
  ITM_melt = 0.0;
}

ITMMassBalance::ITMMassBalance(const Config &config, units::System::Ptr system) {
  m_precip_as_snow     = config.get_flag("surface.itm.interpret_precip_as_snow");
  m_Tmin               = config.get_number("surface.itm.air_temp_all_precip_as_snow");
  m_Tmax               = config.get_number("surface.itm.air_temp_all_precip_as_rain");
  m_refreeze_ice_melt  = config.get_flag("surface.itm.refreeze_ice_melt");
  m_pdd_threshold_temp = config.get_number("surface.itm.positive_threshold_temp");

  m_year_length = units::convert(system, 1.0, "years", "seconds");
  m_n_per_year = static_cast<unsigned int>(config.get_number("surface.itm.max_evals_per_year"));

  m_water_density = config.get_number("constants.fresh_water.density");
  m_ice_density   = config.get_number("constants.ice.density");
  m_albedo_snow   = config.get_number("surface.itm.albedo_snow");
  m_albedo_ice    = config.get_number("surface.itm.albedo_ice");
  m_albedo_land   = config.get_number("surface.itm.albedo_land"); //0.2
  m_albedo_ocean  = config.get_number("surface.itm.albedo_ocean"); // 0.1;
  m_albedo_slope  = config.get_number("surface.itm.albedo_slope"); //-790;

  m_tau_a_slope     = config.get_number("surface.itm.tau_a_slope");
  m_tau_a_intercept = config.get_number("surface.itm.tau_a_intercept");

  m_itm_c      = config.get_number("surface.itm.itm_c");
  m_itm_lambda = config.get_number("surface.itm.itm_lambda");
  m_bm_temp    = config.get_number("surface.itm.background_melting_temp");

  m_L = config.get_number("constants.fresh_water.latent_heat_of_fusion");
  m_solar_constant = config.get_number("surface.itm.solar_constant");

  m_phi = config.get_number("surface.itm.phi") * M_PI / 180.0;

  m_Tmin_refreeze = config.get_number("surface.itm.air_temp_all_refreeze");
  m_Tmax_refreeze = config.get_number("surface.itm.air_temp_no_refreeze");
}


/*! \brief Compute the number of points for temperature and
    precipitation time-series.
 */
unsigned int ITMMassBalance::timeseries_length(double dt) {
  double dt_years = dt / m_year_length;

  return std::max(1U, static_cast<unsigned int>(ceil(m_n_per_year * dt_years)));
}


double ITMMassBalance::CalovGreveIntegrand(double sigma, double TacC) {

  if (sigma == 0) {
    return std::max(TacC, 0.0);
  } else {
    const double Z = TacC / (sqrt(2.0) * sigma);
    return (sigma / sqrt(2.0 * M_PI)) * exp(-Z * Z) + (TacC / 2.0) * erfc(-Z);
  }
}


double ITMMassBalance::albedo(double melt,
                              MaskValue cell_type,
                              double dtseries) {
  // melt has a unit of meters ice equivalent
  // dtseries has a unit of seconds
  if (cell_type == MASK_ICE_FREE_OCEAN) {
    return m_albedo_ocean;
  }

  if (cell_type == MASK_ICE_FREE_BEDROCK) {
    return m_albedo_land;
  }

  double result = m_albedo_snow + m_albedo_slope * melt * m_ice_density / (dtseries);
  return std::max(result, m_albedo_ice);
}


double ITMMassBalance::get_tau_a(double surface_elevation) {
  // transmissivity of the atmosphere, linear fit
  return m_tau_a_intercept + m_tau_a_slope * surface_elevation;
}


double ITMMassBalance::get_h_phi(double phi, double lat, double delta) {
  // calculate the hour angle at which the sun reaches phi (for melting period during the day)
  double input_h_phi         = (sin(phi) - sin(lat) * sin(delta)) / (cos(lat) * cos(delta));
  double input_h_phi_clipped = std::max(-1., std::min(input_h_phi, 1.));
  return acos(input_h_phi_clipped);
}


double ITMMassBalance::get_q_insol(double distance2, double h_phi,
                                   double lat, double delta) {
  if (h_phi == 0) {
    return 0.;
  } else {
    return m_solar_constant * distance2 * (h_phi * sin(lat) * sin(delta) + cos(lat) * cos(delta) * sin(h_phi)) / h_phi;
  }
}

double ITMMassBalance::get_TOA_insol(double distance2, double h0,
                                     double lat, double delta) {
  if (h0 == 0) {
    return 0.;
  } else {
    return m_solar_constant * distance2 * (h0 * sin(lat) * sin(delta) + cos(lat) * cos(delta) * sin(h0)) / M_PI;
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
ITMMassBalance::Melt ITMMassBalance::calculate_melt(double dt_series, double S, double T,
                                                    double surface_elevation, double delta,
                                                    double distance2, double lat,
                                                    double albedo) {
  assert(dt_series > 0.0);

  Melt result;

  const double tau_a = get_tau_a(surface_elevation);
  // if background melting is true, use effective pdd temperatures and do not allow melting below background meltin temp.

  double h_phi            = get_h_phi(m_phi, lat, delta);
  double h0               = get_h_phi(0, lat, delta);
  double quotient_delta_t = h_phi / M_PI;

  result.transmissivity = tau_a;

  double q_insol = get_q_insol(distance2, h_phi, lat, delta);

  result.TOA_insol = get_TOA_insol(distance2, h0, lat, delta);
  result.q_insol   = q_insol;

  double Teff = CalovGreveIntegrand(S, T - m_pdd_threshold_temp);
  if (Teff < 1.e-4) {
    Teff = 0;
  }

  result.T_melt = quotient_delta_t * dt_series / (m_water_density * m_L) * m_itm_lambda * (Teff);

  if (T < m_bm_temp) {
    result.ITM_melt = 0.;
  } else {
    result.ITM_melt =
        quotient_delta_t * dt_series / (m_water_density * m_L) * (tau_a * (1. - albedo) * q_insol + m_itm_c + m_itm_lambda * (Teff));
  }

  result.I_melt = dt_series / (m_water_density * m_L) * (tau_a * (1. - albedo) * q_insol) * quotient_delta_t;
  result.c_melt = dt_series / (m_water_density * m_L) * m_itm_c * quotient_delta_t;

  return result;
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
void ITMMassBalance::get_snow_accumulation(const std::vector<double> &T, std::vector<double> &P) {


  assert(T.size() == P.size());
  const size_t N = T.size();

  // Following \ref Hock2005b we employ a linear transition from Tmin to Tmax
  for (unsigned int i = 0; i < N; i++) {
    // do not allow negative precipitation
    if (P[i] < 0.0) {
      P[i] = 0.0;
      continue;
    }

    if (m_precip_as_snow or T[i] <= m_Tmin) { // T <= Tmin, all precip is snow
      // no change
    } else if (T[i] < m_Tmax) { // linear transition from Tmin to Tmax
      P[i] *= (m_Tmax - T[i]) / (m_Tmax - m_Tmin);
    } else { // T >= Tmax, all precip is rain -- ignore it
      P[i] = 0.0;
    }
  }
}


double ITMMassBalance::refreeze_fraction(double T) {

  if (T <= m_Tmin_refreeze) {
    return 1.0;
  }

  if ((m_Tmin_refreeze < T) and (T <= m_Tmax_refreeze)) {
    return 1.0 / (m_Tmin_refreeze - m_Tmax_refreeze) * T + m_Tmax_refreeze / (m_Tmax_refreeze - m_Tmin_refreeze);
  }

  return 0.0;
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
ITMMassBalance::Changes ITMMassBalance::step(double refreeze_fraction, double thickness, double ITM_melt,
                                             double old_firn_depth, double old_snow_depth, double accumulation) {


  Changes result;


  double
    firn_depth      = old_firn_depth,
    snow_depth      = old_snow_depth,
    max_snow_melted = ITM_melt,
    firn_melted     = 0.0,
    snow_melted     = 0.0;

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
    firn_melted = 0.0;
  } else if (max_snow_melted <= snow_depth) {
    // Some of the snow melted and some is left; in any case, all of
    // the energy available for melt was used up in melting snow.
    snow_melted = max_snow_melted;
    firn_melted = 0.0;
  } else if (max_snow_melted <= firn_depth + snow_depth) {
    // All of the snow is melted but some firn is left; in any case, all of
    // the energy available for melt was used up in melting snow.
    snow_melted = snow_depth;
    firn_melted = max_snow_melted - snow_melted;
  } else {
    // All (firn and snow_depth meters) of snow melted. Excess_melt is
    // available to melt ice.
    firn_melted = firn_depth;
    snow_melted = snow_depth;
  }

  double melt                 = ITM_melt,
      ice_created_by_refreeze = 0.0;


  if (m_refreeze_ice_melt) {
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
  // firn_depth -= accumulation *  config.get_number("surface.pdd.firn_compaction_to_accumulation_ratio");


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
