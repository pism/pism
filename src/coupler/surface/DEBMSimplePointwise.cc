// Copyright (C) 2009--2022 PISM Authors
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

#include "DEBMSimplePointwise.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Context.hh"
#include "pism/util/Time.hh"

namespace pism {
namespace surface {

DEBMSimplePointwise::Changes::Changes() {
  snow_depth = 0.0;
  melt       = 0.0;
  runoff     = 0.0;
  smb        = 0.0;
}

DEBMSimplePointwise::Melt::Melt() {
  temperature_melt = 0.0;
  insolation_melt  = 0.0;
  background_melt  = 0.0;
  total_melt       = 0.0;
}

DEBMSimplePointwise::DEBMSimplePointwise(const Context &ctx) {

  const Config &config = *ctx.config();
  auto system = ctx.unit_system();

  m_time = ctx.time();

  m_precip_as_snow     = config.get_flag("surface.debm_simple.interpret_precip_as_snow");
  m_Tmin               = config.get_number("surface.debm_simple.air_temp_all_precip_as_snow");
  m_Tmax               = config.get_number("surface.debm_simple.air_temp_all_precip_as_rain");
  m_refreeze_ice_melt  = config.get_flag("surface.debm_simple.refreeze_ice_melt");
  m_refreeze_fraction  = config.get_number("surface.debm_simple.refreeze");
  m_positive_threshold_temperature = config.get_number("surface.debm_simple.positive_threshold_temp");

  m_year_length = units::convert(system, 1.0, "years", "seconds");
  m_n_per_year = static_cast<unsigned int>(config.get_number("surface.debm_simple.max_evals_per_year"));

  m_water_density = config.get_number("constants.fresh_water.density");
  m_ice_density   = config.get_number("constants.ice.density");
  m_albedo_snow   = config.get_number("surface.debm_simple.albedo_snow");
  m_albedo_ice    = config.get_number("surface.debm_simple.albedo_ice");
  m_albedo_land   = config.get_number("surface.debm_simple.albedo_land"); //0.2
  m_albedo_ocean  = config.get_number("surface.debm_simple.albedo_ocean"); // 0.1;
  m_albedo_slope  = config.get_number("surface.debm_simple.albedo_slope"); //-790;

  m_tau_a_slope     = config.get_number("surface.debm_simple.tau_a_slope");
  m_tau_a_intercept = config.get_number("surface.debm_simple.tau_a_intercept");

  m_c1      = config.get_number("surface.debm_simple.c1");
  m_c2      = config.get_number("surface.debm_simple.c2");
  m_bm_temp = config.get_number("surface.debm_simple.background_melting_temp");

  m_L              = config.get_number("constants.fresh_water.latent_heat_of_fusion");
  m_solar_constant = config.get_number("surface.debm_simple.solar_constant");

  m_phi = config.get_number("surface.debm_simple.phi", "radian");

  m_constant_eccentricity         = config.get_number("surface.debm_simple.paleo.eccentricity");
  m_constant_perihelion_longitude = config.get_number("surface.debm_simple.paleo.long_peri", "radian");
  m_constant_obliquity            = config.get_number("surface.debm_simple.paleo.obliquity", "radian");

  m_paleo = config.get_flag("surface.debm_simple.paleo.enabled");

  std::string paleo_file = config.get_string("surface.debm_simple.paleo.file");

  if (not paleo_file.empty()) {
    m_use_paleo_file = true;

    m_eccentricity.reset(
        new ScalarForcing(ctx, "surface.debm_simple.paleo", "eccentricity", "", "", "eccentricity of the earth"));

    m_obliquity.reset(
        new ScalarForcing(ctx, "surface.debm_simple.paleo", "obliquity", "radian", "degree", "obliquity of the earth"));

    m_perihelion_longitude.reset(
        new ScalarForcing(ctx, "surface.debm_simple.paleo", "long_peri", "radian", "degree", "longitude of the perihelion relative to the vernal equinox"));
  } else {
    m_use_paleo_file = false;
  }
}


/*! @brief The number of points for temperature and precipitation time-series.
 */
unsigned int DEBMSimplePointwise::timeseries_length(double dt) {
  double dt_years = dt / m_year_length;

  return std::max(1U, static_cast<unsigned int>(ceil(m_n_per_year * dt_years)));
}

/*!
 * The integrand in equation 6 of
 *
 * R. Calov and R. Greve, “A semi-analytical solution for the positive degree-day model
 * with stochastic temperature variations,” Journal of Glaciology, vol. 51, Art. no. 172,
 * 2005.
 *
 */
double DEBMSimplePointwise::CalovGreveIntegrand(double sigma, double TacC) {

  if (sigma == 0) {
    return std::max(TacC, 0.0);
  } else {
    const double Z = TacC / (sqrt(2.0) * sigma);
    return (sigma / sqrt(2.0 * M_PI)) * exp(-Z * Z) + (TacC / 2.0) * erfc(-Z);
  }
}

/*! Albedo parameterized as a function of the melt rate
 *
 * See equation 7 in Zeitz et al.
 *
 * @param[in] melt_rate melt amount (meters (ice equivalent) per second)
 * @param[in] cell_type cell type mask (used to exclude ice free areas)
 */
double DEBMSimplePointwise::albedo(double melt_rate, MaskValue cell_type) {
  // melt has a unit of meters ice equivalent
  //
  // dt has a unit of seconds
  if (cell_type == MASK_ICE_FREE_OCEAN) {
    return m_albedo_ocean;
  }

  if (cell_type == MASK_ICE_FREE_BEDROCK) {
    return m_albedo_land;
  }

  double result = m_albedo_snow + m_albedo_slope * melt_rate * m_ice_density ;
  return std::max(result, m_albedo_ice);
}


/*! Returns atmosphere transmissivity
 *
 * Note: it has no units and acts as a scaling factor.
 *
 * See appendix A2 in Zeitz et al 2021.
 *
 * @param[in] elevation elevation above the geoid (meters)
 */
double DEBMSimplePointwise::atmosphere_transmissivity(double elevation) {
  // transmissivity of the atmosphere (linear fit)
  return m_tau_a_intercept + m_tau_a_slope * elevation;
}


/*!
 * Returns the hour angle at which the sun reaches phi (for melting period during the day)
 *
 * Implements equation (2) in Zeitz et al 2021 (solved for h_{\Phi}).
 *
 * This equation goes back to equation 10 in Krebs-Kanzow 2018.
 *
 * @param[in] phi angle (radians)
 * @param[in] latitude latitude (radians)
 * @param[in] declination solar declination angle (radians)
 */
double DEBMSimplePointwise::get_h_phi(double phi, double latitude, double declination) {
  double cos_h_phi = (sin(phi) - sin(latitude) * sin(declination)) / (cos(latitude) * cos(declination));
  return acos(pism::clip(cos_h_phi, -1.0, 1.0));
}


/*!
 * Returns average top of atmosphere insolation during the daily melt period.
 *
 * Implements equation 5 in Zeitz et al (FIXME -- maybe???)
 *
 * See also 2.2.21 in Liou
 *
 * @param[in] distance2 FIXME
 * @param[in] h_phi FIXME
 * @param[in] lat latitude (radians)
 * @param[in] delta FIXME
 */
double DEBMSimplePointwise::get_q_insol(double distance2, double h_phi, double lat, double delta) {
  if (h_phi == 0) {
    return 0.;
  } else {
    double tmp = (h_phi * sin(lat) * sin(delta) + cos(lat) * cos(delta) * sin(h_phi));
    return m_solar_constant * distance2 * tmp / h_phi;
  }
}

//!
/* compute diurnal melt scheme  by equation (6) by Uta Krebs-Kanzow et al., The Cryosphere, 2018
 *
 * Implements equation (1) in Zeitz et al
 *
 * output in mm water equivalent (FIXME???)
 */
DEBMSimplePointwise::Melt DEBMSimplePointwise::calculate_melt(double time,
                                                              double dt,
                                                              double T_std_deviation,
                                                              double T,
                                                              double surface_elevation,
                                                              double latitude,
                                                              double albedo) {
  assert(dt > 0.0);

  double latitude_rad = (latitude / 180.0) * M_PI;

  double declination = 0.0;
  double distance2   = 0.0;
  if (m_paleo) {
    declination      = solar_declination_paleo(time);
    distance2        = distance_factor_paleo(time);
  } else {
    declination      = solar_declination(time);
    distance2        = distance_factor(time);
  }

  double tau_a   = atmosphere_transmissivity(surface_elevation);
  double h_phi   = get_h_phi(m_phi, latitude_rad, declination);
  double q_insol = get_q_insol(distance2, h_phi, latitude_rad, declination);

  double Teff = CalovGreveIntegrand(T_std_deviation, T - m_positive_threshold_temperature);
  if (Teff < 1.e-4) {
    Teff = 0;
  }

  double A = dt * (h_phi / M_PI / (m_water_density * m_L));

  Melt result;
  result.transmissivity = tau_a;
  result.q_insol        = q_insol;
  result.insolation_melt  = A * (tau_a * (1. - albedo) * q_insol);
  result.temperature_melt = A * m_c1 * Teff;
  result.background_melt  = A * m_c2;
  result.total_melt = result.insolation_melt + result.temperature_melt + result.background_melt;

  if (T < m_bm_temp) {
    result.total_melt = 0.0;
  }

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
 * @param[in] T air temperature
 * @param[in,out] P precipitation rate
 */
void DEBMSimplePointwise::get_snow_accumulation(const std::vector<double> &T, std::vector<double> &P) {

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


//! \brief Compute the surface mass balance at a location from the amount of
//! melted snow and the accumulation amount in a time interval.
/*!
 *
 * - a fraction of the melted snow and ice refreezes, conceptualized
 *   as superimposed ice
 *
 * - the excess of 'ITM_melt' is used to melt both the ice that came from refreeze and
 *   then any ice which is already present. (FIXME - I don't think this is true.)
 */
DEBMSimplePointwise::Changes DEBMSimplePointwise::step(double thickness, double input_melt,
                                                       double old_firn_depth, double old_snow_depth, double accumulation) {
  Changes result;

  double
    firn_depth      = old_firn_depth,
    snow_depth      = old_snow_depth,
    max_snow_melted = input_melt,
    firn_melted     = 0.0,
    snow_melted     = 0.0;

  assert(thickness >= 0);

  // snow depth cannot exceed total thickness
  snow_depth = std::min(snow_depth, thickness);

  assert(snow_depth >= 0);

  // firn depth cannot exceed thickness - snow_depth
  firn_depth = std::min(firn_depth, thickness - snow_depth);

  assert(firn_depth >= 0);

  snow_depth += accumulation;

  if (input_melt <= 0.0) { // The "no melt" case.
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

  double
    melt                    = input_melt,
    ice_created_by_refreeze = 0.0;

  if (m_refreeze_ice_melt) {
    // FIXME: incorrect in the case when all ice melted and there was excess energy that
    // could have melted more
    ice_created_by_refreeze = melt * m_refreeze_fraction;
  } else {
    // Should this only be snow melted?
    ice_created_by_refreeze = (firn_melted + snow_melted) * m_refreeze_fraction;
  }

  snow_depth = std::max(snow_depth - snow_melted, 0.0);
  firn_depth = std::max(firn_depth - firn_melted, 0.0);

  double runoff = melt - ice_created_by_refreeze;
  double smb    = accumulation - runoff;

  result.firn_depth = firn_depth - old_firn_depth;
  result.snow_depth = snow_depth - old_snow_depth;
  result.melt       = melt;
  result.runoff     = runoff;
  result.smb        = thickness + smb >= 0 ? smb : -thickness;

  assert(thickness + result.smb >= 0);

  return result;
}

/*!
 * Eccentricity of the Earth’s orbit (no units).
 */
double DEBMSimplePointwise::eccentricity(double time) {
  if (m_use_paleo_file) {
    return m_eccentricity->value(time);
  }
  return m_constant_eccentricity;
}

/*!
 * Returns the obliquity of the ecliptic in radians.
 */
double DEBMSimplePointwise::obliquity(double time) {
  if (m_use_paleo_file) {
    return m_obliquity->value(time);
  }
  return m_constant_obliquity;
}

/*!
 * Returns the longitude of the perihelion in radians.
 */
double DEBMSimplePointwise::perihelion_longitude(double time) {
  if (m_use_paleo_file) {
    return m_perihelion_longitude->value(time);
  }
  return m_constant_perihelion_longitude;
}


/*!
 * The factor scaling top of atmosphere insolation during the melt period according to the
 * Earth's distance from the Sun.
 *
 * The returned value is `(d_bar / d)^2`, where `d_bar` is the average distance from the
 * Earth to the Sun and `d` is the *current* distance at a given time.
 *
 * Implements equation 2.2.9 from Liou (2002).
 *
 * Liou states: "Note that the factor (a/r)^2 never departs from the unity by more than
 * 3.5%." (`a/r` in Liou is equivalent to `d_bar/d` here.)
 */
double DEBMSimplePointwise::distance_factor(double time) {
  // These coefficients come from Table 2.2 in Liou 2002
  double
    a0 = 1.000110,
    a1 = 0.034221,
    a2 = 0.000719,
    b0 = 0.,
    b1 = 0.001280,
    b2 = 0.000077;

  double t = 2. * M_PI * m_time->year_fraction(time);

  return (a0 + b0 +
          a1 * cos(t) + b1 * sin(t) +
          a2 * cos(2. * t) + b2 * sin(2. * t));
}

/*!
 * Earth declination
 *
 * Implements equation 2.2.10 from Liou (2002)
 */
double DEBMSimplePointwise::solar_declination(double time) {
  // These coefficients come from Table 2.2 in Liou 2002
   double
     a0 = 0.006918,
     a1 = -0.399912,
     a2 = -0.006758,
     a3 = -0.002697,
     b0 = 0.,
     b1 = 0.070257,
     b2 = 0.000907,
     b3 = 0.000148;

  double t = 2. * M_PI * m_time->year_fraction(time);

  return (a0 + b0 +
          a1 * cos(t) + b1 * sin(t) +
          a2 * cos(2. * t) + b2 * sin(2. * t) +
          a3 * cos(3. * t) + b3 * sin(3. * t));
}


/*!
 * Return factor
 *
 * Implements equation A1 in Zeitz et al.
 *
 * See also equation 2.2.5 from Liou (2002).
 */
double DEBMSimplePointwise::distance_factor_paleo(double time) {
  double E   = eccentricity(time);
  double L_p = perihelion_longitude(time);
  double year_fraction = m_time->year_fraction(time);
  double lambda = solar_longitude(year_fraction, E, L_p);

  if (E == 1.0) {
    // protect from division by zero
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid eccentricity value: 1.0");
  }
  return pow((1.0 + E * cos(lambda - L_p)) / (1.0 - E * E), 2);
}


/*!
 * Solar declination
 *
 * Implements equation in the text just above equation A1 in Zeitz et al.
 *
 * See also equation 2.2.4 of Liou (2002).
 */
double DEBMSimplePointwise::solar_declination_paleo(double time) {
  double epsilon = obliquity(time);
  double lambda  = solar_longitude(m_time->year_fraction(time),
                                   eccentricity(time),
                                   perihelion_longitude(time));

  return asin(sin(epsilon) * sin(lambda));
}


/*!
 * Estimates solar longitude at current time in the year.
 *
 * @param[in] year_fraction year fraction (between 0 and 1)
 * @param[in] eccentricity eccentricity of the earth’s orbit
 * @param[in] perihelion_longitude perihelion longitude (radians)
 *
 * Implements equation A2 in Zeitz et al.
 */
double DEBMSimplePointwise::solar_longitude(double year_fraction,
                                            double eccentricity,
                                            double perihelion_longitude) {

  // Shortcuts to make formulas below easier to read:
  double E   = eccentricity;
  double L_p = perihelion_longitude;

  // lambda = 0 at March equinox (80th day of the year)
  double delta_lambda  = 2. * M_PI * (year_fraction - 80. / 365.);
  double beta          = sqrt(1 - E * E);

  double lambda_m = (-2. * ((E / 2. + (pow(E, 3)) / 8.) * (1. + beta) * sin(-L_p) -
                           (pow(E, 2)) / 4. * (1. / 2. + beta) * sin(-2. * L_p) +
                           (pow(E, 3)) / 8. * (1. / 3. + beta) * sin(-3. * L_p)) +
                     delta_lambda);

  return (lambda_m + (2. * E - (pow(E, 3)) / 4.) * sin(lambda_m - L_p) +
          (5. / 4.) * (E * E) * sin(2. * (lambda_m - L_p)) +
          (13. / 12.) * (pow(E, 3)) * sin(3. * (lambda_m - L_p)));
}

} // end of namespace surface
} // end of namespace pism
