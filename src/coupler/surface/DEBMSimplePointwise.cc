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

/*!
 * This class implements dEBM-simple, the simple diurnal energy balance model described in
 *
 * M. Zeitz, R. Reese, J. Beckmann, U. Krebs-Kanzow, and R. Winkelmann, “Impact of the
 * melt–albedo feedback on the future evolution of the Greenland Ice Sheet with
 * PISM-dEBM-simple,” The Cryosphere, vol. 15, Art. no. 12, Dec. 2021.
 *
 * See also
 *
 * U. Krebs-Kanzow, P. Gierz, and G. Lohmann, “Brief communication: An ice surface melt
 * scheme including the diurnal cycle of solar radiation,” The Cryosphere, vol. 12, Art.
 * no. 12, Dec. 2018.
 *
 * and chapter 2 of
 *
 * K. N. Liou, Introduction to Atmospheric Radiation. Elsevier Science & Technology Books, 2002.
 *
 */
namespace pism {
namespace surface {

namespace details {
// Disable clang-tidy warnings about "magic numbers":
// NOLINTBEGIN(readability-magic-numbers)

/*!
 * The integrand in equation 6 of
 *
 * R. Calov and R. Greve, “A semi-analytical solution for the positive degree-day model
 * with stochastic temperature variations,” Journal of Glaciology, vol. 51, Art. no. 172,
 * 2005.
 *
 * @param[in] sigma standard deviation of daily variation of near-surface air temperature (Kelvin)
 * @param[in] temperature near-surface air temperature in "degrees Kelvin above the melting point"
 */
static double CalovGreveIntegrand(double sigma, double temperature) {

  if (sigma == 0) {
    return std::max(temperature, 0.0);
  }

  double Z = temperature / (sqrt(2.0) * sigma);
  return (sigma / sqrt(2.0 * M_PI)) * exp(-Z * Z) + (temperature / 2.0) * erfc(-Z);
}

/*!
 * The hour angle (radians) at which the sun reaches the solar angle `phi`
 *
 * Implements equation 11 in Krebs-Kanzow et al solved for h_phi.
 *
 * Equation 2 in Zeitz et al should be equivalent but misses "acos(...)".
 *
 * The return value is in the range [0, pi].
 *
 * @param[in] phi angle (radians)
 * @param[in] latitude latitude (radians)
 * @param[in] declination solar declination angle (radians)
 */
static double hour_angle(double phi, double latitude, double declination) {
  double cos_h_phi = ((sin(phi) - sin(latitude) * sin(declination)) /
                      (cos(latitude) * cos(declination)));
  return acos(pism::clip(cos_h_phi, -1.0, 1.0));
}

/*!
 * Solar longitude (radians) at current time in the year.
 *
 * @param[in] year_fraction year fraction (between 0 and 1)
 * @param[in] eccentricity eccentricity of the earth’s orbit (no units)
 * @param[in] perihelion_longitude perihelion longitude (radians)
 *
 * Implements equation A2 in Zeitz et al.
 */
static double solar_longitude(double year_fraction,
                              double eccentricity,
                              double perihelion_longitude) {

  // Shortcuts to make formulas below easier to read:
  double E   = eccentricity;
  double E2  = E * E;
  double E3  = E * E * E;
  double L_p = perihelion_longitude;

  // Note: lambda = 0 at March equinox (80th day of the year)
  const double equinox_day_number = 80.0;
  double delta_lambda  = 2.0 * M_PI * (year_fraction - equinox_day_number / 365.0);
  double beta          = sqrt(1.0 - E2);

  double lambda_m = (-2.0 * ((E / 2.0 + E3 / 8.0) * (1.0 + beta) * sin(-L_p) -
                             E2 / 4.0 * (1.0 / 2.0 + beta) * sin(-2.0 * L_p) +
                             E3 / 8.0 * (1.0 / 3.0 + beta) * sin(-3.0 * L_p)) +
                     delta_lambda);

  return (lambda_m +
          (2.0 * E - E3 / 4.0) * sin(lambda_m - L_p) +
          (5.0 / 4.0)   * E2 * sin(2.0 * (lambda_m - L_p)) +
          (13.0 / 12.0) * E3 * sin(3.0 * (lambda_m - L_p)));
}

/*!
 * The unit-less factor scaling top of atmosphere insolation according to the Earth's
 * distance from the Sun.
 *
 * The returned value is `(d_bar / d)^2`, where `d_bar` is the average distance from the
 * Earth to the Sun and `d` is the *current* distance at a given time.
 *
 * Implements equation 2.2.9 from Liou (2002).
 *
 * Liou states: "Note that the factor (a/r)^2 never departs from the unity by more than
 * 3.5%." (`a/r` in Liou is equivalent to `d_bar/d` here.)
 */
static double distance_factor_present_day(double year_fraction) {
  // These coefficients come from Table 2.2 in Liou 2002
  double
    a0 = 1.000110,
    a1 = 0.034221,
    a2 = 0.000719,
    b0 = 0.,
    b1 = 0.001280,
    b2 = 0.000077;

  double t = 2. * M_PI * year_fraction;

  return (a0 + b0 +
          a1 * cos(t) + b1 * sin(t) +
          a2 * cos(2. * t) + b2 * sin(2. * t));
}

/*!
 * The unit-less factor scaling top of atmosphere insolation according to the Earth's
 * distance from the Sun. This is the "paleo" version used when the trigonometric
 * expansion (equation 2.2.9 in Liou 2002) is not valid.
 *
 * Implements equation A1 in Zeitz et al.
 *
 * See also equation 2.2.5 from Liou (2002).
 */
static double distance_factor_paleo(double eccentricity,
                                    double perihelion_longitude,
                                    double solar_longitude) {
  double E = eccentricity;

  if (E == 1.0) {
    // protect from division by zero
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid eccentricity value: 1.0");
  }

  return pow((1.0 + E * cos(solar_longitude - perihelion_longitude)) / (1.0 - E * E), 2);
}

/*!
 * Solar declination (radian)
 *
 * Implements equation 2.2.10 from Liou (2002)
 */
static double solar_declination_present_day(double year_fraction) {
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

  double t = 2. * M_PI * year_fraction;

  return (a0 + b0 +
          a1 * cos(t) + b1 * sin(t) +
          a2 * cos(2. * t) + b2 * sin(2. * t) +
          a3 * cos(3. * t) + b3 * sin(3. * t));
}

/*!
 * Solar declination (radians). This is the "paleo" version used when
 * the trigonometric expansion (equation 2.2.10 in Liou 2002) is not valid.
 *
 * The return value is in the range [-pi/2, pi/2].
 *
 * Implements equation in the text just above equation A1 in Zeitz et al.
 *
 * See also equation 2.2.4 of Liou (2002).
 */
static double solar_declination_paleo(double obliquity,
                                      double solar_longitude) {
  return asin(sin(obliquity * sin(solar_longitude)));
}

/*!
 * Average top of atmosphere insolation (rate) during the daily melt period, in W/m^2.
 *
 * This should be equation 5 in Zeitz et al or equation 12 in Krebs-Kanzow et al, but both
 * of these miss a factor of Delta_t (day length in seconds) in the numerator.
 *
 * To confirm this, see the derivation of equation 2.2.21 in Liou and note that
 *
 * omega = 2 * pi (radian/day)
 *
 * or
 *
 * omega = (2 * pi / 86400) (radian/second).
 *
 * The correct equation should say
 *
 * S_Phi = A * B^2 * (h_phi * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(h_phi)),
 *
 * where
 *
 * A = (S0 * Delta_t) / (Delta_t_Phi * pi),
 * B = d_bar / d.
 *
 * Note that we do not know Delta_t_phi but we can use equation 2 in Zeitz et al (or
 * equation 11 in Krebs-Kanzow et al) to get
 *
 * Delta_t_phi = h_phi * Delta_t / pi.
 *
 * This gives
 *
 * S_Phi = C * B^2 * (h_phi * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(h_phi))
 *
 * with
 *
 * C = (S0 * Delta_t * pi) / (h_phi * Delta_t * pi)
 *
 * or
 *
 * C = S0 / h_phi.
 *
 * @param[in] solar constant solar constant, W/m^2
 * @param[in] distance_factor square of the ratio of the mean sun-earth distance to the current sun-earth distance (no units)
 * @param[in] hour_angle hour angle (radians) when the sun reaches the critical angle Phi
 * @param[in] latitude latitude (radians)
 * @param[in] declination declination (radians)
 *
 */
static double insolation(double solar_constant,
                         double distance_factor,
                         double hour_angle,
                         double latitude,
                         double declination) {
  if (hour_angle == 0) {
    return 0.0;
  }

  return ((solar_constant / hour_angle) * distance_factor *
          (hour_angle * sin(latitude) * sin(declination) +
           cos(latitude) * cos(declination) * sin(hour_angle)));
}

// NOLINTEND(readability-magic-numbers)
} // end of namespace details

DEBMSimplePointwise::Changes::Changes() {
  snow_depth = 0.0;
  melt       = 0.0;
  runoff     = 0.0;
  smb        = 0.0;
}

DEBMSimpleMelt::DEBMSimpleMelt() {
  temperature_melt = 0.0;
  insolation_melt  = 0.0;
  background_melt  = 0.0;
  total_melt       = 0.0;
  insolation       = 0.0;
}

DEBMSimplePointwise::DEBMSimplePointwise(const Context &ctx) {

  const Config &config = *ctx.config();

  m_time = ctx.time();

  m_L                              = config.get_number("constants.fresh_water.latent_heat_of_fusion");
  m_albedo_ice                     = config.get_number("surface.debm_simple.albedo_ice");
  m_albedo_ocean                   = config.get_number("surface.debm_simple.albedo_ocean");
  m_albedo_slope                   = config.get_number("surface.debm_simple.albedo_slope");
  m_albedo_snow                    = config.get_number("surface.debm_simple.albedo_snow");
  m_bm_temp                        = config.get_number("surface.debm_simple.background_melting_temp");
  m_c1                             = config.get_number("surface.debm_simple.c1");
  m_c2                             = config.get_number("surface.debm_simple.c2");
  m_constant_eccentricity          = config.get_number("surface.debm_simple.paleo.eccentricity");
  m_constant_obliquity             = config.get_number("surface.debm_simple.paleo.obliquity", "radian");
  m_constant_perihelion_longitude  = config.get_number("surface.debm_simple.paleo.perihelion_longitude", "radian");
  m_paleo                          = config.get_flag("surface.debm_simple.paleo.enabled");
  m_phi                            = config.get_number("surface.debm_simple.phi", "radian");
  m_positive_threshold_temperature = config.get_number("surface.debm_simple.positive_threshold_temp");
  m_refreeze_fraction              = config.get_number("surface.debm_simple.refreeze");
  m_refreeze_ice_melt              = config.get_flag("surface.debm_simple.refreeze_ice_melt");
  m_solar_constant                 = config.get_number("surface.debm_simple.solar_constant");
  m_transmissivity_intercept       = config.get_number("surface.debm_simple.tau_a_intercept");
  m_transmissivity_slope           = config.get_number("surface.debm_simple.tau_a_slope");

  m_ice_density   = config.get_number("constants.ice.density");
  m_water_density = config.get_number("constants.fresh_water.density");

  assert(m_albedo_slope < 0.0);
  assert(m_ice_density > 0.0);

  std::string paleo_file = config.get_string("surface.debm_simple.paleo.file");

  if (not paleo_file.empty()) {
    m_use_paleo_file = true;

    m_eccentricity.reset(
        new ScalarForcing(ctx, "surface.debm_simple.paleo", "eccentricity", "", "", "eccentricity of the earth"));

    m_obliquity.reset(
        new ScalarForcing(ctx, "surface.debm_simple.paleo", "obliquity", "radian", "degree", "obliquity of the earth"));

    m_perihelion_longitude.reset(
        new ScalarForcing(ctx, "surface.debm_simple.paleo", "perihelion_longitude", "radian", "degree", "longitude of the perihelion relative to the vernal equinox"));
  } else {
    m_use_paleo_file = false;
  }
}

/*! Albedo parameterized as a function of the melt rate
 *
 * See equation 7 in Zeitz et al.
 *
 * @param[in] melt_rate melt rate (meters (liquid water equivalent) per second)
 * @param[in] cell_type cell type mask (used to exclude ice free areas)
 */
double DEBMSimplePointwise::albedo(double melt_rate, MaskValue cell_type) const {
  if (cell_type == MASK_ICE_FREE_OCEAN) {
    return m_albedo_ocean;
  }

  assert(melt_rate >= 0.0);

  double result = m_albedo_snow + m_albedo_slope * melt_rate * m_ice_density ;
  return std::max(result, m_albedo_ice);
}

/*! Atmosphere transmissivity (no units; acts as a scaling factor)
 *
 * See appendix A2 in Zeitz et al 2021.
 *
 * @param[in] elevation elevation above the geoid (meters)
 */
double DEBMSimplePointwise::atmosphere_transmissivity(double elevation) const {
  return m_transmissivity_intercept + m_transmissivity_slope * elevation;
}

/* Melt amount (in m water equivalent) and its components over the time step `dt`
 *
 * Implements equation (1) in Zeitz et al.
 *
 * See also the equation (6) in Krebs-Kanzow et al.
 *
 * @param[in] time current time (seconds)
 * @param[in] dt time step length (seconds)
 * @param[in] T_std_deviation standard deviation of the near-surface air temperature (Kelvin)
 * @param[in] T near-surface air temperature (Kelvin)
 * @param[in] surface_elevation surface elevation (meters)
 * @param[in] latitude latitude (degrees north)
 * @param[in] albedo current albedo (fraction)
 */
DEBMSimpleMelt DEBMSimplePointwise::melt(double time,
                                         double dt,
                                         double T_std_deviation,
                                         double T,
                                         double surface_elevation,
                                         double latitude,
                                         double albedo) const {
  assert(dt > 0.0);

  const double degrees_to_radians = M_PI / 180.0;
  double latitude_rad = latitude * degrees_to_radians;

  double declination     = 0.0;
  double distance_factor = 0.0;
  double year_fraction = m_time->year_fraction(time);
  if (m_paleo) {
    double eccentricity         = this->eccentricity(time);
    double perihelion_longitude = this->perihelion_longitude(time);

    double solar_longitude = details::solar_longitude(year_fraction,
                                                      eccentricity,
                                                      perihelion_longitude);
    {
      double obliquity = this->obliquity(time);
      declination = details::solar_declination_paleo(obliquity, solar_longitude);
    }

    distance_factor = details::distance_factor_paleo(eccentricity,
                                                     perihelion_longitude,
                                                     solar_longitude);
  } else {
    declination          = details::solar_declination_present_day(year_fraction);
    distance_factor      = details::distance_factor_present_day(year_fraction);
  }

  double transmissivity = atmosphere_transmissivity(surface_elevation);
  double h_phi          = details::hour_angle(m_phi, latitude_rad, declination);
  double insolation     = details::insolation(m_solar_constant,
                                              distance_factor,
                                              h_phi,
                                              latitude_rad,
                                              declination);

  double Teff = details::CalovGreveIntegrand(T_std_deviation,
                                             T - m_positive_threshold_temperature);
  const double eps = 1.0e-4;
  if (Teff < eps) {
    Teff = 0;
  }

  // Note that in the line below we replace "Delta_t_Phi / Delta_t" with "h_Phi / pi". See
  // equations 1 and 2 in Zeitz et al.
  double A = dt * (h_phi / M_PI / (m_water_density * m_L));

  DEBMSimpleMelt result;

  result.insolation       = insolation;
  result.insolation_melt  = A * (transmissivity * (1.0 - albedo) * insolation);
  result.temperature_melt = A * m_c1 * Teff;
  result.background_melt  = A * m_c2;

  double total_melt = (result.insolation_melt + result.temperature_melt +
                       result.background_melt);
  // this model should not produce negative melt rates
  result.total_melt = std::max(total_melt, 0.0);

  if (T < m_bm_temp) {
    result.total_melt = 0.0;
  }

  return result;
}

/*! @brief Compute the surface mass balance at a location from the amount of melted snow
 *  and the solid accumulation amount in a time interval.
 *
 * - a fraction of the melted snow and ice refreezes, conceptualized
 *   as superimposed ice
 */
DEBMSimplePointwise::Changes DEBMSimplePointwise::step(double ice_thickness,
                                                       double max_melt,
                                                       double old_snow_depth,
                                                       double accumulation) const {
  Changes result;

  double
    snow_depth      = old_snow_depth,
    snow_melted     = 0.0,
    ice_melted      = 0.0;

  assert(ice_thickness >= 0);

  // snow depth cannot exceed total ice_thickness
  snow_depth = std::min(snow_depth, ice_thickness);

  assert(snow_depth >= 0);

  snow_depth += accumulation;

  if (max_melt <= 0.0) { // The "no melt" case.
    snow_melted = 0.0;
    ice_melted  = 0.0;
  } else if (max_melt <= snow_depth) {
    // Some of the snow melted and some is left; in any case, all of the energy available
    // for melt was used up in melting snow.
    snow_melted = max_melt;
    ice_melted  = 0.0;
  } else {
    // All (snow_depth meters) of snow melted. Excess melt is available to melt ice.
    snow_melted = snow_depth;
    ice_melted  = std::min(max_melt - snow_melted, ice_thickness);
  }

  double ice_created_by_refreeze = m_refreeze_fraction * snow_melted;
  if (m_refreeze_ice_melt) {
    ice_created_by_refreeze += m_refreeze_fraction * ice_melted;
  }

  snow_depth = std::max(snow_depth - snow_melted, 0.0);

  double total_melt = (snow_melted + ice_melted);
  double runoff     = total_melt - ice_created_by_refreeze;
  double smb        = accumulation - runoff;

  result.snow_depth = snow_depth - old_snow_depth;
  result.melt       = total_melt;
  result.runoff     = runoff;
  result.smb        = ice_thickness + smb >= 0 ? smb : -ice_thickness;

  assert(ice_thickness + result.smb >= 0);

  return result;
}

/*!
 * Eccentricity of the Earth’s orbit (no units).
 */
double DEBMSimplePointwise::eccentricity(double time) const {
  if (m_use_paleo_file) {
    return m_eccentricity->value(time);
  }
  return m_constant_eccentricity;
}

/*!
 * Returns the obliquity of the ecliptic in radians.
 */
double DEBMSimplePointwise::obliquity(double time) const {
  if (m_use_paleo_file) {
    return m_obliquity->value(time);
  }
  return m_constant_obliquity;
}

/*!
 * Returns the longitude of the perihelion in radians.
 */
double DEBMSimplePointwise::perihelion_longitude(double time) const {
  if (m_use_paleo_file) {
    return m_perihelion_longitude->value(time);
  }
  return m_constant_perihelion_longitude;
}

} // end of namespace surface
} // end of namespace pism
