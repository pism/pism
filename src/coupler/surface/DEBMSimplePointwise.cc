// Copyright (C) 2009--2026 PISM Authors
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
#include <cmath>

#include "pism/coupler/surface/DEBMSimplePointwise.hh"
#include "pism/util/Config.hh"
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

// Disable clang-tidy warnings about "magic numbers":
// NOLINTBEGIN(readability-magic-numbers)

/*!
 * The integrand in equation 6 of
 *
 * R. Calov and R. Greve, “A semi-analytical solution for the positive degree-day model
 * with stochastic temperature variations,” Journal of Glaciology, vol. 51, Art. no. 172,
 * 2005.
 *
 * @param[in] sigma standard deviation of daily variation of near-surface air temperature (kelvin)
 * @param[in] temperature near-surface air temperature in "kelvin above the melting point"
 */
double DEBMSimplePointwise::CalovGreveIntegrand(double sigma, double temperature) {

  if (sigma == 0) {
    return std::max(temperature, 0.0);
  }

  double Z = temperature / (sqrt(2.0) * sigma);
  return (sigma / sqrt(2.0 * M_PI)) * exp(-Z * Z) + (temperature / 2.0) * erfc(-Z);
}

/*!
 * The hour angle (radians) at which the sun reaches the solar altitude angle `phi`
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
double DEBMSimplePointwise::hour_angle(double phi, double latitude, double declination) {
  double cos_h_phi = ((sin(phi) - sin(latitude) * sin(declination)) /
                      (cos(latitude) * cos(declination)));
  return acos(pism::clip(cos_h_phi, -1.0, 1.0));
}


/*!
 * Average top of atmosphere insolation (rate) during the daily melt period, in W/m^2.
 *
 * This should be equation 5 in Zeitz et al or equation 12 in Krebs-Kanzow et al, but both
 * of these miss a factor of Delta_t (day length in seconds) in the numerator. (Note that
 * in equation 5 of Zeitz et al the units of the LHS and the RHS should match, but [LHS] =
 * W/m^2 and [RHS] = W/(m^2 s).)
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
 * @param[in] solar_constant solar constant, W/m^2
 * @param[in] distance_factor square of the ratio of the mean sun-earth distance to the current sun-earth distance (no units)
 * @param[in] hour_angle hour angle (radians) when the sun reaches the critical angle Phi
 * @param[in] latitude latitude (radians)
 * @param[in] declination declination (radians)
 *
 */
double DEBMSimplePointwise::insolation_rate(double solar_constant, double distance_factor,
                                            double hour_angle, double latitude,
                                            double declination) {
  if (hour_angle == 0) {
    return 0.0;
  }

  return ((solar_constant / hour_angle) * distance_factor *
          (hour_angle * sin(latitude) * sin(declination) +
           cos(latitude) * cos(declination) * sin(hour_angle)));
}

// NOLINTEND(readability-magic-numbers)

DEBMSimpleChanges::DEBMSimpleChanges() {
  snow_depth = 0.0;
  melt       = 0.0;
  runoff     = 0.0;
  smb        = 0.0;
}

DEBMSimpleMelt::DEBMSimpleMelt() {
  temperature_melt = 0.0;
  insolation_melt  = 0.0;
  offset_melt  = 0.0;
  total_melt       = 0.0;
}

DEBMSimplePointwise::DEBMSimplePointwise(const Config &config) : m_transmissivity(config) {

  m_L                              = config.get_number("constants.fresh_water.latent_heat_of_fusion");
  m_albedo_min                     = config.get_number("surface.debm_simple.albedo_min");
  m_albedo_slope                   = config.get_number("surface.debm_simple.albedo_slope");
  m_albedo_max                     = config.get_number("surface.debm_simple.albedo_max");
  m_melt_threshold_temp            = config.get_number("surface.debm_simple.melting_threshold_temp");
  m_melt_c1                        = config.get_number("surface.debm_simple.c1");
  m_melt_c2                        = config.get_number("surface.debm_simple.c2");
  m_phi                            = config.get_number("surface.debm_simple.phi", "radian");
  m_positive_threshold_temperature = config.get_number("surface.debm_simple.positive_threshold_temp");
  m_refreeze_fraction              = config.get_number("surface.debm_simple.refreeze");
  m_refreeze_ice_melt              = config.get_flag("surface.debm_simple.refreeze_ice_melt");
  m_solar_constant                 = config.get_number("surface.debm_simple.solar_constant");

  m_ice_density   = config.get_number("constants.ice.density");
  m_water_density = config.get_number("constants.fresh_water.density");

  assert(m_albedo_slope < 0.0);
  assert(m_ice_density > 0.0);

}

/*! Albedo parameterized as a function of the melt rate
 *
 * See equation 7 in Zeitz et al.
 *
 * @param[in] melt_rate melt rate (meters (liquid water equivalent) per second)
 */
double DEBMSimplePointwise::albedo(double melt_rate) const {
  assert(melt_rate >= 0.0);

  return std::max(m_albedo_max + m_albedo_slope * melt_rate * m_ice_density, //
                  m_albedo_min);
}

/*! Atmosphere transmissivity (no units; acts as a scaling factor)
 *
 * See appendix A2 in Zeitz et al 2021.
 *
 * @param[in] elevation elevation above the geoid (meters)
 */
double DEBMSimplePointwise::atmosphere_transmissivity(double elevation) const {
  return m_transmissivity(elevation);
}


/*!
 * Compute top of atmosphere insolation to report as a diagnostic quantity.
 *
 * Do not use this in the model itself: doing so will make it slower because that way we'd
 * end up computing hour_angle more than once.
 */
double DEBMSimplePointwise::insolation_diagnostic(double declination, double distance_factor,
                                                  double latitude_degrees) const {
  const double degrees_to_radians = M_PI / 180.0;
  double latitude_rad = latitude_degrees * degrees_to_radians;

  double h_phi = hour_angle(m_phi, latitude_rad, declination);

  return insolation_rate(m_solar_constant, distance_factor, h_phi, latitude_rad, declination);
}

/* Melt amount (in m water equivalent) and its components over the time step `dt`
 *
 * Implements equation (1) in Zeitz et al.
 *
 * See also the equation (6) in Krebs-Kanzow et al.
 *
 * @param[in] time current time (seconds)
 * @param[in] dt time step length (seconds)
 * @param[in] T_std_deviation standard deviation of the near-surface air temperature (kelvin)
 * @param[in] T near-surface air temperature (kelvin)
 * @param[in] surface_elevation surface elevation (meters)
 * @param[in] latitude latitude (degrees north)
 * @param[in] albedo current albedo (fraction)
 */
DEBMSimpleMelt DEBMSimplePointwise::melt(double declination,
                                         double distance_factor,
                                         double dt,
                                         double T_std_deviation,
                                         double T,
                                         double surface_elevation,
                                         double latitude,
                                         double albedo) const {
  // Compute the top-of-atmosphere insolation energy reaching the surface over `dt`.
  double E_toa = insolation_energy(declination, distance_factor, latitude, dt);

  double E_surface = atmosphere_transmissivity(surface_elevation) * E_toa;

  return melt_from_insolation(declination, latitude, E_surface, dt, T_std_deviation, T,
                              surface_elevation, albedo);
}

/*!
 * Analytic top-of-atmosphere insolation *energy* (J/m^2) reaching the surface over the
 * time step `dt`.
 *
 * This is the average insolation rate S_Phi (W/m^2) during the daily melt period times the
 * length of that melt period `dt`, i.e. `S_Phi * dt * (h_phi / pi)`.
 *
 * @param[in] declination solar declination (radians)
 * @param[in] distance_factor square of the ratio of the mean sun-earth distance to the current sun-earth distance (no units)
 * @param[in] latitude latitude (degrees)
 * @param[in] dt time-step length (seconds)
 */
double DEBMSimplePointwise::insolation_energy(double declination, double distance_factor,
                                              double latitude, double dt) const {
  const double degrees_to_radians = M_PI / 180.0;
  double latitude_rad = latitude * degrees_to_radians;

  double h_phi  = hour_angle(m_phi, latitude_rad, declination);
  double S_phi  = insolation_rate(m_solar_constant, distance_factor, h_phi, latitude_rad, declination);

  return S_phi * dt * (h_phi / M_PI);
}

/* Melt amount (in m water equivalent) and its components over the time step `dt`, given the
 * insolation *energy* (J/m^2) reaching the surface.
 *
 * Implements equation (1) in Zeitz et al; see also equation (6) in Krebs-Kanzow et al.
 *
 * @param[in] declination solar declination (radians), used for the daily melt-period length
 * @param[in] latitude latitude (degrees north)
 * @param[in] insolation_energy insolation energy reaching the surface over `dt` (J/m^2)
 * @param[in] dt time step length (seconds)
 * @param[in] T_std_deviation standard deviation of the near-surface air temperature (kelvin)
 * @param[in] T near-surface air temperature (kelvin)
 * @param[in] surface_elevation surface elevation (meters)
 * @param[in] albedo current albedo (fraction)
 */
DEBMSimpleMelt DEBMSimplePointwise::melt_from_insolation(double declination,
                                                         double latitude,
                                                         double insolation_energy,
                                                         double dt,
                                                         double T_std_deviation,
                                                         double T,
                                                         double surface_elevation,
                                                         double albedo) const {
  assert(dt > 0.0);

  const double degrees_to_radians = M_PI / 180.0;
  double latitude_rad = latitude * degrees_to_radians;

  double h_phi          = hour_angle(m_phi, latitude_rad, declination);

  double Teff = CalovGreveIntegrand(T_std_deviation,
                                    T - m_positive_threshold_temperature);
  const double eps = 1.0e-4;
  if (Teff < eps) {
    Teff = 0;
  }

  // Note that in the line below we replace "Delta_t_Phi / Delta_t" with "h_Phi / pi". See
  // equations 1 and 2 in Zeitz et al. This melt-period weight applies to the temperature-
  // and offset-driven melt terms.
  double A = dt * (h_phi / M_PI / (m_water_density * m_L));

  DEBMSimpleMelt result;

  // insolation_energy is the energy reaching the surface over dt; albedo converts it to
  // absorbed energy, and dividing by rho_w * L gives melt in meters of water equivalent.
  result.insolation_melt  = ((1.0 - albedo) * insolation_energy) /
                            (m_water_density * m_L);
  result.temperature_melt = A * m_melt_c1 * Teff;
  result.offset_melt      = A * m_melt_c2;

  double total_melt = (result.insolation_melt + result.temperature_melt +
                       result.offset_melt);

  // this model should not produce negative melt rates
  result.total_melt = std::max(total_melt, 0.0);

  if (T < m_melt_threshold_temp) {
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
DEBMSimpleChanges DEBMSimplePointwise::step(double ice_thickness, double max_melt,
                                            double old_snow_depth, double accumulation) const {
  DEBMSimpleChanges result;

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

DEBMSimpleAtmosphereTransmissivity::DEBMSimpleAtmosphereTransmissivity(const Config &config) {
  m_slope = config.get_number("surface.debm_simple.tau_a_slope");
  m_intercept = config.get_number("surface.debm_simple.tau_a_intercept");
}

/*! Atmosphere transmissivity (no units; acts as a scaling factor)
 *
 * See appendix A2 in Zeitz et al 2021.
 *
 * @param[in] surface_elevation elevation above the geoid (meters)
 */
double DEBMSimpleAtmosphereTransmissivity::operator()(double surface_elevation) const {
  return m_intercept + m_slope * surface_elevation;

}

} // end of namespace surface
} // end of namespace pism
