/* Copyright (C) 2026 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <cmath>

#include "pism/util/OrbitalParameters.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/ScalarForcing.hh"
#include "pism/util/Context.hh"
#include "pism/util/Time.hh"

namespace pism {
/*!
 * Solar longitude (radians) at current time in the year.
 *
 * @param[in] year_fraction year fraction (between 0 and 1)
 * @param[in] eccentricity eccentricity of the earth’s orbit (no units)
 * @param[in] perihelion_longitude perihelion longitude (radians) in the geocentric
 *            ecliptic coordinate system
 *
 * Implements equation A2 in Zeitz et al. or equivalent equations in Berger1978 (section
 * 3).
 */
double OrbitalParameterCalculator::solar_longitude(double year_fraction, double eccentricity,
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
 * The unit-less factor scaling top of atmosphere insolation according to the earth's
 * distance from the sun.
 *
 * The returned value is `(d_bar / d)^2`, where `d_bar` is the average distance from the
 * earth to the sun and `d` is the *current* distance at a given time.
 *
 * Implements equation 2.2.9 from Liou (2002).
 *
 * Liou states: "Note that the factor (a/r)^2 never departs from the unity by more than
 * 3.5%." (`a/r` in Liou is equivalent to `d_bar/d` here.)
 *
 * This quantity is equal to `1 / R^2`, where R is the sun-earth distance in units of the
 * semi-major axis of the earth's orbit.
 */
double OrbitalParameterCalculator::distance_factor_present_day(double year_fraction) {
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
 * The unit-less factor scaling top of atmosphere insolation according to the earth's
 * distance from the sun. This is the "paleo" version used when the trigonometric
 * expansion (equation 2.2.9 in Liou 2002) is not valid.
 *
 * Implements equation A1 in Zeitz et al.
 *
 * See also equation 2.2.5 from Liou (2002).
 *
 * This quantity is equal to `1 / R^2`, where R is the sun-earth distance in units of the
 * semi-major axis of the earth's orbit.
 *
 * @param[in] eccentricity eccentricity of the earth's orbit
 * @param[in] true_anomaly true anomaly of the earth in the heliocentric ecliptic coordinate system
 */
double OrbitalParameterCalculator::distance_factor_paleo(double eccentricity, double true_anomaly) {
  double E = eccentricity;

  if (E == 1.0) {
    // protect from division by zero
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid eccentricity value: 1.0");
  }

  return pow((1.0 + E * cos(true_anomaly)) / (1.0 - E * E), 2);
}

/*!
 * Solar declination (radian)
 *
 * Implements equation 2.2.10 from Liou (2002)
 */
double OrbitalParameterCalculator::solar_declination_present_day(double year_fraction) {
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
double OrbitalParameterCalculator::solar_declination_paleo(double obliquity,
                                                  double solar_longitude) {
  return asin(sin(obliquity) * sin(solar_longitude));
}

OrbitalParameterCalculator::OrbitalParameterCalculator(const Context &ctx) {
  const auto &config = *ctx.config();

  m_time = ctx.time();

  m_constant_eccentricity          = config.get_number("surface.debm_simple.paleo.eccentricity");
  m_constant_obliquity             = config.get_number("surface.debm_simple.paleo.obliquity", "radian");
  m_constant_perihelion_longitude  = config.get_number("surface.debm_simple.paleo.perihelion_longitude", "radian");
  m_paleo                          = config.get_flag("surface.debm_simple.paleo.enabled");

  std::string paleo_file = config.get_string("surface.debm_simple.paleo.file");

  if (not paleo_file.empty()) {
    m_eccentricity.reset(new ScalarForcing(ctx, "surface.debm_simple.paleo", "eccentricity", "1",
                                           "1", "eccentricity of the earth"));

    m_obliquity.reset(new ScalarForcing(ctx, "surface.debm_simple.paleo", "obliquity", "radian",
                                        "degree", "obliquity of the earth"));

    m_perihelion_longitude.reset(
        new ScalarForcing(ctx, "surface.debm_simple.paleo", "perihelion_longitude", "radian",
                          "degree", "longitude of the perihelion relative to the vernal equinox, "
                          "in the geocentric ecliptic coordinate system"));
  }
}

OrbitalParameterCalculator::~OrbitalParameterCalculator() = default;

OrbitalParameters OrbitalParameterCalculator::compute(double time) const {
  double solar_declination = 0.0;
  double distance_factor   = 0.0;

  double year_fraction = m_time->year_fraction(time);
  if (m_paleo) {
    double eccentricity                    = this->eccentricity(time);
    double geocentric_perihelion_longitude = this->perihelion_longitude(time);

    double sun_true_longitude =
        this->solar_longitude(year_fraction, eccentricity, geocentric_perihelion_longitude);

    solar_declination = solar_declination_paleo(obliquity(time), sun_true_longitude);

    // Note: here sun_true_longitude is the true longitude of the sun in the geocentric
    // ecliptic coordinate system. The perihelion longitude is in the same coordinate
    // system.
    //
    // In the heliocentric ecliptic system
    //
    // earth_true_longitude              = sun_true_longitude - 180 degrees
    // heliocentric_perihelion_longitude = geocentric_perihelion_longitude - 180 degrees
    //
    // So earth_true_anomaly = earth_true_longitude - heliocentric_perihelion_longitude
    //                       = sun_true_longitude - geocentric_perihelion_longitude
    //
    double earth_true_anomaly = sun_true_longitude - geocentric_perihelion_longitude;

    distance_factor = distance_factor_paleo(eccentricity, earth_true_anomaly);
  } else {
    solar_declination = solar_declination_present_day(year_fraction);
    distance_factor   = distance_factor_present_day(year_fraction);
  }

  return { solar_declination, distance_factor };
}

/*!
 * Eccentricity of the earth’s orbit (no units).
 */
double OrbitalParameterCalculator::eccentricity(double time) const {
  if (m_eccentricity != nullptr) {
    return m_eccentricity->value(time);
  }
  return m_constant_eccentricity;
}

/*!
 * Returns the obliquity of the ecliptic in radians.
 */
double OrbitalParameterCalculator::obliquity(double time) const {
  if (m_obliquity != nullptr) {
    return m_obliquity->value(time);
  }
  return m_constant_obliquity;
}

/*!
 * Returns the longitude of the perihelion (radians) in the geocentric ecliptic coordinate
 * system.
 */
double OrbitalParameterCalculator::perihelion_longitude(double time) const {
  if (m_perihelion_longitude != nullptr) {
    double L_p = remainder(m_perihelion_longitude->value(time), 2.0 * M_PI);
    if (L_p < 0.0) {
      L_p = L_p + 2 * M_PI;
    }
    return L_p;
  }
  return m_constant_perihelion_longitude;
}

} // end of namespace pism
