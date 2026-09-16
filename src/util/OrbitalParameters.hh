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

#ifndef PISM_ORBITALPARAMETERS_H
#define PISM_ORBITALPARAMETERS_H

#include <memory>

namespace pism {

class Context;
class Time;
class ScalarForcing;

struct OrbitalParameters {
  // Solar declination, radians
  double solar_declination;
  // Square of the ratio of the mean sun-earth distance to the current sun-earth distance
  // (d_bar / d)^2
  double distance_factor;
};

class OrbitalParameterCalculator {
public:
  OrbitalParameterCalculator(const Context &ctx);
  ~OrbitalParameterCalculator();

  OrbitalParameters compute(double time) const;

  static double solar_longitude(double year_fraction, double eccentricity,
                                double perihelion_longitude);
  static double distance_factor_present_day(double year_fraction);
  static double distance_factor_paleo(double eccentricity, double true_anomaly);
  static double solar_declination_present_day(double year_fraction);
  static double solar_declination_paleo(double obliquity,
                                        double solar_longitude);

  double eccentricity(double time) const;
  double obliquity(double time) const;
  double perihelion_longitude(double time) const;

private:
  std::unique_ptr<ScalarForcing> m_eccentricity;
  std::unique_ptr<ScalarForcing> m_obliquity;
  std::unique_ptr<ScalarForcing> m_perihelion_longitude;

  std::shared_ptr<const Time> m_time;

  bool m_paleo;

  double m_constant_eccentricity;
  double m_constant_perihelion_longitude;
  double m_constant_obliquity;
};

} // end of namespace pism

#endif /* PISM_ORBITALPARAMETERS_H */
