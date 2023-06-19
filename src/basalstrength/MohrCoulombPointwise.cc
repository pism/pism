/* Copyright (C) 2019, 2023 PISM Authors
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
#include <cmath>                // pow, tan, atan

#include "pism/basalstrength/MohrCoulombPointwise.hh"

namespace pism {

MohrCoulombPointwise::MohrCoulombPointwise(Config::ConstPtr config) {
  m_W_till_max                   = config->get_number("hydrology.tillwat_max");
  m_till_cohesion                = config->get_number("basal_yield_stress.mohr_coulomb.till_cohesion");
  m_reference_effective_pressure = config->get_number("basal_yield_stress.mohr_coulomb.till_reference_effective_pressure");
  m_reference_void_ratio         = config->get_number("basal_yield_stress.mohr_coulomb.till_reference_void_ratio");
  m_compressibility_coefficient  = config->get_number("basal_yield_stress.mohr_coulomb.till_compressibility_coefficient");
}

double MohrCoulombPointwise::effective_pressure(double delta,
                                                double P_overburden,
                                                double water_thickness) const {

  double
    s      = water_thickness / m_W_till_max,
    N0     = m_reference_effective_pressure,
    N_till = (N0 * pow(delta * P_overburden / N0, s) *
              pow(10.0, (m_reference_void_ratio / m_compressibility_coefficient) * (1.0 - s)));

  return std::min(P_overburden, N_till);
}

double MohrCoulombPointwise::yield_stress(double delta,
                                          double P_overburden,
                                          double water_thickness,
                                          double phi) const {

  double N_till = effective_pressure(delta, P_overburden, water_thickness);

  return m_till_cohesion + N_till * tan((M_PI / 180.0) * phi);
}

double MohrCoulombPointwise::till_friction_angle(double delta,
                                                 double P_overburden,
                                                 double water_thickness,
                                                 double yield_stress) const {

  double N_till = effective_pressure(delta, P_overburden, water_thickness);

  return 180.0 / M_PI * atan((yield_stress - m_till_cohesion) / N_till);

}

} // end of namespace pism
