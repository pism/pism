/* Copyright (C) 2018 PISM Authors
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

namespace pism {

class Config;

namespace ocean {

class PicopPhysics {
public:
  PicopPhysics(const Config &config);

  double characteristic_freezing_point(const double s_a, const double z_gl) const;

  double effective_heat_exchange_coefficient(const double t_a, const double t_f_gl, const double alpha) const;

  double length_scaling(const double t_a, const double t_f_gl, const double GammaTS, const double alpha) const;
  
  double geometric_scaling(const double GammaTS, const double alpha) const;

  double dimensionless_coordinate(const double z_b, const double z_gl, const double l) const;

  double dimensionless_melt_curve(const double X_hat) const;

  double melt_function(const double t_a, const double t_f_gl, const double g_alpha) const;
  
  double melt_rate(const double M, const double X_hat) const;
  
private:

  // drag coefficient
  double m_Cd, m_Cd12;
  double m_GammaT, m_Cd12GammaT;
  // heat exchange parameters
  double m_Cd12GammaTS0, m_gamma1, m_gamma2;
  // freezing point coefficients
  double m_lambda1, m_lambda2, m_lambda3;
  // entrainment coefficient
  double m_E0, m_M0;
  // scaling coefficients
  double  m_x0;

};

} // end of namespace ocean
} // end of namespace pism
