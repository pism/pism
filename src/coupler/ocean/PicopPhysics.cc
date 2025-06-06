/* Copyright (C) 2018, 2019, 2023, 2025 PISM Authors
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
#include <cmath> // sqrt
#include <cassert>              // assert
#include <vector>
#include <algorithm>            // std::min
#include <iostream>

#include "pism/util/Config.hh"
#include "pism/coupler/ocean/PicopPhysics.hh"

namespace pism {
namespace ocean {

PicopPhysics::PicopPhysics(const Config &config) {
  
    m_E0             = config.get_number("ocean.picop.entrainment_coefficient");
    m_Cd             = config.get_number("ocean.picop.drag_coefficient");
    m_Cd12GammaT     = config.get_number("ocean.picop.turbulent_heat_exchange_coefficient");
    m_Cd12GammaTS0   = config.get_number("ocean.picop.heat_exchange_parameter");
    m_gamma1         = config.get_number("ocean.picop.heat_exchange_parameter_1");
    m_gamma2         = config.get_number("ocean.picop.heat_exchange_parameter_2");
    m_lambda1        = config.get_number("ocean.picop.freezing_point_salinity_coefficient");
    m_lambda2        = config.get_number("ocean.picop.freezing_point_offset", "kelvin");
    m_lambda3        = config.get_number("ocean.picop.freezing_point_depth_coefficient", "kelvin m^-1");
    m_M0             = config.get_number("ocean.picop.melt_rate_parameter", "m s^-1 degree_Celsius^-2");
    m_x0             = config.get_number("ocean.picop.dimensionless_scaling_factor");

    m_Cd12 = sqrt(m_Cd);
    m_GammaT = m_Cd12GammaT / m_Cd12;  // (1)
}

//! equation 4 in the PICOP paper.
double PicopPhysics::characteristic_freezing_point(const double s_a, const double z) const {
  // in K * g /kg + K  + K / m
  return m_lambda1 * s_a + m_lambda2 + m_lambda3 * z;
}

//! equation 5 in the PICOP paper.
double PicopPhysics::effective_heat_exchange_coefficient(const double t_a, const double t_f_gl, const double alpha) const {
  const double t2 = (t_a - t_f_gl) / m_lambda3;
  const double num = m_E0 * sin(alpha);
  const double denom = (m_Cd12GammaTS0 + m_E0 * sin(alpha));
  return m_GammaT * (m_gamma1 + m_gamma2 * t2 * num / denom);
}

//! equation 6 in the PICOP paper.
double PicopPhysics::geometric_scaling(const double GammaTS, const double alpha) const {
  // 1 * 1 * 1
  const double g1 = sqrt((sin(alpha) / (m_Cd + m_E0 * sin(alpha))));
  const double g2 = sqrt((m_E0 * sin(alpha)) / (m_Cd12 * GammaTS + m_E0 * sin(alpha)));
  const double g3 = (m_E0 * sin(alpha)) / (m_Cd12 * GammaTS + m_E0 * sin(alpha));
  return g1 * g2 * g3;
}

//! equation 7 in the PICOP paper.
double PicopPhysics::length_scaling(const double t_a, const double t_f_gl, const double GammaTS, const double alpha) const {
  // K * K^-1 * m * 1
  std::cout << t_f_gl << std::endl;
  const double l1 = (t_a - t_f_gl) / m_lambda3;
  const double l2 = m_x0 * m_Cd12 * GammaTS  + m_E0 * sin(alpha);
  const double l3 =  m_x0 * (m_Cd12 * GammaTS + m_E0 * sin(alpha));
  return l1 * l2 / l3;
}

//! equation 8 in the PICOP paper.
double PicopPhysics::dimensionless_coordinate(const double z_b, const double z_gl,const  double l) const {
  return (z_b - z_gl) / l;
}

//! equation in Corrigendum of Lazerome et al 2018
double PicopPhysics::dimensionless_melt_curve(const double X_hat) const {
    const std::vector<double> ps = {
        1.371330075095435e-01,   // p0
        5.527656234709359e+01,   // p1
        -8.951812433987858e+02,  // p2
        8.927093637594877e+03,   // p3
        -5.563863123811898e+04,  // p4
        2.218596970948727e+05,   // p5
        -5.820015295669482e+05,  // p6
        1.015475347943186e+06,   // p7
        -1.166290429178556e+06,  // p8
        8.466870335320488e+05,   // p9
        -3.520598035764990e+05,  // p10
        6.387953795485420e+04    // p11
    };

    double result = 0.0;
    for (size_t k = 0; k < ps.size(); ++k) {
        result += ps[k] * pow(X_hat, static_cast<int>(k));
    }

    return result;
}

//! equation 9 in the PICOP paper.
double PicopPhysics::melt_rate(const double M, const double X_hat) const {
  // 1 * m s^-1
  return dimensionless_melt_curve(X_hat) * M;
}

//! equation 10 in the PICOP paper.
double PicopPhysics::melt_function(const double t_a, const double t_f_gl, const double g_alpha) const {
  // m s^-1  degC^-2 * 1 * deg_C^2
  return m_M0 * g_alpha * pow(t_a - t_f_gl, 2);
}

} // end of namespace ocean
} // end of namespace pism
