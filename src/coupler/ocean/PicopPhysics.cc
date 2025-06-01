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

#include "pism/coupler/ocean/PicopPhysics.hh"

#include "pism/util/Config.hh"

namespace pism {
namespace ocean {

PicopPhysics::PicopPhysics(const Config &config) {
  
    m_T_pmp           = config.get_number("constants.fresh_water.melting_point_temperature");

    m_E0             = config.get_number("ocean.picop.entrainment_coefficient");
    m_Cd             = config.get_number("ocean.picop.drag_coefficient");
    m_Cd12GammaT     = config.get_number("ocean.picop.turbulent_heat_exchange_coefficient");
    m_Cd12GammaTS0   = config.get_number("ocean.picop.heat_exchange_parameter");
    m_gamma1         = config.get_number("ocean.picop.heat_exchange_parameter_1");
    m_gamma2         = config.get_number("ocean.picop.heat_exchange_parameter_2");
    m_lambda1        = config.get_number("ocean.picop.freezing_point_salinity_coefficient");
    m_lambda2        = config.get_number("ocean.picop.freezing_point_offset");
    m_lambda3        = config.get_number("ocean.picop.freezing_point_depth_coefficient");
    m_M0             = config.get_number("ocean.picop.melt_rate_parameter");
    m_x0             = config.get_number("ocean.picop.dimensionless_scaling_factor");
  
    m_Cd12 = sqrt(m_Cd);
    m_GammaT = m_Cd12GammaT / m_Cd12;
}

//! equation 4 in the PICOP paper.
double PicopPhysics::characteristic_freezing_poing(double s_a, double z) const {
  // in m/s
  return m_T_pmp + m_lambda1 * s_a + m_lambda2 + m_lambda3 * z;
}

//! equation 5 in the PICOP paper.
double PicopPhysics::effective_heat_exchange_coefficient(double t_a, double t_f_gl, double alpha) const {
  return m_GammaT * (m_gamma1 + m_gamma2 * (t_a - t_f_gl) / m_lambda3 * m_E0 * sin(alpha) / (m_Cd12GammaTS0 + m_E0 * sin(alpha)));
}

//! equation 6 in the PICOP paper.
double PicopPhysics::geometric_scaling(double GammaTS, double alpha) const {
  return pow((sin(alpha) / (m_Cd + m_E0 * sin(alpha))), 0.5) * pow((sin(alpha) / (m_Cd12 * GammaTS * sin(alpha))), 0.5);
}

//! equation 7 in the PICOP paper.
double PicopPhysics::length_scaling(double t_a, double t_f_gl, double GammaTS, double alpha) const {
  return (t_a - t_f_gl) / m_lambda3 * (m_x0 * m_Cd12 * GammaTS  + m_E0 * sin(alpha)) / m_x0 * (m_Cd12 * GammaTS + m_E0 * sin(alpha));
}

//! equation 8 in the PICOP paper.
double PicopPhysics::dimensionless_coordinate(double z_b, double z_gl, double l) const {
  return (z_b - z_gl) / l;
}

//! equation in Corrigendum of Lazerome et al 2018
double PicopPhysics::dimensionless_melt_curve(double X_hat) const {
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
        result += ps[k] * std::pow(X_hat, static_cast<int>(k));
    }

    return result;
}

//! equation 9 in the PICOP paper.
double PicopPhysics::melt_rate(double M, double X_hat) const {
  return dimensionless_melt_curve(X_hat) * M;
}

//! equation 10 in the PICOP paper.
double PicopPhysics::melt_function(double t_a, double s_a, double z_gl, double g_alpha) const {
  return m_M0 * g_alpha * pow(t_a - characteristic_freezing_poing(s_a, z_gl), 2);
}

} // end of namespace ocean
} // end of namespace pism
