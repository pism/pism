/* Copyright (C) 2025, 2026 Andy Aschwanden, Constantine Khrulev
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

    beta_S         = config.get_number("ocean.picop.haline_concentration_coefficient");
    beta_T         = config.get_number("ocean.picop.thermal_expansion_coefficient");
    E0             = config.get_number("ocean.picop.entrainment_coefficient");
    c_p            = config.get_number("constants.fresh_water.specific_heat_capacity");
    Cd             = config.get_number("ocean.picop.drag_coefficient");
    CdT            = config.get_number("ocean.picop.turbulent_heat_exchange_coefficient");
    CdTS0          = config.get_number("ocean.picop.heat_exchange_parameter");
    g              = config.get_number("constants.standard_gravity");
    gamma1         = config.get_number("ocean.picop.heat_exchange_parameter_1");
    gamma2         = config.get_number("ocean.picop.heat_exchange_parameter_2");
    lambda1        = config.get_number("ocean.picop.freezing_point_salinity_coefficient");
    lambda2        = config.get_number("ocean.picop.freezing_point_offset", "kelvin");
    lambda3        = config.get_number("ocean.picop.freezing_point_depth_coefficient");
    L_fw           = config.get_number("constants.fresh_water.latent_heat_of_fusion");
    M0             = config.get_number("ocean.picop.melt_rate_parameter", "m s^-1 kelvin^-2");
    x0             = config.get_number("ocean.picop.dimensionless_scaling_factor");
    power_alpha    = config.get_number("ocean.picop.power_alpha");
    power_beta     = config.get_number("ocean.picop.power_beta");
    YT = CdT / sqrt(Cd);
}


//! equation 4 in the PICOP paper.
double PicopPhysics::characteristic_freezing_point(const double s_a, const double z) const {
  // in K * g /kg + K  + K / m
  return lambda1 * s_a + lambda2 + lambda3 * z;
}

//! equation 2 in the PICOP/q_sg paper.
double PicopPhysics::freezing_point_depth(const double t_a, const double s_a, const double z) const {
  // in m
  return (t_a - lambda1 * s_a + lambda2 + lambda3 * z) / lambda3;
}

//! equation 5 in the PICOP paper.
double PicopPhysics::effective_heat_exchange_coefficient(const double t_a, const double t_f_gl, const double alpha) const {

  return YT * (gamma1 + gamma2 * (((t_a-t_f_gl) * E0 * sin(alpha)) / (lambda3 * (CdTS0 + E0 * sin(alpha)))));
}

double PicopPhysics::geometric_scaling(const double Gamma_TS, const double alpha) const {
  const double CdTS = sqrt(Cd) * Gamma_TS;
  const double E0_sin_alpha = E0 * sin(alpha);
  const double G1 = sqrt(sin(alpha) / (Cd + E0_sin_alpha));
  const double G2 = sqrt(CdTS / (CdTS + E0_sin_alpha));
  const double G3 = E0_sin_alpha / (CdTS + E0_sin_alpha);
  return G1 * G2 * G3;
}

//! equation 7 in the PICOP paper.
double PicopPhysics::length_scaling(const double t_a, const double t_f_gl, const double Gamma_TS, const double alpha) const {
  const double CdTS = sqrt(Cd) * Gamma_TS;
  const double L1 = (t_a - t_f_gl) / lambda3;
  const double L2 = x0 * CdTS  + E0 * sin(alpha);
  const double L3 =  x0 * (CdTS + E0 * sin(alpha));
  return L1 * L2 / L3;
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
  // M = M0 * g(alpha) * (T_a - T_f)^beta. The thermal-forcing exponent beta (Eq. 10 uses 2,
  // the Antarctic plume value) is configurable via ocean.picop.power_beta; Cai et al. (2017)
  // find beta ~ 1.2 for Petermann. Clamp the base at 0 so a non-integer exponent is safe (and
  // gives no melt when the ambient water is at/below the local freezing point).
  return M0 * g_alpha * pow(std::max(0.0, t_a - t_f_gl), power_beta);
}

//! Equation 13 in Pelle et al (2023)
double PicopPhysics::fresh_water_melt_rate(const double q_sg,
                                           const double s_a,
                                           const double t_a,
                                           const double Gamma_TS,
                                           const double t_f_gl,
                                           const double alpha) const {
  // No subglacial discharge (or a non-finite discharge value) => no discharge-driven
  // plume, hence no discharge melt. Written as "not (> 0)" so a NaN q_sg returns 0 too.
  if (not (q_sg > 0.0)) {
    return 0.0;
  }

  // thermal forcing lambda3 * dTf_gl (K) -- same convention as melt_function()
  const double thermal_forcing = t_a - t_f_gl;

  // delta_rho_i = beta_S * S_a - beta_T * (lambda3 * dTf_gl), with S = 0 (fresh discharge).
  // The discharge plume only upwells (and melts) while it is buoyant, i.e. delta_rho_i > 0.
  const double delta_rho_i = beta_S * s_a - beta_T * thermal_forcing;
  if (not (delta_rho_i > 0.0)) {
    return 0.0;
  }

  const double CdTS = sqrt(Cd) * Gamma_TS;
  const double E0_sin_alpha = E0 * sin(alpha);
  const double G1 = sqrt(sin(alpha) / (Cd + E0_sin_alpha));
  const double G2 = sqrt(CdTS / (CdTS + E0_sin_alpha));

  // Pelle et al. (2023), Eq. 13 (in m s^-1; the paper's 31536000 m/s->m/yr factor is
  // omitted since PISM works in SI). The prefactor is c_p / L_fw = 1 / (L / c_p).
  //
  // The 1/3 power in Eqs. 13-14 is the plume-velocity exponent; it applies to the whole
  // buoyancy-and-geometry group G2 * (g q_sg delta_rho), which the paper (and Eq. 14) keeps
  // under a single cube root. It is configurable via ocean.picop.power_alpha (default 1/3, at
  // which pow(G2 * g q_sg delta_rho, 1/3) == cbrt(G2) * cbrt(g q_sg delta_rho), the original).
  //
  // pow() is safe: G2 = sqrt(...) >= 0, g > 0, and q_sg > 0 and delta_rho_i > 0 are guaranteed
  // above, so the base is non-negative (no NaN as pow() would give for a negative base).
  return (c_p / L_fw) * CdTS0 * G1
         * pow(G2 * g * q_sg * delta_rho_i, power_alpha) * thermal_forcing;
}

//! Governing length scale 5L' (m), Pelle et al. (2023), Eq. 15.
/*!
 * Eqs. 13 and 15 share the same geometric/thermal block and are reciprocal in it, so
 * they collapse to  5L' = 5 * q_sg / m_fw  (the paper's leading 500 = 5 * L_fw/c_p is
 * carried by 1/m_fw, whose prefactor is c_p/L_fw). Reusing m_fw guarantees consistency
 * with Eq. 13 and keeps the full slope dependence (no "sin(alpha) ~ 0.01" assumption).
 *
 * @param[in] q_sg subglacial discharge flux at the outflow (m^2 s^-1)
 * @param[in] m_fw discharge melt rate from fresh_water_melt_rate() (m s^-1)
 */
double PicopPhysics::governing_length_scale(const double q_sg, const double m_fw) const {
  if (m_fw <= 0.0) {
    return 0.0;
  }
  return 5.0 * q_sg / m_fw;
}

} // end of namespace ocean
} // end of namespace pism
