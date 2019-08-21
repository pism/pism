/* Copyright (C) 2015, 2016 PISM Authors
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
#include <stdexcept>
#include <gsl/gsl_math.h>       // M_PI

#include "GoldsbyKohlstedt.hh"

namespace pism {
namespace rheology {

// Goldsby-Kohlstedt (forward) ice flow law

GoldsbyKohlstedt::GoldsbyKohlstedt(const std::string &prefix,
                                   const Config &config, EnthalpyConverter::Ptr ec)
  : FlowLaw(prefix, config, ec) {
  m_name = "Goldsby-Kohlstedt / Paterson-Budd (hybrid)";

  m_V_act_vol      = -13.e-6;   // m^3/mol
  m_d_grain_size   = 1.0e-3;    // m  (see p. ?? of G&K paper)
  //--- dislocation creep ---
  m_disl_crit_temp = 258.0;     // Kelvin
  //disl_A_cold    = 4.0e5;                  // MPa^{-4.0} s^{-1}
  //disl_A_warm    = 6.0e28;                 // MPa^{-4.0} s^{-1}
  m_disl_A_cold    = 4.0e-19;   // Pa^{-4.0} s^{-1}
  m_disl_A_warm    = 6.0e4;     // Pa^{-4.0} s^{-1} (GK)
  m_disl_n         = 4.0;       // stress exponent
  m_disl_Q_cold    = 60.e3;     // J/mol Activation energy
  m_disl_Q_warm    = 180.e3;    // J/mol Activation energy (GK)
  //--- grain boundary sliding ---
  m_gbs_crit_temp  = 255.0;     // Kelvin
  //gbs_A_cold     = 3.9e-3;                  // MPa^{-1.8} m^{1.4} s^{-1}
  //gbs_A_warm     = 3.e26;                   // MPa^{-1.8} m^{1.4} s^{-1}
  m_gbs_A_cold     = 6.1811e-14; // Pa^{-1.8} m^{1.4} s^{-1}
  m_gbs_A_warm     = 4.7547e15; // Pa^{-1.8} m^{1.4} s^{-1}
  m_gbs_n          = 1.8;       // stress exponent
  m_gbs_Q_cold     = 49.e3;     // J/mol Activation energy
  m_gbs_Q_warm     = 192.e3;    // J/mol Activation energy
  m_p_grain_sz_exp = 1.4;       // from Peltier
  //--- easy slip (basal) ---
  //basal_A        = 5.5e7;                      // MPa^{-2.4} s^{-1}
  m_basal_A        = 2.1896e-7; // Pa^{-2.4} s^{-1}
  m_basal_n        = 2.4;       // stress exponent
  m_basal_Q        = 60.e3;     // J/mol Activation energy
  //--- diffusional flow ---
  m_diff_crit_temp = 258.0;     // when to use enhancement factor
  m_diff_V_m       = 1.97e-5;   // Molar volume (m^3/mol)
  m_diff_D_0v      = 9.10e-4;   // Preexponential volume diffusion (m^2 second-1)
  m_diff_Q_v       = 59.4e3;    // activation energy, vol. diff. (J/mol)
  m_diff_D_0b      = 5.8e-4;    // preexponential grain boundary coeff.
  m_diff_Q_b       = 49.e3;     // activation energy, g.b. (J/mol)
  m_diff_delta     = 9.04e-10;  // grain boundary width (m)
}

double GoldsbyKohlstedt::flow_impl(double stress, double E,
                                   double pressure, double grainsize) const {
  double temp = m_EC->temperature(E, pressure);
  return flow_from_temp(stress, temp, pressure, grainsize);
}

double GoldsbyKohlstedt::averaged_hardness_impl(double, int,
                                                const double *,
                                                const double *) const {

  throw std::runtime_error("double GoldsbyKohlstedt::averaged_hardness is not implemented");

#ifndef __GNUC__
  return 0;
#endif
}

double GoldsbyKohlstedt::hardness_impl(double enthalpy, double pressure) const {

  // We use the Paterson-Budd relation for the hardness parameter. It would be nice if we didn't
  // have to, but we currently need ice hardness to compute the strain heating. See
  // SIAFD::compute_volumetric_strain_heating().
  double
    T_pa = m_EC->pressure_adjusted_temperature(enthalpy, pressure),
    A = softness_paterson_budd(T_pa);

  return pow(A, m_hardness_power);
}

double GoldsbyKohlstedt::softness_impl(double , double) const {
  throw std::runtime_error("double GoldsbyKohlstedt::softness is not implemented");

#ifndef __GNUC__
  return 0;
#endif
}

/*!
  This is the (forward) Goldsby-Kohlstedt flow law.  See:
  D. L. Goldsby & D. L. Kohlstedt (2001), "Superplastic deformation
  of ice: experimental observations", J. Geophys. Res. 106(M6), 11017-11030.
*/
double GoldsbyKohlstedt::flow_from_temp(double stress, double temp,
                                        double pressure, double gs) const {
  double eps_diff, eps_disl, eps_basal, eps_gbs, diff_D_b;

  if (fabs(stress) < 1e-10) {
    return 0;
  }
  const double T = temp + (m_beta_CC_grad / (m_rho * m_standard_gravity)) * pressure;
  const double pV = pressure * m_V_act_vol;
  const double RT = m_ideal_gas_constant * T;
  // Diffusional Flow
  const double diff_D_v = m_diff_D_0v * exp(-m_diff_Q_v/RT);
  diff_D_b = m_diff_D_0b * exp(-m_diff_Q_b/RT);
  if (T > m_diff_crit_temp) {
    diff_D_b *= 1000; // Coble creep scaling
  }
  eps_diff = 42 * m_diff_V_m *
    (diff_D_v + M_PI * m_diff_delta * diff_D_b / gs) / (RT*(gs*gs));
  // Dislocation Creep
  if (T > m_disl_crit_temp) {
    eps_disl = m_disl_A_warm * pow(stress, m_disl_n-1) * exp(-(m_disl_Q_warm + pV)/RT);
  } else {
    eps_disl = m_disl_A_cold * pow(stress, m_disl_n-1) * exp(-(m_disl_Q_cold + pV)/RT);
  }
  // Basal Slip
  eps_basal = m_basal_A * pow(stress, m_basal_n-1) * exp(-(m_basal_Q + pV)/RT);
  // Grain Boundary Sliding
  if (T > m_gbs_crit_temp) {
    eps_gbs = m_gbs_A_warm * (pow(stress, m_gbs_n-1) / pow(gs, m_p_grain_sz_exp)) *
      exp(-(m_gbs_Q_warm + pV)/RT);
  } else {
    eps_gbs = m_gbs_A_cold * (pow(stress, m_gbs_n-1) / pow(gs, m_p_grain_sz_exp)) *
      exp(-(m_gbs_Q_cold + pV)/RT);
  }

  return eps_diff + eps_disl + (eps_basal * eps_gbs) / (eps_basal + eps_gbs);
}


/*****************
THE NEXT PROCEDURE REPEATS CODE; INTENDED ONLY FOR DEBUGGING
*****************/
GKparts GoldsbyKohlstedt::flowParts(double stress, double temp, double pressure) const {
  double gs, eps_diff, eps_disl, eps_basal, eps_gbs, diff_D_b;
  GKparts p;

  if (fabs(stress) < 1e-10) {
    p.eps_total=0.0;
    p.eps_diff=0.0; p.eps_disl=0.0; p.eps_gbs=0.0; p.eps_basal=0.0;
    return p;
  }
  const double T = temp + (m_beta_CC_grad / (m_rho * m_standard_gravity)) * pressure;
  const double pV = pressure * m_V_act_vol;
  const double RT = m_ideal_gas_constant * T;
  // Diffusional Flow
  const double diff_D_v = m_diff_D_0v * exp(-m_diff_Q_v/RT);
  diff_D_b = m_diff_D_0b * exp(-m_diff_Q_b/RT);
  if (T > m_diff_crit_temp) {
    diff_D_b *= 1000; // Coble creep scaling
  }
  gs = m_d_grain_size;
  eps_diff = 14 * m_diff_V_m *
    (diff_D_v + M_PI * m_diff_delta * diff_D_b / gs) / (RT*(gs*gs));
  // Dislocation Creep
  if (T > m_disl_crit_temp) {
    eps_disl = m_disl_A_warm * pow(stress, m_disl_n-1) * exp(-(m_disl_Q_warm + pV)/RT);
  } else {
    eps_disl = m_disl_A_cold * pow(stress, m_disl_n-1) * exp(-(m_disl_Q_cold + pV)/RT);
  }
  // Basal Slip
  eps_basal = m_basal_A * pow(stress, m_basal_n-1) * exp(-(m_basal_Q + pV)/RT);
  // Grain Boundary Sliding
  if (T > m_gbs_crit_temp) {
    eps_gbs = m_gbs_A_warm * (pow(stress, m_gbs_n-1) / pow(gs, m_p_grain_sz_exp)) *
      exp(-(m_gbs_Q_warm + pV)/RT);
  } else {
    eps_gbs = m_gbs_A_cold * (pow(stress, m_gbs_n-1) / pow(gs, m_p_grain_sz_exp)) *
      exp(-(m_gbs_Q_cold + pV)/RT);
  }

  p.eps_diff=eps_diff;
  p.eps_disl=eps_disl;
  p.eps_basal=eps_basal;
  p.eps_gbs=eps_gbs;
  p.eps_total=eps_diff + eps_disl + (eps_basal * eps_gbs) / (eps_basal + eps_gbs);
  return p;
}
/*****************/

GoldsbyKohlstedtStripped::GoldsbyKohlstedtStripped(const std::string &prefix,
                                                   const Config &config,
                                                   EnthalpyConverter::Ptr ec)
  : GoldsbyKohlstedt(prefix, config, ec) {
  m_name = "Goldsby-Kohlstedt / Paterson-Budd (hybrid, simplified)";

  m_d_grain_size_stripped = 3.0e-3;  // m; = 3mm  (see Peltier et al 2000 paper)
}


double GoldsbyKohlstedtStripped::flow_from_temp(double stress, double temp, double pressure, double) const {
  // note value of gs is ignored
  // note pressure only effects the temperature; the "P V" term is dropped
  // note no diffusional flow
  double eps_disl, eps_basal, eps_gbs;

  if (fabs(stress) < 1e-10) {
    return 0;
  }
  const double T = temp + (m_beta_CC_grad / (m_rho * m_standard_gravity)) * pressure;
  const double RT = m_ideal_gas_constant * T;
  // NO Diffusional Flow
  // Dislocation Creep
  if (T > m_disl_crit_temp) {
    eps_disl = m_disl_A_warm * pow(stress, m_disl_n-1) * exp(-m_disl_Q_warm/RT);
  } else {
    eps_disl = m_disl_A_cold * pow(stress, m_disl_n-1) * exp(-m_disl_Q_cold/RT);
  }
  // Basal Slip
  eps_basal = m_basal_A * pow(stress, m_basal_n-1) * exp(-m_basal_Q/RT);
  // Grain Boundary Sliding
  if (T > m_gbs_crit_temp) {
    eps_gbs = m_gbs_A_warm *
      (pow(stress, m_gbs_n-1) / pow(m_d_grain_size_stripped, m_p_grain_sz_exp)) *
      exp(-m_gbs_Q_warm/RT);
  } else {
    eps_gbs = m_gbs_A_cold *
      (pow(stress, m_gbs_n-1) / pow(m_d_grain_size_stripped, m_p_grain_sz_exp)) *
      exp(-m_gbs_Q_cold/RT);
  }

  return eps_disl + (eps_basal * eps_gbs) / (eps_basal + eps_gbs);
}


} // end of namespace rheology
} // end of namespace pism
