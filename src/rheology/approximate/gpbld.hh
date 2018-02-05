/* Copyright (C) 2015-2017 PISM Authors
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

#include "pism/external/vdt/vdtMath.h" // fast_exp

#include "enthalpy_converter.h"

#include "inverse_cbrt.h"       // inverse_cbrt

#ifndef _GPBLD_APPROXIMATE_H_
#define _GPBLD_APPROXIMATE_H_

/* flow law constants */
struct gpbld_constants {
  /* Ideal gas constant, Joule / (mol Kelvin) */
  double ideal_gas_constant;
  /* Paterson-Budd cold case activation energy, Pascal-3 second-1 */
  double A_cold;
  /* A_cold^(-1/3) */
  double A_cold_inv_cuberoot;
  /* Paterson-Budd Q_cold, Joule / mol */
  double Q_cold;
  /* Paterson-Budd critical temperature, Kelvin */
  double T_critical;
  /* Paterson-Budd warm case activation energy, Pascal-3 second-1 */
  double A_warm;
  /* A_warm^(-1/3) */
  double A_warm_inv_cuberoot;
  /* Paterson-Budd Q_warm, Joule / mol */
  double Q_warm;
  /* Glen-Paterson-Budd-Lliboutry-Duval softness parameter, pure number */
  double water_fraction_coeff;
  /* Melting point of pure water, Kelvin */
  double T_melting;
  /* Maximum observed liquid water fraction. */
  double water_frac_observed_limit;
};

static const struct gpbld_constants gpbld = {
  /* .ideal_gas_constant        = */ 8.31441,
  /* .A_cold                    = */ 3.610e-13,
  /* .A_cold_inv_cuberoot       = */ 14044.2192037997098,
  /* .Q_cold                    = */ 60000.00000,
  /* .T_critical                = */ 263.15000,
  /* .A_warm                    = */ 1730.00000,
  /* .A_warm_inv_cuberoot       = */ 0.083301207912519429,
  /* .Q_warm                    = */ 139000.00000,
  /* .water_fraction_coeff      = */ 181.25,
  /* .T_melting                 = */ 273.15,
  /* .water_frac_observed_limit = */ 0.01,
};

/* Returns constants used by the implementation. */
inline struct gpbld_constants gpbld_get_constants() {
  return gpbld;
}

/* Paterson-Budd flow law softness as a function of pressure-adjusted
   temperature. This approximation is valid within the range (200,
   273.15) Kelvin. This code performs a range reduction and then
   approximates the exponent using a polynomial. The relative error is
   around 1e-6. */
inline double paterson_budd_softness(double T_pa) {
  double
    QR = gpbld.Q_cold / gpbld.ideal_gas_constant,
    A  = gpbld.A_cold;

  if (T_pa >= gpbld.T_critical) {
    QR = gpbld.Q_warm / gpbld.ideal_gas_constant;
    A  = gpbld.A_warm;
  }

  return A * vdt::fast_exp(-QR / T_pa);
}

/* Glen-Paterson-Budd-Lliboutry-Duval softness (temperate case). */
inline double gpbld_softness_temperate(double E, double E_cts, double P) {
  double omega = enth_water_fraction(E, E_cts, P);
  if (omega > gpbld.water_frac_observed_limit) {
    omega = gpbld.water_frac_observed_limit;
  }

  return paterson_budd_softness(gpbld.T_melting) * (1.0 + gpbld.water_fraction_coeff * omega);
}

/* Glen-Paterson-Budd-Lliboutry-Duval softness as a function of
   enthalpy and pressure. */
inline double gpbld_softness(double E, double P) {
  const double E_cts = enth_enthalpy_cts(P);
  const double softness_cold = paterson_budd_softness(enth_pressure_adjusted_temperature(E, P));
  const double softness_temp = gpbld_softness_temperate(E, E_cts, P);

  if (E < E_cts) {
    return softness_cold;
  } else {
    return softness_temp;
  }
}

inline double paterson_budd_hardness(double T_pa) {

  double
    QR = gpbld.Q_cold / gpbld.ideal_gas_constant / 3.0,
    A  = gpbld.A_cold_inv_cuberoot;

  if (T_pa >= gpbld.T_critical) {
    QR = gpbld.Q_warm / gpbld.ideal_gas_constant / 3.0;
    A  = gpbld.A_warm_inv_cuberoot;
  }

  return A * vdt::fast_exp(QR / T_pa);
}

inline double gpbld_hardness_temperate(double E, double E_cts, double P) {
  double omega = enth_water_fraction(E, E_cts, P);
  if (omega > gpbld.water_frac_observed_limit) {
    omega = gpbld.water_frac_observed_limit;
  }

  double C = inverse_cbrt(1.0 + gpbld.water_fraction_coeff * omega);

  return paterson_budd_hardness(gpbld.T_melting) * C;
}

/* Glen-Paterson-Budd-Lliboutry-Duval hardness as a function of
   enthalpy and pressure. */
inline double gpbld_hardness(double E, double P) {
  const double E_cts = enth_enthalpy_cts(P);
  const double hardness_cold = paterson_budd_hardness(enth_pressure_adjusted_temperature(E, P));
  const double hardness_temp = gpbld_hardness_temperate(E, E_cts, P);

  if (E < E_cts) {
    return hardness_cold;
  } else {
    return hardness_temp;
  }
}

/* Glen-Paterson-Budd-Lliboutry-Duval flow function, optimized for the
   Glen exponent n == 3. */
inline double gpbld_flow(double stress, double E, double P) {
  return gpbld_softness(E, P) * (stress * stress);
}

// Same as gpbld_flow, but with all the function calls (except for vdt::fast_exp())
// inlined to help less capable optimizing compilers.
inline double gpbld_flow_inlined(double stress, double E, double P) {
  double softness = 0.0;
  {
    const double
      T_melting = enth.T_melting - enth.beta * P,
      T_pa      = ((1.0 / enth.c_i) * E + enth.T_0) - T_melting + enth.T_melting,
      E_cts     = enth.c_i * (T_melting - enth.T_0);

    double softness_cold = 0.0;
    {
      double
        QR = gpbld.Q_cold / gpbld.ideal_gas_constant,
        A  = gpbld.A_cold;

      if (T_pa >= gpbld.T_critical) {
        QR = gpbld.Q_warm / gpbld.ideal_gas_constant;
        A  = gpbld.A_warm;
      }

      softness_cold = A * vdt::fast_exp(-QR / T_pa);
    }

    double softness_temp = 0.0;
    {
      const double L = enth.L0 + (enth.c_w - enth.c_i) * (T_melting - enth.T_melting);

      double omega = (E - E_cts) / L;

      if (omega > gpbld.water_frac_observed_limit) {
        omega = gpbld.water_frac_observed_limit;
      }

      return softness_cold * (1.0 + gpbld.water_fraction_coeff * omega);
    }

    if (E < E_cts) {
      softness = softness_cold;
    } else {
      softness = softness_temp;
    }
  }

  return softness * (stress * stress);
}


#endif /* _GPBLD_APPROXIMATE_H_ */
