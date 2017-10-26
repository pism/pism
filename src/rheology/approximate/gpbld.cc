/* Copyright (C) 2015, 2017 PISM Authors
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

#include "gpbld.hh"

#include "pism/external/vdt/vdtMath.h"        // fast_exp

#include "inverse_cbrt.h"       // inverse_cbrt

/* include the .c file, not the header */
#include "enthalpy_converter.c"

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

struct gpbld_constants gpbld_get_constants() {
  return gpbld;
}

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

void paterson_budd_hardness_n(double *T_pa,
                              unsigned int n, double *result) {
#pragma ivdep
  for (unsigned int k = 0; k < n; ++k) {
    result[k] = paterson_budd_hardness(T_pa[k]);
  }
}


/* Glen-Paterson-Budd-Lliboutry-Duval softness (temperate case). */
inline double gpbld_softness_temperate(double E, double E_cts, double P) {
  double omega = enth_water_fraction(E, E_cts, P);
  if (omega > gpbld.water_frac_observed_limit) {
    omega = gpbld.water_frac_observed_limit;
  }

  return paterson_budd_softness(gpbld.T_melting) * (1.0 + gpbld.water_fraction_coeff * omega);
}

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

inline double gpbld_hardness_temperate(double E, double E_cts, double P) {
  double omega = enth_water_fraction(E, E_cts, P);
  if (omega > gpbld.water_frac_observed_limit) {
    omega = gpbld.water_frac_observed_limit;
  }

  double C = inverse_cbrt(1.0 + gpbld.water_fraction_coeff * omega);

  return paterson_budd_hardness(gpbld.T_melting) * C;
}

double gpbld_hardness(double E, double P) {
  const double E_cts = enth_enthalpy_cts(P);
  const double hardness_cold = paterson_budd_hardness(enth_pressure_adjusted_temperature(E, P));
  const double hardness_temp = gpbld_hardness_temperate(E, E_cts, P);

  if (E < E_cts) {
    return hardness_cold;
  } else {
    return hardness_temp;
  }
}

void gpbld_hardness_n(const double *E, const double *P,
                      unsigned int n, double *result) {
#pragma ivdep
  for (unsigned int k = 0; k < n; ++k) {
    // we could call gpbld_hardness(E[k], P[k]), but clang thinks that inlining it is too costly
    const double E_cts = enth_enthalpy_cts(P[k]);
    const double hardness_cold = paterson_budd_hardness(enth_pressure_adjusted_temperature(E[k], P[k]));
    const double hardness_temp = gpbld_hardness_temperate(E[k], E_cts, P[k]);

    if (E[k] < E_cts) {
      result[k] = hardness_cold;
    } else {
      result[k] = hardness_temp;
    }
  }
}

inline double gpbld_flow(double stress, double E, double P) {
  return gpbld_softness(E, P) * (stress * stress);
}

void gpbld_flow_n(const double *stress, const double *E, const double *P,
                  unsigned int n, double *result) {
#pragma ivdep
  for (unsigned int k = 0; k < n; ++k) {
    result[k] = gpbld_flow(stress[k], E[k], P[k]);
  }
}
