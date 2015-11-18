/* Copyright (C) 2015 PISM Authors
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

#include <math.h>

#include "gpbld.h"

/* include the .c file, not the header */
#include "enthalpy_converter.c"
#include "exp6.h"

static const struct gpbld_constants gpbld = {
  .ideal_gas_constant        = 8.31441,
  .A_cold                    = 3.610e-13,
  .Q_cold                    = 60000.00000,
  .T_critical                = 263.15000,
  .A_warm                    = 1730.00000,
  .Q_warm                    = 139000.00000,
  .water_fraction_coeff      = 181.25,
  .T_min                     = 200,
  .T_melting                 = 273.15,
  .water_frac_observed_limit = 0.01,
};

struct gpbld_constants gpbld_get_constants() {
  return gpbld;
}

/* A Helper function used to compute range reduction constants (cold
   case). */
static double z_cold(double T) {
  return (-gpbld.Q_cold/gpbld.ideal_gas_constant) / T;
}

/* A Helper function used to compute range reduction constants (warm
   case).*/
static double z_warm(double T) {
  return (-gpbld.Q_warm/gpbld.ideal_gas_constant) / T;
}

double paterson_budd_softness(double T_pa) {
  /* these constants will be pre-computed */
  const double
    T_min      = gpbld.T_min,
    T_max      = gpbld.T_melting,
    shift_cold = 0.5 * (z_cold(T_min) + z_cold(gpbld.T_critical)),
    C_cold     = gpbld.A_cold * exp(shift_cold),
    shift_warm = 0.5 * (z_warm(gpbld.T_critical) + z_warm(T_max)),
    C_warm     = gpbld.A_warm * exp(shift_warm);

  double
    QR    = gpbld.Q_cold / gpbld.ideal_gas_constant / 8.0,
    shift = shift_cold / 8.0,
    C     = C_cold;

  if (T_pa >= gpbld.T_critical) {
    QR    = gpbld.Q_warm / gpbld.ideal_gas_constant / 8.0;
    shift = shift_warm / 8.0;
    C     = C_warm;
  }

  double E = exp6(-QR / T_pa - shift);
  /* raise to the power 8 by doubling 3 times */
  E = E * E;
  E = E * E;
  E = E * E;

  return C * E;
}

/* Glen-Paterson-Budd-Lliboutry-Duval softness (temperate case). */
double gpbld_softness_temperate(double E, double E_cts, double P) {
  double omega = enth_water_fraction(E, E_cts, P);
  if (omega > gpbld.water_frac_observed_limit) {
    omega = gpbld.water_frac_observed_limit;
  }

  return paterson_budd_softness(gpbld.T_melting) * (1.0 + gpbld.water_fraction_coeff * omega);
}

double gpbld_softness(double E, double P) {
  const double E_cts = enth_enthalpy_cts(P);
  const double softness_cold = paterson_budd_softness(enth_pressure_adjusted_temperature(E, P));
  const double softness_temp = gpbld_softness_temperate(E, E_cts, P);

  if (E < E_cts) {
    return softness_cold;
  } else {
    return softness_temp;
  }
}

double gpbld_flow(double stress, double E, double P) {
  return gpbld_softness(E, P) * (stress * stress);
}

void gpbld_flow_n(double *restrict stress, double *restrict E, double *restrict P,
                  unsigned int n, double *restrict result) {
  for (unsigned int k = 0; k < n; ++k) {
    result[k] = gpbld_flow(stress[k], E[k], P[k]);
  }
}
