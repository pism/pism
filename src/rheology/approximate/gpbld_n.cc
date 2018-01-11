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

#include "gpbld.hh"

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

void gpbld_flow_n(const double *stress, const double *E, const double *P,
                  unsigned int n, double *result) {
#pragma ivdep
  for (unsigned int k = 0; k < n; ++k) {
    result[k] = gpbld_flow_inlined(stress[k], E[k], P[k]);
  }
}
