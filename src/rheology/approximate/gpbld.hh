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

#ifndef _GPBLD_APPROXIMATE_H_
#define _GPBLD_APPROXIMATE_H_

#ifdef __cplusplus
extern "C" {
#endif

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

/* Paterson-Budd flow law softness as a function of pressure-adjusted
   temperature. This approximation is valid within the range (200,
   273.15) Kelvin. This code performs a range reduction and then
   approximates the exponent using a polynomial. The relative error is
   around 1e-6. */
double paterson_budd_softness(double T_pa);

/* Glen-Paterson-Budd-Lliboutry-Duval softness as a function of
   enthalpy and pressure. */
double gpbld_softness(double E, double P);

/* Glen-Paterson-Budd-Lliboutry-Duval hardness as a function of
   enthalpy and pressure. */
double gpbld_hardness(double E, double P);

/* Glen-Paterson-Budd-Lliboutry-Duval hardness as a function of
   enthalpy and pressure. */
void gpbld_hardness_n(const double *E, const double *P, unsigned int n, double *result);

/* Glen-Paterson-Budd-Lliboutry-Duval flow function, optimized for the
   Glen exponent n == 3. */
double gpbld_flow(double stress, double E, double P);

/* Glen-Paterson-Budd-Lliboutry-Duval flow function for a column of
   ice, optimized for the Glen exponent n == 3. */
void gpbld_flow_n(const double *stress, const double *E, const double *P,
                  unsigned int n, double *result);

/* Returns constants used by the implementation. */
struct gpbld_constants gpbld_get_constants();

#ifdef __cplusplus
}
#endif

#endif /* _GPBLD_APPROXIMATE_H_ */
