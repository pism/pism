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

#ifndef _C_ENTHALPY_CONVERTER_H_
#define _C_ENTHALPY_CONVERTER_H_

#ifdef __cplusplus
extern "C" {
#endif

struct enth_constants {
  /* Melting temperature of water, Kelvin */
  double T_melting;
  /* Specific heat capacity of ice, Joule / (kg Kelvin) */
  double c_i;
  /* Specific heat capacity of water, Joule / (kg Kelvin) */
  double c_w;
  /* Enthalpy reference temperature, Kelvin */
  double T_0;
  /* Clausius-Clapeyron constant, Kelvin / Pascal */
  double beta;
  /* Latent heat of fusion of water, Joule / kg */
  double L0;
};

/* Melting temperature as a function of pressure. */
double enth_melting_temperature(double P);

/* Temperature as a function of enthalpy, cold case. */
double enth_temp_cold(double E);

/* Pressure-adjusted temperature. */
double enth_pressure_adjusted_temperature(double E, double P);

/* Cold-temperate transition enthalpy. */
double enth_enthalpy_cts(double P);

/* Latent heat of fusion as a function of melting temperature. */
double enth_L(double T_pm);

/* Liquid water fraction as a function of enthalpy, enthalpy at the
   cold-temperate transition, and pressure.*/
double enth_water_fraction(double E, double E_cts, double P);

/* Returns constants used by the enthalpy converter code. */
struct enth_constants enth_get_constants();

#ifdef __cplusplus
}
#endif

#endif /* _C_ENTHALPY_CONVERTER_H_ */
