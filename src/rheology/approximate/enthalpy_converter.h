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

#ifndef _C_ENTHALPY_CONVERTER_H_
#define _C_ENTHALPY_CONVERTER_H_

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

static const struct enth_constants enth = {
  /* .T_melting = */ 273.15,
  /* .c_i = */ 2009.0,
  /* .c_w = */ 4170,
  /* .T_0 = */ 223.15,
  /* .beta = */ 7.9e-8,
  /* .L0 = */ 3.34e5
};

/* Melting temperature as a function of pressure. */
inline double enth_melting_temperature(double P) {
  return enth.T_melting - enth.beta * P;
}

/* Temperature as a function of enthalpy, cold case. */
inline double enth_temp_cold(double E) {
  return (1.0 / enth.c_i) * E + enth.T_0;
}

/* Pressure-adjusted temperature. */
inline double enth_pressure_adjusted_temperature(double E, double P) {
  return enth_temp_cold(E) - enth_melting_temperature(P) + enth.T_melting;
}

/* Cold-temperate transition enthalpy. */
inline double enth_enthalpy_cts(double P) {
  return enth.c_i * (enth_melting_temperature(P) - enth.T_0);
}

/* Latent heat of fusion as a function of melting temperature. */
inline double enth_L(double T_pm) {
  return enth.L0 + (enth.c_w - enth.c_i) * (T_pm - enth.T_melting);
}

/* Liquid water fraction as a function of enthalpy, enthalpy at the
   cold-temperate transition, and pressure.*/
inline double enth_water_fraction(double E, double E_cts, double P) {
  return (E - E_cts) / enth_L(enth_melting_temperature(P));
}

/* Returns constants used by the enthalpy converter code. */
inline struct enth_constants enth_get_constants() {
  return enth;
}

#endif /* _C_ENTHALPY_CONVERTER_H_ */
