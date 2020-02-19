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

#ifndef GPBLD_N_H
#define GPBLD_N_H

/* Glen-Paterson-Budd-Lliboutry-Duval hardness as a function of
   enthalpy and pressure. */
void gpbld_hardness_n(const double *E, const double *P, unsigned int n, double *result);

/* Glen-Paterson-Budd-Lliboutry-Duval flow function for a column of
   ice, optimized for the Glen exponent n == 3. */
void gpbld_flow_n(const double *stress, const double *E, const double *P,
                  unsigned int n, double *result);

#endif /* GPBLD_N_H */
