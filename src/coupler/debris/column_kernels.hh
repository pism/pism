// Copyright (C) 2026 Andy Aschwanden and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef PISM_DEBRIS_COLUMN_KERNELS_HH
#define PISM_DEBRIS_COLUMN_KERNELS_HH

#include <vector>

namespace pism {
namespace debris {

/*!
 * Bookkeeping of englacial debris in one ice column.
 *
 * The englacial debris state is stored as *mass per unit area* `m[k]` (kg m^-2) in control
 * volumes around the levels `z[k]` of PISM's vertical grid: volume `k` extends from
 * `zi[k]` to `zi[k+1]`, where `zi[0] = 0`, `zi[k] = (z[k-1] + z[k]) / 2`, and `zi[Mz]` is
 * unbounded. The top volume is truncated at the ice thickness `H`.
 *
 * This is deliberately free of PISM, PETSc and MPI so that it can be unit-tested with a
 * plain executable.
 */
namespace column {

//! Interfaces `zi` (size `z.size() + 1`) of the control volumes around the levels `z`.
std::vector<double> interfaces(const std::vector<double> &z);

//! Thickness of the part of volume `k` that is below the ice surface `H`.
double volume_thickness(const std::vector<double> &zi, int k, double H);

//! Index of the volume containing the height `h` (the top-most volume if `h` is above the
//! last interface).
int volume_index(const std::vector<double> &zi, double h);

//! Convert mass per unit area to concentration (kg m^-3) in a column of height `H`.
void mass_to_concentration(const std::vector<double> &zi, double H,
                           const double *m, double *C);

//! Convert concentration (kg m^-3) to mass per unit area in a column of height `H`.
void concentration_to_mass(const std::vector<double> &zi, double H,
                           const double *C, double *m);

//! Integrated mass per unit area of the column (sum of `m`).
double total_mass(int Mz, const double *m);

/*!
 * Remove the ice between `H_new` and `H_old` (`H_new < H_old`) from the top of the column,
 * assuming uniform concentration within each volume.
 *
 * @return mass per unit area removed
 */
double remove_top(const std::vector<double> &zi, double H_old, double H_new, double *m);

/*!
 * Add ice between `H_old` and `H_new` (`H_new > H_old`) at the top of the column,
 * containing the mass per unit area `M_add` distributed uniformly over the added ice.
 *
 * If `H_new <= H_old` the mass is put in the volume containing `H_old`.
 */
void add_top(const std::vector<double> &zi, double H_old, double H_new, double M_add,
             double *m);

/*!
 * Remove the ice between `0` and `dH` (`dH > 0`) at the bottom of the column and shift the
 * rest of the column down (basal melt).
 *
 * @return mass per unit area removed
 */
double remove_bottom(const std::vector<double> &zi, double H, double dH, double *m);

/*!
 * Move mass located above the ice surface `H` into the top-most volume that is (at least
 * partially) below `H`.
 *
 * @return mass per unit area moved
 */
double fold_above(const std::vector<double> &zi, double H, double *m);

/*!
 * Largest time step allowed by the explicit vertical CFL condition in one column:
 * `min_k dz_k / |w_k|` over the volumes below `H`.
 *
 * @return infinity if the column is at rest
 */
double vertical_dt_max(const std::vector<double> &zi, double H, const double *w);

} // end of namespace column
} // end of namespace debris
} // end of namespace pism

#endif /* PISM_DEBRIS_COLUMN_KERNELS_HH */
