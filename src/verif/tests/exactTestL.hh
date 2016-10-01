/* Copyright (C) 2016 PISM Authors
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

#ifndef EXACTTESTL_HH
#define EXACTTESTL_HH

#include <vector>
#include <cstddef>
/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exactTestL is a C++ implementation of an isothermal "exact" solution on a
! no-flat bed described in section 2.3 of an incomplete preprint
!
!    Ed Bueler (March 2006) "Equilibrium ice sheets solve variational
!       inequalities"
!
! in this case the exact solution requires solving an ODE numerically.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

struct ExactLParameters {
  ExactLParameters(size_t n);
  std::vector<double> H, a, b;
};

ExactLParameters exactL(const std::vector<double> &r);

#endif /* EXACTTESTL_HH */
