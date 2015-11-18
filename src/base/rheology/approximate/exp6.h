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

#ifndef _EXP_APPROXIMATION_H_
#define _EXP_APPROXIMATION_H_

/* A polynomial approximation of f(x) = exp(x) in the range (-0.6,
   0.6). Coefficients of this 6-degree polynomial were computed using
   NumPy polyfit with NumPy exp() evaluated at Chebyshev nodes for
   (-0.6, 0.6). This approximation has the relative error of about
   1e-6. */
double exp6(double x) {
  static const unsigned int p_len = 7;
  static const double p[] = {
    0.00140459224886816855769333667325327041908167302608489990234375,
    0.00845911740260831858384538151085507706739008426666259765625000,
    0.04166383582482915959310787457070546224713325500488281250000000,
    0.16664398312747433217317905018717283383011817932128906250000000,
    0.50000012745147282000601762774749659001827239990234375000000000,
    1.00000102139571334092238430457655340433120727539062500000000000,
    0.99999999999999922284388276239042170345783233642578125000000000,
  };
  double result = p[0];
  for (unsigned int j = 1; j < p_len; ++j) {
#ifdef FP_FAST_FMA
    result = fma(result, x, p[j]);
#else
    result = result * x + p[j];
#endif
  }
  return result;
}

#endif /* _EXP_APPROXIMATION_H_ */
