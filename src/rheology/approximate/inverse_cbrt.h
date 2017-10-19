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

#ifndef _INVERSE_CBRT_H_
#define _INVERSE_CBRT_H_

/* Evaluate a polynomial of degree n-1 with coefficients p[j] at a
 * point x. Uses Horner's rule and assumes that coefficients are
 * stored with the highest degree first.
 */
static inline double polyval(const double *p, unsigned int n, double x) {
  double result = p[0];
  for (unsigned int j = 1; j < n; ++j) {
    result = result * x + p[j];
  }
  return result;
}

static const unsigned int p_inverse_cbrt_n = 7;
static const double p_inverse_cbrt[] = {
   2.053062934487385072374054928445730183739215135574340820312500e-03,
  -2.920177205984430396368267679463315289467573165893554687500000e-02,
   1.756105120169964839416110180536634288728237152099609375000000e-01,
  -5.816529863373526287873005458095576614141464233398437500000000e-01,
   1.162814102843798158559707189851906150579452514648437500000000e+00,
  -1.481389224164499474056810868205502629280090332031250000000000e+00,
   1.751738939472610612213543390680570155382156372070312500000000e+00};

/* Approximate x^(-1/3) on the interval [1, 3] using a 6-th degree
 * polynomial approximation followed by two iterations of the Newton's
 * method. The relative error of the result is less than twice the
 * machine epsilon.
 *
 * Note that range reduction is not needed here.
 *
 * This code can be vectorized, unlike a call to pow(x, -1.0/3.0).
 */
inline double inverse_cbrt(double z) {
  static const double one_third = 1.0 / 3.0;
  double z_over_3 = one_third * z;

  double x = polyval(p_inverse_cbrt, p_inverse_cbrt_n, z);
  /* Group multiplications so that (x * x) and (x * z_over_3) can be
     computed in parallel. */
  x -= x * ((x * x) * (x * z_over_3) - one_third);
  x -= x * ((x * x) * (x * z_over_3) - one_third);

  return x;
}

#endif /* _INVERSE_CBRT_H_ */
