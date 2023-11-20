/*
   Copyright (C) 2004-2006, 2014, 2016, 2023 Jed Brown and Ed Bueler

   This file is part of PISM.

   PISM is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 3 of the License, or (at your option) any later
   version.

   PISM is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public License
   along with PISM; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <stdio.h>
#include <math.h>
#include "pism/verification/tests/exactTestH.h"

static const double SperA = 31556926.0;  /* seconds per year; 365.2422 days */

int exactH_old(const double f, const double tIN, const double r,
               double *H, double *M) {

  const double n = 3.0;
  const double H0 = 3600.0, R0=750000.0;
  /* t0 = (beta/Gamma) * pow((2n+1)/((n+1)(1-f)),n) * (pow(R0,n+1)/pow(H0,2n+1))
     when beta=2; */
  double t0 = (15208.0 / pow(1-f,n)) * SperA;
  /* t0 = 40033 years; for test C with isostasy f = rho_ice/rho_rock with
     rho_ice = 910 and rho_rock = 3300 kg/m^3 */
  double lambda, alpha, beta, t0post, Rmargin;
  double t;

  t = tIN;

  if (t < t0) { /* t <= t0: version of test C */
    lambda = 5.0;
    alpha = -1.0;  /* alpha=(2-(n+1)*lambda)/(5*n+3) */
    beta = 2.0;  /* beta=(1+(2*n+1)*lambda)/(5*n+3) */
  } else { /* t >= t0: version of test B */
    t0post = (t0 / 2.0) * (1.0 / 18.0);  /* reset t and t0 */
    t = t - t0 + t0post; /* reset to Halfar w. f */
    t0 = t0post;
    lambda = 0.0;
    alpha = 1.0 / 9.0;  /* alpha=(2-(n+1)*lambda)/(5*n+3)=1/9 */
    beta = 1.0 / 18.0;  /* beta=(1+(2*n+1)*lambda)/(5*n+3)=1/18 */
  }

  Rmargin = R0 * pow(t / t0, beta);
  if (r < Rmargin)
    *H = H0 * pow(t / t0, -alpha)
            * pow(1.0 - pow(pow(t / t0, -beta) * (r / R0), (n + 1.0) / n),
                  n / (2.0 * n + 1.0));
  else
    *H = 0.0;

  if (t > 0.1*SperA)
    *M = (lambda / t) * (*H);
  else {  /* when less than 0.1 year, avoid division by time */
    Rmargin = R0 * pow(0.1*SperA / t0, beta);
    if (r < Rmargin)
      *M = lambda * H0 / t0;  /* constant value in disc of Rmargin radius */
    else
      *M = 0.0;
  }

  return 0;
}

struct TestHParameters exactH(const double f, const double t, const double r) {
  struct TestHParameters result;
  result.error_code = exactH_old(f, t, r, &result.H, &result.M);
  return result;
}
