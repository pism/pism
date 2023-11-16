/*
   Copyright (C) 2010, 2014, 2016, 2023 Ed Bueler
  
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

#include <math.h>
#include "pism/verification/tests/exactTestN.h"

#define secpera  31556926.0    /* seconds per year; 365.2422 days */
#define g        9.81
#define rho      910.0         /* ice density; kg m-3 */
#define rhow     1028.0        /* sea water density; kg m-3 */
#define n        3.0           /* Glen power */

struct TestNConstants exactNConstants(void) {
  double s = 0.0;
  struct TestNConstants result;

  /* geometry */
  result.H0   = 3000.0;
  result.L0   = 500.0e3;
  result.xc   = 0.9 * (result.L0);

  /* mass balance */
  result.a     = 0.003 / secpera;   /* s-1; mass balance gradient with elevation */
  result.H_ela = (result.H0) / 1.5;       /* m;  H0 = 1.5 H_ela  exactly */

  /* sliding */
  /* s m-1; choose k so that eqn (24) gives our L0 */
  result.k = 9.0 * (result.H_ela) / ((result.a) * (result.L0) * (result.L0));

  /* grounded calving front boundary condition, imposed at xc = .9 L0, determines
     constant vertically-integrated longitudinal stress T; see (2.12) in Schoof (2006);
     treats Hc = H(xc) as exactly at flotation */
  s = (result.xc) / (result.L0);
  result.H_xc = (result.H0) * (1.0 - s * s);
  result.T_xc = 0.5 * (1.0 - rho / rhow) * rho * g * (result.H_xc) * (result.H_xc);

  return result;
}

struct TestNParameters exactN(double x) {

  double q = 0.0, hxx = 0.0, ux = 0.0;
  const struct TestNConstants c = exactNConstants();
  struct TestNParameters result;
  result.error_code = 0;

  if (x < 0.0) {
    result.error_code = 1;
    return result;
  }

  if (x > c.L0) {
    result.error_code = 2;
    return result;
  }

  q   = (1.0 / n) - 1.0;              /* a useful power */
  hxx = - 2.0 * c.H0 / (c.L0 * c.L0); /* constant concavity of h(x) */
  ux  = - hxx / c.k;                  /* constant strain rate */

  result.H = c.H0 * (1.0 - (x / c.L0) * (x / c.L0));  /* eqn (23) in Bodvardsson */

  result.h_x = hxx * x;
  
  result.u = - (result.h_x) / c.k; /* eqn (10) in Bodvardson, once SIA is dropped */
  
  result.M = c.a * ((result.H) - c.H_ela); /* page 6 in Bodvardsson, just before eqn (23) */

  result.B = c.T_xc / (2.0 * (result.H) * pow(fabs(ux),q) * ux); /* Bueler interpretation */

  result.beta = c.k * rho * g * (result.H);

  return result;
}


