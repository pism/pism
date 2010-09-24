/*
   Copyright (C) 2010 Ed Bueler
  
   This file is part of PISM.
  
   PISM is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2 of the License, or (at your option) any later
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
#include "exactTestN.h"

#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a > b ? b : a)

#define pi       3.1415926535897931
#define secpera  31556926.0    /* seconds per year; 365.2422 days */
#define g        9.81
#define rho      910.0         /* ice density; kg/m^3 */
#define n        3.0           /* Glen power */
#define barB     1.9e8         /* vertical average hardness; Pa s^(1/3);
                                  from MacAyeal et al 1996 */

#define L        300.0e3       /* m */
#define b0       1000.0        /* m */
#define h0       3500.0        /* m */
#define umax     800.0/secpera /* m s-1 */


int exactN(double t, double x,
           double *u, double *ux, double *h, double *b, double *hx,
           double *taud, double *taub, double *tau11,
           double *C) {

  double H, ulow;
  double AA, BB, aaa, bbb;  /* factors in h and in hx */
  double tau11x; /* longitudinal stress gradient */
  const double s = x / L;

  if (t < 0.0) { return 1; }
  if (fabs(x) > 400.0e3) { return 2; }

  *u = umax * pow(fabs(s), 5.0) * s;

  *b = b0 * cos(0.6 * pi * x / L);
  
  AA = h0 + 10.0 * sin(pi * fabs(s)) * (t / secpera);
  BB = pow(1.0 - 0.99 * fabs(s),0.5);
  *h = AA * BB;

  H = (*h) - (*b);
  *ux = (umax / L) * pow(fabs(s), 5.0);
  *tau11 = 2.0 * barB * H * pow(fabs(*ux),(1.0/n)-1.0) * (*ux);

  aaa = (10.0 * pi / L) * cos(pi * fabs(s)) * (t / secpera);
  bbb = 0.5 * pow(1.0 - 0.99 * fabs(s),-0.5) * (- 0.99 / L);
  *hx = aaa * BB + AA * bbb;
  if (x < 0.0) {
    *hx = - (*hx);
  } else if (fabs(x) < 1.0) {
    *hx = 0.0;
  }

  *taud = rho * g * H * (*hx);

/* FIXME */
  tau11x = 0.0;
  
  ulow = 0.001 / secpera;
  if (fabs(*u) > ulow) {
    *C = ( tau11x - *taud ) / (*u);
  } else {
    *C = ( tau11x - *taud ) / ulow;
  }
  *taub = - (*C) * (*u);

  return 0;
}


