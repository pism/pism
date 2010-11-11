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

#define secpera  31556926.0    /* seconds per year; 365.2422 days */
#define g        9.81
#define rho      910.0         /* ice density; kg/m^3 */
#define n        3.0           /* Glen power */

int geometry_exactN(double *H0, double *L0, double *xc) {
  *H0 = 3000.0;
  *L0 = 500.0e3;
  *xc = 0.9 * (*L0);
  return 0;
}

#define Hela     (H0 / 1.5)                   /* m;  H0 = 1.5 Hela  exactly */
#define hxx      (- 2.0 * H0 / (L0 * L0))     /* constant concavity of h(x) */
#define q        ((1.0 / n) - 1.0)            /* a useful power */

#define a        (0.001 / secpera)            /* s-1; mass balance gradient with elevation */
#define k        (9.0 * Hela / (a * L0 * L0)) /* s m-1; choose k so that eqn (24) gives our L0 */
#define ux       (- hxx / k)                  /* constant strain rate */

/* grounded calving front boundary condition, imposed at xc = .9 L0, determines
   constant vertically-integrated longitudinal stress T; see (2.12) in Schoof (2006);
   treats Hc = H(xc) as exactly at flotation */
#define rhow     1028.0

#define Hc       (H0 * (1.0 - (xc / L0) * (xc / L0)))
#define Ttau     (0.5 * (1.0 - rho / rhow) * rho * g * Hc * Hc)

int exactN(double x, double *h, double *hx, double *u, double *M, double *A) {

  double H0, L0, xc;
  
  geometry_exactN(&H0, &L0, &xc);

  if (x < 0.0) { return 1; }
  if (x > L0) { return 2; }

  *h = H0 * (1.0 - (x / L0) * (x / L0));  /* eqn (23) in Bodvardsson */

  *hx = hxx * x;
  
  *u = - (*hx) / k;                       /* eqn (10) in Bodvardson, once SIA is dropped */
  
  *M = a * ((*h) - Hela);                 /* page 6 in Bodvardsson, just before eqn (23) */

  *A = pow(2.0 * (*h) * pow(fabs(ux),q) * ux / Ttau, n); /* Bueler interpretation */

  return 0;
}


