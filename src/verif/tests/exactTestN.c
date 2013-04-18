/*
   Copyright (C) 2010 Ed Bueler
  
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
#include "exactTestN.h"

#define secpera  31556926.0    /* seconds per year; 365.2422 days */
#define g        9.81
#define rho      910.0         /* ice density; kg m-3 */
#define rhow     1028.0        /* sea water density; kg m-3 */
#define n        3.0           /* Glen power */

int params_exactN(double *H0, double *L0, double *xc,
                  double *a, double *Hela, double *k,
                  double *H_xc, double *T_xc) {
  double s;

  /* geometry */
  *H0   = 3000.0;
  *L0   = 500.0e3;
  *xc   = 0.9 * (*L0);

  /* mass balance */
  *a    = 0.003 / secpera;   /* s-1; mass balance gradient with elevation */
  *Hela = (*H0) / 1.5;       /* m;  H0 = 1.5 Hela  exactly */

  /* sliding */
  *k    = 9.0 * (*Hela) / ((*a) * (*L0) * (*L0)); /* s m-1; choose k so that eqn (24) gives our L0 */

  /* grounded calving front boundary condition, imposed at xc = .9 L0, determines
     constant vertically-integrated longitudinal stress T; see (2.12) in Schoof (2006);
     treats Hc = H(xc) as exactly at flotation */
  s = (*xc) / (*L0);
  *H_xc = (*H0) * (1.0 - s * s);
  *T_xc = 0.5 * (1.0 - rho / rhow) * rho * g * (*H_xc) * (*H_xc);

  return 0;
}


int exactN(double x, double *H, double *hx, double *u, double *M, double *B, double *beta) {

  double H0, L0, xc, a, Hela, k, Hc, Tc;
  double q, hxx, ux;
  
  params_exactN(&H0, &L0, &xc, &a, &Hela, &k, &Hc, &Tc);

  if (x < 0.0) { return 1; }
  if (x > L0) { return 2; }

  q   = (1.0 / n) - 1.0;           /* a useful power */
  hxx = - 2.0 * H0 / (L0 * L0);    /* constant concavity of h(x) */
  ux  = - hxx / k;                 /* constant strain rate */

  *H = H0 * (1.0 - (x / L0) * (x / L0));  /* eqn (23) in Bodvardsson */

  *hx = hxx * x;
  
  *u = - (*hx) / k;                       /* eqn (10) in Bodvardson, once SIA is dropped */
  
  *M = a * ((*H) - Hela);                 /* page 6 in Bodvardsson, just before eqn (23) */

  *B = Tc / ( 2.0 * (*H) * pow(fabs(ux),q) * ux ); /* Bueler interpretation */

  *beta = k * rho * g * (*H);

  return 0;
}


