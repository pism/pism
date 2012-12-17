/*
   Copyright (C) 2012 Ed Bueler

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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include "exactTestP.h"

#define SperA    31556926.0    /* seconds per year; 365.2422 days */
#define g        9.81          /* m s-2 */
#define rhoi     910.0         /* kg m-3 */
#define rhow     1000.0        /* kg m-3 */

/* major model parameters: */
#define Aglen    3.1689e-24    /* Pa-3 s-1 */
#define K        1.0e-2        /* m s-1 */
#define Wr       1.0           /* m */
#define c1       0.500         /* m-1 */
#define c2       0.040         /* [pure] */

/* model regularizations */
#define E0       1.0           /* m */
#define Y0       0.001         /* m */

/* specific to exact solution */
#define Phi0     0.020         /* m a-1 */
#define h0       500.0         /* m */
#define v0       100.0 / SperA /* m s-1 */
#define R1       5000.0        /* m */


int funcP(double r, const double W[], double f[], void *params) {
  /* Computes RHS f(r,W) for differential equation:
      dW
      -- = f(r,W) = FIXME
      dr
  Compare doublediff.m at https://github.com/bueler/hydrolakes
  Assumes Glen power n=3.
  */

  double sb, dsb, CC, CD, CZ, zz, dPo, tmp1, c0, vphi0, numer, denom;

  if (params == NULL) {} /* quash warning "unused parameters" */

  if ((r >= 0.0) && (r <= L)) {
    if (r < R1) {
      sb  = 0.0;
      dsb = 0.0;
    } else {
      CC  = c1 / (c2 * Aglen);
      /* vb = v0 * (r - R1).^5 / (R0-R1)^5   and   sb = (CC * vb)^(1/3) */
      CZ  = pow(CC * v0, 1.0/3.0);
      zz  = pow((r - R1) / (R0 - R1), 1.0/3.0);
      sb  = CZ * pow(zz,5.0);
      CD  = (5.0 * CZ) / (3.0 * (R0 - R1));
      dsb = CD * zz * zz;
    }
    dPo   = - (2.0 * rhoi * g * h0 / (R0*R0)) * r;
    numer = dsb * (W[0] + Y0) * (Wr - W[0]);
    tmp1  = pow(W[0] + Y0,4.0/3.0) * pow(Wr - W[0],2.0/3.0);
    c0    = K / (rhow * g);
    vphi0 = Phi0 / (2 * c0);
    numer = numer - ( vphi0 * r / W[0] + dPo) * tmp1;
    denom = (1.0/3.0) * (Wr + Y0) * sb + rhow * g * tmp1;
    f[0] = numer / denom;
    return GSL_SUCCESS;
  } else {
    f[0] = 0.0;  /* place-holder */
    return TESTP_R_OUT_OF_RANGE;
  }
}


double initialconditionW() {
  /* in notes: return value is W_c(L^-) */
  double hL, vbL, PoL, sbL;
  hL  = h0 * (1.0 - (L/R0) * (L/R0));
  vbL = v0 * pow( (L - R1) / (R0-R1) ,5.0);
  PoL = rhoi * g * hL;
  sbL = pow( c1 * vbL / (c2 * Aglen) ,1.0/3.0);
  return (pow(sbL,3.0) * Wr - pow(PoL,3.0) * Y0) / (pow(sbL,3.0) + pow(PoL,3.0));
}


/* combination EPS_ABS = 1e-12, EPS_REL=0.0, method = 1 = RK Cash-Karp
   is believed to be predictable and accurate */
int getW(double *r, int N, double *W,
         const double EPS_ABS, const double EPS_REL, const int ode_method) {
   /* solves ODE for W(r), the exact soluiton
      r and W must be allocated vectors of length N; r[] must be decreasing */

   int i, count;
   const gsl_odeiv_step_type* T;
   int status = TESTP_NOT_DONE;
   double rr, step;

   gsl_odeiv_step*    s;
   gsl_odeiv_control* c;
   gsl_odeiv_evolve*  e;
   gsl_odeiv_system   sys = {funcP, NULL, 1, NULL};  /* Jac-free method and no params */

   /* check first: we have a list, r is decreasing, r is in range [0,R0] */
   if (N < 1) return TESTP_NO_LIST;
   for (i = 1; i<N; i++) {
     if (r[i] > r[i-1]) return TESTP_NOT_DECREASING;
     if (r[i] < 0.0)    return TESTP_R_OUT_OF_RANGE;
     if (r[i] > L)     return TESTP_R_OUT_OF_RANGE;
   }

   /* setup for GSL ODE solver; following step choices don't need Jacobian,
      but should we chose one that does?  */
   switch (ode_method) {
     case 1:
       T = gsl_odeiv_step_rkck;
       break;
     case 2:
       T = gsl_odeiv_step_rk2;
       break;
     case 3:
       T = gsl_odeiv_step_rk4;
       break;
     case 4:
       T = gsl_odeiv_step_rk8pd;
       break;
     default:
       printf("INVALID ode_method in getW() for Test P: must be 1,2,3,4\n");
       return TESTP_INVALID_METHOD;
   }

   s = gsl_odeiv_step_alloc(T, (size_t)1);     /* one scalar ode */
   c = gsl_odeiv_control_y_new(EPS_ABS,EPS_REL);
   e = gsl_odeiv_evolve_alloc((size_t)1);    /* one scalar ode */

   /* initial conditions: (r,W) = (R0,W_c(L^-));  r decreases from L toward 0 */
   rr = L;
   for (count = 0; count < N; count++) {
     /* generally use value at end of last interval as initial guess */
     W[count] = (count == 0) ? initialconditionW() : W[count-1];
     while (rr > r[count]) {
       step = r[count] - rr;
       status = gsl_odeiv_evolve_apply(e, c, s, &sys, &rr, r[count], &step, &W[count]);
       if (status != GSL_SUCCESS)   break;
     }
   }

   gsl_odeiv_evolve_free(e);
   gsl_odeiv_control_free(c);
   gsl_odeiv_step_free(s);
   return status;
}


int exactP(double r, double *h, double *magvb, double *W,
           const double EPS_ABS, const double EPS_REL, const int ode_method) {

  if (r > L)   return TESTP_R_OUT_OF_RANGE;
  if (r < 0.0) return TESTP_R_OUT_OF_RANGE;

  *h = h0 * (1.0 - (r/R0) * (r/R0));
  if (r > R1)
    *magvb = v0 * pow((r - R1)/(R0 - R1),5.0);

  return getW(&r,1,W,EPS_ABS,EPS_REL,ode_method);
}

#if 0
int exactP_list(double *r, int N, double *h, double *magvb, double *W, 
           const double EPS_ABS, const double EPS_REL, const int ode_method);
  /* N values r[0] > r[1] > ... > r[N-1]  (decreasing)
     assumes r, h, magvb, W are allocated length N arrays  */

  double *W;
  int stat, i;

  W = (double *) malloc((size_t)N * sizeof(double)); /* temporary arrays */

  /* combination EPS_ABS = 1e-12, EPS_REL=0.0, method = 1 = RK Cash-Karp
     believed to be predictable and accurate */
FIXME  stat = getW(r,N,FIXME,1.0e-12,0.0,1);
  if (stat != GSL_SUCCESS) {
    return stat;
  }

  for (i = 0; i < N; i++) {
    h[i] = FIXME
    magvb[i] = FIXME
  }

  free(W);
  return 0;
}
#endif

