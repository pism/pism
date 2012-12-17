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
#include "exactTestL.h"

#define pi       3.1415926535897931
#define SperA    31556926.0    /* seconds per year; 365.2422 days */

#define Phi0     0.020         /* m a-1 */
#define h0         500.0       /* m */
#define R0       25000.0       /* m */
#define R1        5000.0       /* m */
#define v0         100.0       /* m a-1 */

#define g        9.81          /* m s-2 */
#define rhoi     910.0         /* kg m-3 */
#define nglen    3.0           /* [pure] */

int funcP(double r, const double W[], double f[], void *params) {
  /*
  RHS f(r,W) for differential equation:
      dW
      -- = f(r,W) = FIXME
      dr
  */

  /*const double A = 1.0e-16 / SperA;  = 3.17e-24  1/(Pa^3 s); EISMINT I flow law */

  if (params == NULL) {} /* quash warning "unused parameters" */

  if ((r >= 0.0) && (r <= R0)) {
    /*const double freq = z0 * pi / L;*/
    f[0] =  FIXME;
  } else {
    f[0] = 0.0;  /* no changes outside of defined interval */
  }
  return GSL_SUCCESS;
}


/* combination EPS_ABS = 1e-12, EPS_REL=0.0, method = 1 = RK Cash-Karp
   is believed to be predictable and accurate */
int getW(double *r, int N, FIXME double *u,
         const double EPS_ABS, const double EPS_REL, const int ode_method) {
   /* solves ODE for FIXME
      r and u must be allocated vectors of length N; r[] must be decreasing */
FIXME
   int i, k, count;
   const gsl_odeiv_step_type* T;
   int status = TESTL_NOT_DONE;
   double rr, step;

   gsl_odeiv_step*    s;
   gsl_odeiv_control* c;
   gsl_odeiv_evolve*  e;
   gsl_odeiv_system   sys = {funcL, NULL, 1, NULL};  /* Jac-free method and no params */

   /* check first: we have a list, and r is decreasing */
   if (N < 1) return TESTL_NO_LIST;
   for (i = 1; i<N; i++) {  if (r[i] > r[i-1]) return TESTL_NOT_DECREASING;  }

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
       printf("INVALID ode_method in getU(): must be 1,2,3,4\n");
       return TESTL_INVALID_METHOD;
   }

   s = gsl_odeiv_step_alloc(T, (size_t)1);     /* one scalar ode */
   c = gsl_odeiv_control_y_new(EPS_ABS,EPS_REL);
   e = gsl_odeiv_evolve_alloc((size_t)1);    /* one scalar ode */


   /* outside of ice cap, u = 0 */
   k = 0;
   while (r[k] >= L) {
     u[k] = 0.0;
     k++;
     if (k == N) return GSL_SUCCESS;
   }

   /* initial conditions: (r,u) = (L,0);  r decreases from L */
   rr = L;
   for (count = k; count < N; count++) {
     /* generally use value at end of last interval as initial guess */
     u[count] = (count == 0) ? 0.0 : u[count-1];
     while (rr > r[count]) {
       step = r[count] - rr;
       status = gsl_odeiv_evolve_apply(e, c, s, &sys, &rr, r[count], &step, &u[count]);
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

  double h[1] = { 0.0 };
  double magvb[1] = { 0.0 };
  double W[1] = { 0.0 };
  double s = r / R0;

  if (r >= R0) return 0;

  *h = h0 * (1.0 - s * s);
  if (r > R1)
    *magvb = v0 * pow((r - R1)/(R0 - R1),5.0);

FIXME  getW(&r,1,FIXME,EPS_ABS,EPS_REL,ode_method);
  return 0;
}


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

