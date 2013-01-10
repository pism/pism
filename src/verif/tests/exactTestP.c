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
#include <gsl/gsl_odeiv2.h>
#include "exactTestP.h"

#define SperA    31556926.0    /* seconds per year; 365.2422 days */
#define g        9.81          /* m s-2 */
#define rhoi     910.0         /* kg m-3 */
#define rhow     1000.0        /* kg m-3 */

/* major model parameters: */
#define Aglen    3.1689e-24    /* Pa-3 s-1 */
#define K        0.01          /* m s-1 */
#define Wr       1.0           /* m */
#define c1       0.500         /* m-1 */
#define c2       0.040         /* [pure] */

/* model regularizations */
#define E0       1.0           /* m */
#define Y0       0.001         /* m */

/* specific to exact solution */
#define Phi0     (0.20 / SperA) /* m s-1 */
#define h0       500.0         /* m */
#define v0       (100.0 / SperA) /* m s-1 */
#define R1       5000.0        /* m */


int getsb(double r, double *sb, double *dsbdr) {
  double CC, CZ, CD, zz;
  if (r < R1) {
    *sb    = 0.0;
    *dsbdr = 0.0;
  } else {
    CC     = pow( (c1 * v0) / (c2 * Aglen * pow((TESTP_L - R1),5.0)) , (1.0/3.0) );
    *sb    = CC * pow(r - R1, (5.0/3.0));
    *dsbdr = (5.0/3.0) * CC * pow(r - R1, (2.0/3.0));
  }
  return 0;
}


double criticalW(double r) {
  double h = h0 * (1.0 - (r/TESTP_R0) * (r/TESTP_R0)),
         Po = rhoi * g * h,
         sb, dsb, sbcube, Pocube;
  getsb(r,&sb,&dsb);
  sbcube = sb * sb * sb;
  Pocube = Po * Po * Po;
  return ((sbcube * Wr - Pocube * Y0) / (sbcube + Pocube));
}


int funcP(double r, const double W[], double f[], void *params) {
  /* Computes RHS f(r,W) for differential equation as given in dampnotes.pdf
  at https://github.com/bueler/hydrolakes:
      dW
      -- = f(r,W)
      dr
  Compare doublediff.m.  Assumes Glen power n=3.
  */

  double sb, dsb, dPo, tmp1, c0, vphi0, numer, denom;

  if (params == NULL) {} /* quash warning "unused parameters" */

  if (r < 0.0) {
    f[0] = 0.0;  /* place-holder */
    return TESTP_R_NEGATIVE;
  } else if (r > TESTP_L) {
    f[0] = 0.0;
    return GSL_SUCCESS;
  } else {
    getsb(r,&sb,&dsb);
    c0    = K / (rhow * g);
    vphi0 = Phi0 / (2 * c0);
    dPo   = - (2.0 * rhoi * g * h0 / (TESTP_R0*TESTP_R0)) * r;
    tmp1  = pow(W[0] + Y0,4.0/3.0) * pow(Wr - W[0],2.0/3.0);
    numer = dsb * (W[0] + Y0) * (Wr - W[0]);
    numer = numer - ( vphi0 * r / W[0] + dPo ) * tmp1;
    denom = (1.0/3.0) * (Wr + Y0) * sb + rhow * g * tmp1;
    f[0] = numer / denom;
    return GSL_SUCCESS;
  }
}


/* Computes initial condition W(r=L) = W_c(L^-). */
double initialconditionW() {
  double hL, PoL, sbL;
  hL  = h0 * (1.0 - (TESTP_L/TESTP_R0) * (TESTP_L/TESTP_R0));
  PoL = rhoi * g * hL;
  sbL = pow( c1 * v0 / (c2 * Aglen), 1.0/3.0);
  return (pow(sbL,3.0) * Wr - pow(PoL,3.0) * Y0) / (pow(sbL,3.0) + pow(PoL,3.0));
}


double psteady(double W, double magvb, double Po) {
  double sbcube, frac, P;
  sbcube = c1 * fabs(magvb) / (c2 * Aglen);
  frac = (W < Wr) ? (Wr - W) / (W + Y0) : 0.0;
  P = Po - pow(sbcube * frac, 1.0/3.0);
  if (P < 0.0)  P = 0.0;
  return P;
}


/* Solves ODE for W(r), the exact solution.  Input r[] and output W[] must be
allocated vectors of length N.  Input r[] must be decreasing.  The combination
EPS_ABS = 1e-12, EPS_REL=0.0, method = RK Dormand-Prince O(8)/O(9)
is believed for now to be predictable and accurate.  Note hstart is negative
so that the ODE solver does negative steps.  Assumes
   0 <= r[N-1] <= r[N-2] <= ... <= r[1] <= r[0] <= L.                            */
int getW(double *r, int N, double *W,
         const double EPS_ABS, const double EPS_REL, const int ode_method) {
   int i, count;
   int status = TESTP_NOT_DONE;
   double rr, hstart;
   const gsl_odeiv2_step_type* Tpossible[4];
   const gsl_odeiv2_step_type *T;
   gsl_odeiv2_system sys = {funcP, NULL, 1, NULL};  /* Jac-free method and no params */
   gsl_odeiv2_driver *d;

   /* setup for GSL ODE solver; these choices don't need Jacobian */
   Tpossible[0] = gsl_odeiv2_step_rk8pd;
   Tpossible[1] = gsl_odeiv2_step_rk2;
   Tpossible[2] = gsl_odeiv2_step_rkf45;
   Tpossible[3] = gsl_odeiv2_step_rkck;
   if ((ode_method > 0) && (ode_method < 5))
     T = Tpossible[ode_method-1];
   else {
     printf("INVALID ode_method in getW(): must be 1,2,3,4\n");
     return TESTP_INVALID_METHOD;
   }

   hstart = -1000.0;
   d = gsl_odeiv2_driver_alloc_y_new(&sys, T, hstart, EPS_ABS, EPS_REL);

   /* initial conditions: (r,W) = (L,W_c(L^-));  r decreases from L toward 0 */
   rr = TESTP_L;
   for (count = 0; count < N; count++) {
     /* except at start, use value at end of last interval as initial value for subinterval */
     W[count] = (count == 0) ? initialconditionW() : W[count-1];
     while (rr > r[count]) {
       status = gsl_odeiv2_driver_apply(d, &rr, r[count], &(W[count]));
       if (status != GSL_SUCCESS) {
         printf("gsl_odeiv2_driver_apply() returned status = %d\n",status);
         break;
       }
       if (W[count] > Wr) {
         return TESTP_W_EXCEEDS_WR;
       } else if (W[count] < criticalW(r[count])) {
         return TESTP_W_BELOW_WCRIT;
       }
     }
   }

   gsl_odeiv2_driver_free(d);
   return status;
}


int exactP(double r, double *h, double *magvb, double *Wcrit, double *W, double *P,
           const double EPS_ABS, const double EPS_REL, const int ode_method) {

  int status;
  if (r < 0.0) return TESTP_R_NEGATIVE;
  if (r > TESTP_L) {
    *h     = 0.0;
    *magvb = 0.0;
    *Wcrit = 0.0;
    *W     = 0.0;
    *P     = 0.0;
    return 0;
  }

  *h = h0 * (1.0 - (r/TESTP_R0) * (r/TESTP_R0));
  if (r > R1)
    *magvb = v0 * pow((r - R1)/(TESTP_L - R1),5.0);
  else
    *magvb = 0.0;
  *Wcrit = criticalW(r);

  if (r == TESTP_L) {
    *W = initialconditionW();
    *P = psteady(*W, *magvb, rhoi * g * (*h));
    return 0;
  } else {
    status = getW(&r,1,W,EPS_ABS,EPS_REL,ode_method);
    if (status) {
      *P = 0;
      return status;
    } else {
      *P = psteady(*W, *magvb, rhoi * g * (*h));
      return 0;
    }
  }
}



int exactP_list(double *r, int N, double *h, double *magvb, double *Wcrit, double *W, double *P,
                const double EPS_ABS, const double EPS_REL, const int ode_method) {

  int i, M, status;
  /* check first: we have a list, r is decreasing, r is in range [0,L) */
  if (N < 1) return TESTP_NO_LIST;
  for (i = 0; i<N; i++) {
    if ((i > 0) && (r[i] > r[i-1])) {
      printf("error:  r[%d] = %.18e > r[%d] = %.18e\n",i,r[i],i-1,r[i-1]);
      printf("   ... so list not decreasing ... returning error code ...\n");
      return TESTP_LIST_NOT_DECREASING;
    }
    if (r[i] < 0.0)  return TESTP_R_NEGATIVE;
  }

  M = 0;
  while (r[M] > TESTP_L) {
      h[M]     = 0.0;
      magvb[M] = 0.0;
      Wcrit[M] = 0.0;
      W[M]     = 0.0;
      P[M]     = 0.0;
      M++;
  }

  for (i = M; i<N; i++) {
    h[i] = h0 * (1.0 - (r[i]/TESTP_R0) * (r[i]/TESTP_R0));
    if (r[i] > R1)
      magvb[i] = v0 * pow((r[i] - R1)/(TESTP_L - R1),5.0);
    else
      magvb[i] = 0.0;
    Wcrit[i] = criticalW(r[i]);
  }

  status = getW(&(r[M]),N-M,&(W[M]),EPS_ABS,EPS_REL,ode_method);

  if (status) {
    for (i = M; i<N; i++)
      P[i] = 0.0;
    return status;
  } else {
    for (i = M; i<N; i++)
      P[i] = psteady(W[i], magvb[i], rhoi * g * h[i]);
    return 0;
  }

}


int error_message_testP(int status) {
  switch (status) {
    case TESTP_R_NEGATIVE:
      printf("error in Test P: r < 0\n");
      break;
    case TESTP_W_EXCEEDS_WR:
      printf("error in Test P: W > W_r\n");
      break;
    case TESTP_W_BELOW_WCRIT:
      printf("error in Test P: W < W_crit\n");
      break;
    case TESTP_INVALID_METHOD:
      printf("error in Test P: invalid choice for ODE method\n");
      break;
    case TESTP_NOT_DONE:
      printf("error in Test P: ODE integrator not done\n");
      break;
    case TESTP_NO_LIST:
      printf("error in Test P: no list of r values at input to exactP_list()\n");
      break;
    case TESTP_LIST_NOT_DECREASING:
      printf("error in Test P: input list of r values to exactP_list() is not decreasing\n");
      break;
    default:
      if (status > 0) printf("unknown error status in Test P\n");
  }
  return 0;
}

