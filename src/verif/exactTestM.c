/*
   Copyright (C) 2008 Ed Bueler
  
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
#include "exactTestM.h"

#define pi       3.1415926535897931
#define SperA    31556926.0    /* seconds per year; 365.2422 days */
#define g        9.81
#define rho      910.0
#define n        3.0           /* Glen power */

#define Rg       300.0e3       /* m;    300 km */
#define Rc       700.0e3       /* m;    700 km */
#define H0       500.0         /* m */


FIXME ALL FROM HERE DOWN; SEE exactM.py  !!

int func(double r, const double u[], double f[], void *params) {
  /*
  RHS for differential equation:
      du                  5/8   / a_0  r  (L^2 - r^2) \ 1/3
      -- = - (8/3) b'(r) u    - |---------------------|       
      dr                        \ 2 L^2 \tilde\Gamma  /
  */
  
  const double Lsqr = L * L;
  const double a0 = 0.3 / SperA;   /* m/s;  i.e. 0.3 m/a */
  const double A = 1.0e-16 / SperA;  /* = 3.17e-24  1/(Pa^3 s); EISMINT I flow law parameter */
  const double Gamma = 2 * pow(rho * g,n) * A / (n+2);
  const double tilGamma = Gamma * pow(n,n) / (pow(2.0 * n + 2.0, n));
  const double C = a0 / (2.0 * Lsqr * tilGamma);
  if ((r >= 0.0) && (r <= L)) {
    const double freq = z0 * pi / L;
    const double bprime = b0 * freq * sin(freq * r);
    f[0] =  - (8.0/3.0) * bprime * pow(u[0], 5.0/8.0)
            - pow(C * r * (Lsqr - r * r), 1.0/3.0);
  } else {
    f[0] = 0.0;  /* no changes outside of defined interval */
  }
  return GSL_SUCCESS;
}


#define NOT_DONE       8966
#define NOT_DECREASING 8967
#define INVALID_METHOD 8968

/* combination EPS_ABS = 1e-12, EPS_REL=0.0, method = 1 = RK Cash-Karp
   is believed to be predictable and accurate */
int getU(double *r, int N, double *u, 
         const double EPS_ABS, const double EPS_REL, const int ode_method) {
   /* solves ODE for u(r)=H(r)^{8/3}, 0 <= r <= L, for test L
      r and u must be allocated vectors of length N; r[] must be decreasing */

   /* check r is decreasing first */
   int i;
   for (i = 1; i<N; i++) {
     if (r[i] > r[i-1]) {
       printf("r[] not decreasing in getU()\n");
       return NOT_DECREASING;
     }
   }

   /* setup for GSL ODE solver; following step choices don't need Jacobian,
      but should we chose one that does?  */
	const gsl_odeiv_step_type* T;
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
       return INVALID_METHOD;
   }
	gsl_odeiv_step* s = gsl_odeiv_step_alloc(T, 1);     /* one scalar ode */
	gsl_odeiv_control* c = gsl_odeiv_control_y_new(EPS_ABS,EPS_REL);
	gsl_odeiv_evolve* e = gsl_odeiv_evolve_alloc(1);    /* one scalar ode */
	gsl_odeiv_system sys = {func, NULL, 1, NULL};  /* Jac-free method and no params */

   /* initial conditions: (r,u) = (L,0);  r decreases from L */
	double rr = L, step;
   int status = NOT_DONE;
   int count  = 0;
   for (count = 0; count < N; count++) {
     if (count == 0) {
       u[count] = 0;
     } else {
       u[count] = u[count-1];  /* use value at end of last interval as initial */
     }
     while (rr > r[count]) {
       step = r[count] - rr;
       status = gsl_odeiv_evolve_apply(e, c, s, &sys, &rr, r[count], &step, &u[count]);
	    if (status != GSL_SUCCESS)   break;
	  }
     /* printf ("%d  %d  %d  %.5e %.5e\n", status, GSL_SUCCESS, count, r[count], u[count]); */
	}

	gsl_odeiv_evolve_free(e);
	gsl_odeiv_control_free(c);
	gsl_odeiv_step_free(s);
   return status;
}


/*int exactL(double r, double *H, double *b, double *a, 
           const double EPS_ABS, const double EPS_REL, const int ode_method) {*/
int exactM(double r,
           double *alpha,
           const double EPS_ABS, const double EPS_REL, const int ode_method) {

  double u[1] = { 0.0 };
  getU(&r,1,u,EPS_ABS,EPS_REL,ode_method);
  *alpha = 
  return 0;
}


