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

#define H0       500.0         /* m */
#define Rg       300.0e3       /* m;    300 km */
#define Rc       600.0e3       /* m;    600 km */


int funcM(double r, const double alpha[], double f[], void *params) {
  /*  RHS for differential equation:
      dalpha   
      ------ =       
        dr     
  */
  if ((r > Rg) && (r < Rc)) {
    f[0] = 0.0;  /* FIXME!! */
  } else {
    f[0] = 0.0;  /* no changes outside of defined interval */
  }
  return GSL_SUCCESS;
}



#define NOT_DONE       8966
#define INVALID_METHOD 8968
#define NEGATIVE_R     8969

/* combination EPS_ABS = 1e-12, EPS_REL=0.0, method = 1 = RK Cash-Karp
   is believed to be predictable and accurate ??? */
/* returns GSL_SUCCESS=0 if success */
int exactM(double r,
           double *alpha,
           const double EPS_ABS, const double EPS_REL, const int ode_method) {

   double ug = 100.0 / SperA;  /* velocity across grounding line is 100 m/a */

   if (r < 0) {
     return NEGATIVE_R;  /* only nonnegative radial coord allowed */
   } else if (r <= Rg/4.0) {
     *alpha = 0.0;  /* zero velocity near center */
     return GSL_SUCCESS;
   } else if (r <= Rg) {
     /* smooth transition from alpha=0 to alpha=ug in   Rg/4 < r <= Rg  */
     double ratio = (r - 0.25 * Rg) / 0.75 * Rg;
     *alpha = (ug / 2.0) * (1.0 - cos(pi * ratio));   
     return GSL_SUCCESS;
   } else if (r >= Rc) {
     *alpha = 0.0;  /* zero velocity beyond calving front */
     return GSL_SUCCESS;
   }
   
   /* need to solve ODE to find alpha, so setup for GSL ODE solver  */
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
       printf("INVALID ode_method in exactM(): must be 1,2,3,4\n");
       return INVALID_METHOD;
   }
   gsl_odeiv_step* s = gsl_odeiv_step_alloc(T, 1);     /* one scalar ode */
   gsl_odeiv_control* c = gsl_odeiv_control_y_new(EPS_ABS,EPS_REL);
   gsl_odeiv_evolve* e = gsl_odeiv_evolve_alloc(1);    /* one scalar ode */
   gsl_odeiv_system sys = {funcM, NULL, 1, NULL};  /* Jac-free method and no params */

   /* initial conditions: (r,alf) = (Rg,ug);  r increases */
   double rr = Rg; 
   double myalf = ug;
   printf (" r (km)        alpha (m/a)\n");
   printf ("%12.5e %12.5e\n", rr/1000.0, myalf * SperA);
   double step;
   int status = NOT_DONE;
   while (rr < r) {
     step = r - rr;  /* try to get to solution in one step; trust stepping algorithm */
     status = gsl_odeiv_evolve_apply(e, c, s, &sys, &rr, r, &step, &myalf);
     if (status != GSL_SUCCESS)   break;
     printf ("%12.5e %12.5e\n", rr/1000.0, myalf * SperA);
   }

   gsl_odeiv_evolve_free(e);
   gsl_odeiv_control_free(c);
   gsl_odeiv_step_free(s);

   *alpha = myalf;
   return status;
}

