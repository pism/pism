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

#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a > b ? b : a)

#define pi       3.1415926535897931
#define SperA    31556926.0    /* seconds per year; 365.2422 days */
#define g        9.81
#define rho      910.0         /* ice density; kg/m^3 */
#define rhow     1028.0        /* sea water density; kg/m^3 */
#define n        3.0           /* Glen power */
#define barB     3.7e8         /* strength of shelf; Pa s^(1/3); as in Schoof 2006;
                                  compare 1.9e8 from MacAyeal et al 1996 */
#define H0       500.0         /* m */
#define Rg       300.0e3       /* m;    300 km */
#define Rc       600.0e3       /* m;    600 km */


double F_M(double x, double alpha, double r, double Q) {
  const double 
     aor = alpha / r,
     DD = x * x + x * aor + pow(aor, 2.0); 
  return Q * pow(DD, 1./3.) - 2.0 * r * x - alpha;
}


double dF_M(double x, double alpha, double r, double Q) {
  const double 
     aor = alpha / r,
     DD = x * x + x * aor + pow(aor, 2.0); 
  return (1. / 3.) * Q * pow(DD, - 2./3.) * (2.0 * x + aor) - 2.0 * r;
}


int funcM_ode_G(double r, const double alpha[], double f[], void *params) {
  /*   RHS G for differential equation:
          alpha' = G(alpha,r)      
     but where we solve this equation to find alpha':
          F(alpha',alpha,r) = 0 
     heuristic: guess is about 1/7 th of solution to a nearby problem;
     no range checking on r, so use away from zero */
  const double Q = (1.0 - rho / rhow) * rho * g * Rc * H0 / (2.0 * barB),
               guess = 0.15 * (  pow(Q/r,n) - alpha[0]/r  );
  /* in Python (exactM.py):  f[0] = fsolve(F_M,guess,args=(alpha[0],r));
     we could call GSL to find root, but hand-coding Newton's is easier */
  double old = guess, new;
  int i;
  for (i = 1; i < 100; i++) {
    new = old - F_M(old,alpha[0],r,Q) / dF_M(old,alpha[0],r,Q);
    if (fabs((new-old)/old) < 1.0e-12)   break;
    old = new;
  }
  if (i >= 90)
    printf("exactTestM WARNING: Newton iteration not converged in funcM_ode_G!\n");
  f[0] = new;
  return GSL_SUCCESS;
}


#define NOT_DONE       8966
#define INVALID_METHOD 8968
#define NEGATIVE_R     8969

/* combination EPS_ABS = 1e-12, EPS_REL=0.0, method = 1 = RK Cash-Karp
 is believed to be predictable and accurate; returns GSL_SUCCESS=0 if success */
int exactM(double r,
           double *alpha, double *Drr,
           const double EPS_ABS, const double EPS_REL, const int ode_method) {

   double ug = 100.0 / SperA;  /* velocity across grounding line is 100 m/a */

   if (r < 0) {
     return NEGATIVE_R;  /* only nonnegative radial coord allowed */
   } else if (r <= Rg/4.0) {
     *alpha = 0.0;  /* zero velocity near center */
     *Drr = 0.0;
     return GSL_SUCCESS;
   } else if (r <= Rg) {
     /* quartic transition from alpha=0 to alpha=ug in   Rg/4 < r <= Rg;
        f(r) w: f(Rg/4)=f'(Rg/4)=0 and f(Rg)=ug and f(Rg) = DrrRg         */
     double DrrRg;
     funcM_ode_G(Rg, &ug, &DrrRg, NULL);  /* first get Drr = alpha' at Rg where alpha=ug */
     /* printf("DrrRg=%e (1/a)\n",DrrRg*SperA); */
     const double xx = r - 0.25 * Rg,
                  AA = 0.75 * Rg,
                  aa = (DrrRg - 2.0 * ug / AA) / (2.0 * AA * sqrt(ug)),
                  yy = aa * (xx - AA) + sqrt(ug) / AA;
     *alpha = xx * xx * yy * yy;
     *Drr = 2.0 * xx * yy * (yy + xx * aa);
     return GSL_SUCCESS;
   } else if (r >= Rc + 1.0) {
     *alpha = 0.0;  /* zero velocity beyond calving front */
     *Drr = 0.0;
     return GSL_SUCCESS;
   }
   
   /* need to solve ODE to find alpha, so setup for GSL ODE solver  */
   const gsl_odeiv_step_type* T;
   switch (ode_method) {
     case 1:
       T = gsl_odeiv_step_rkck; /* RK Cash-Karp */
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
   gsl_odeiv_system sys = {funcM_ode_G, NULL, 1, NULL};  /* Jac-free method and no params */

   /* initial conditions: (r,alf) = (Rg,ug);  r increases */
   double rr = Rg; 
   double myalf = ug;
   /* printf (" r (km)        alpha (m/a)\n");
      printf (" %11.5e   %11.5e\n", rr/1000.0, myalf * SperA); */
   double step;
   int status = NOT_DONE;
   while (rr < r) {
     /* step = r - rr;  try to get to solution in one step; trust stepping algorithm */
     step = MIN(r-rr,20.0e3);
     status = gsl_odeiv_evolve_apply(e, c, s, &sys, &rr, r, &step, &myalf);
     if (status != GSL_SUCCESS)   break;
     /* printf (" %11.5e   %11.5e\n", rr/1000.0, myalf * SperA); */
   }

   gsl_odeiv_evolve_free(e);
   gsl_odeiv_control_free(c);
   gsl_odeiv_step_free(s);

   *alpha = myalf;
   funcM_ode_G(r, alpha, Drr, NULL);
   return status;
}

