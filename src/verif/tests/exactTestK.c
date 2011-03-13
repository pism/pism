/*
   Copyright (C) 2007, 2011 Ed Bueler
  
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
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "exactTestK.h"

#define pi             3.1415926535897931
#define SperA          31556926.0   /* seconds per year; 365.2422 days */

#define c_p_ICE        2009.0       /* J/(kg K)  specific heat capacity of ice */
#define rho_ICE        910.0        /* kg/(m^3)  density of ice */
#define k_ICE          2.10         /* J/(m K s) = W/(m K)  thermal conductivity of ice */
#define c_p_BRdefault  1000.0       /* J/(kg K)  specific heat capacity of bedrock */
#define rho_BRdefault  3300.0       /* kg/(m^3)  density of bedrock */
#define k_BRdefault    3.0          /* J/(m K s) = W/(m K)  thermal conductivity of bedrock */

#define H0             3000.0       /* m */
#define B0             1000.0       /* m */
#define Ts             223.15       /* K */
#define G              0.042        /* W/(m^2) */
#define phi            0.0125       /* K/m */

#define Nsum           30           /* number of terms in eigenfunction expansion; the exact
                                       solution is deliberately chosen to have finite expansion */


int exactK(const double t, const double z, double *TT, double *FF, const int bedrockIsIce_p) {
  int k;
  int belowB0;
  double ZZ, alpha, lambda, beta, my_gamma, XkSQR, Xk,
         theta, dthetakdz, P, dPdz,
         Ck, I1, I2, aH, bB, mI, mR;
  double c_p_BR, rho_BR, k_BR;
  /* following constants were produced by calling print_alpha_k(30) (below) */
  double alf[Nsum] = {3.350087528822397e-04, 1.114576827617396e-03, 1.953590840303518e-03,
                      2.684088585781064e-03, 3.371114869333445e-03, 4.189442265117592e-03,
                      5.008367405382524e-03, 5.696044031764593e-03, 6.425563506942886e-03,
                      7.264372872913219e-03, 8.044853066396166e-03, 8.714877612414516e-03,
                      9.493529164160654e-03, 1.033273985210279e-02, 1.106421822502108e-02,
                      1.175060460132703e-02, 1.256832682090360e-02, 1.338784224692084e-02,
                      1.407617951778051e-02, 1.480472324161026e-02, 1.564331999062109e-02,
                      1.642470780103220e-02, 1.709475346624607e-02, 1.787248418996684e-02,
                      1.871188358061674e-02, 1.944434477688470e-02, 2.013010181370026e-02,
                      2.094721145334310e-02, 2.176730968036079e-02, 2.245631776169424e-02};

  if (bedrockIsIce_p) {
    c_p_BR = c_p_ICE;
    rho_BR = rho_ICE;
    k_BR = k_ICE;
    for (k = 0; k < Nsum; k++) { /* overwrite alpha_k with ice-meets-ice values; see preprint */
      alf[k] = (2.0 * k + 1.0) * pi / (2.0 * (H0 + B0));
    }
  } else {
    c_p_BR = c_p_BRdefault;
    rho_BR = rho_BRdefault;
    k_BR = k_BRdefault;
  }
  if (z > H0) {
    *TT = Ts;
    return 0;
  }
  belowB0 = (z < -B0);

  ZZ = sqrt((rho_BR * c_p_BR * k_ICE) / (rho_ICE * c_p_ICE * k_BR));
  mI = (G / k_ICE) - phi;     mR = (G / k_BR) - phi;
  /* DEBUG: printf("ZZ = %10e, mI = %10e, mR = %10e\n", ZZ,mI,mR); */
  *TT = 0.0;
  *FF = 0.0;
  for (k = Nsum-1; k >= 0; k--) {
    /* constants only having to do with eigenfunctions; theta = theta_k(z) is the
       normalized eigenfunction */ 
    alpha = alf[k];
    beta = ZZ * alpha;
    my_gamma = sin(alpha * H0) / cos(beta * B0);
    XkSQR = (rho_BR * c_p_BR * my_gamma * my_gamma * B0 + rho_ICE * c_p_ICE * H0) / 2.0;
    Xk = sqrt(XkSQR);
    /* theta = ( (z > 0) ? sin(alpha * (H0 - z)) : my_gamma * cos(beta * (B0 + z)) ) / Xk; */
    theta = (z > 0) ? sin(alpha * (H0 - z))
                    : my_gamma * cos(beta * (B0 + z)); 
    theta /= Xk;
    dthetakdz = (z > 0) ? - alpha * cos(alpha * (H0 - z)) 
                        : - beta * my_gamma * sin(beta * (B0 + z));
    dthetakdz /= Xk;
    lambda = (k_ICE * alpha * alpha) / (rho_ICE * c_p_ICE);
    /* DEBUG: printf("k = %3d:  alpha = %10e, Xk = %10e, theta = %10e, dthetakdz = %10e, lambda = %10e,\n",
           k,alpha,Xk,theta,dthetakdz,lambda); */
    /* constants involved in computing the expansion coefficients */
    aH = alpha * H0;            bB = beta * B0;
    I1 = - mI * (sin(aH) - aH * cos(aH)) / (alpha * alpha);
    I2 = mR * (cos(bB) - 1.0 + bB * sin(bB)) / (beta * beta)
         - (B0 * mR + H0 * mI) * sin(bB) / beta;
    Ck = (rho_ICE * c_p_ICE * I1 + rho_BR * c_p_BR * my_gamma * I2) / Xk;
    /* add the term to the expansion */
    *TT += Ck * exp(- lambda * t) * theta;
    *FF += - ((z > 0) ? k_ICE : k_BR) * Ck * exp(- lambda * t) * dthetakdz;
    /* DEBUG: printf("          I1 = %10e, I2 = %10e, Ck = %10e, term = %10f\n",
           I1,I2,Ck, Ck * exp(- lambda * t) * theta ); */
  }
  /* P = (z >= 0) ? (z / k_ICE) - (H0 / k_ICE) : (z / k_BR) - (H0 / k_ICE); */
  P = (z / ((z > 0) ? k_ICE : k_BR)) - (H0 / k_ICE);
  dPdz = 1.0 / ((z > 0) ? k_ICE : k_BR);
  *TT += Ts - G * P;
  *FF += ((z > 0) ? k_ICE : k_BR) * G * dPdz;

  return belowB0;

}


#define COMPUTE_ALPHA 0
#if COMPUTE_ALPHA

#define ALPHA_RELTOL   1.0e-14
#define ITER_MAXED_OUT 999

/* parameters needed for root problem: */
struct coscross_params {
  double Afrac, HZBsum, HZBdiff;
};

/* the root problem is to make this function zero: */
double coscross(double alpha, void *params) {
  struct coscross_params *p = (struct coscross_params *) params;
  return cos(p->HZBsum * alpha) - p->Afrac * cos(p->HZBdiff * alpha);
}

/* compute the first N roots alpha_k of the equation
     ((A-1)/(A+1)) cos((H - Z B) alpha) = cos((H + Z B) alpha)
where H and B are heights and A, Z are defined in terms of material 
constants */
int print_alpha_k(const int N) {
  int status, iter, k, max_iter = 200;
  double Z, A;
  double alpha, alpha_lo, alpha_hi, temp_lo;
  const gsl_root_fsolver_type *solvT;
  gsl_root_fsolver *solv;
  gsl_function F;
  struct coscross_params params;
  
  Z = sqrt((rho_BR * c_p_BR * k_ICE) / (rho_ICE * c_p_ICE * k_BR));
  A = (k_BR / k_ICE) * Z;
  params.Afrac   = (A - 1.0) / (A + 1.0);     
  params.HZBsum  = H0 + Z * B0;
  params.HZBdiff = H0 - Z * B0;
     
  F.function = &coscross;
  F.params = &params;
  solvT = gsl_root_fsolver_brent;  /* faster than bisection but still bracketing */
  solv = gsl_root_fsolver_alloc(solvT);

  for (k = 0; k < N; k++) {
    /* these numbers bracket exactly one solution */
    alpha_lo = (double(k) * pi) / params.HZBsum;
    alpha_hi = (double(k + 1) * pi) / params.HZBsum;
    gsl_root_fsolver_set(solv, &F, alpha_lo, alpha_hi);
     
    iter = 0;
    do {
      iter++;
      status = gsl_root_fsolver_iterate(solv);
      alpha = gsl_root_fsolver_root(solv);
      alpha_lo = gsl_root_fsolver_x_lower(solv);
      alpha_hi = gsl_root_fsolver_x_upper(solv);
      temp_lo = (alpha_lo > 0) ? alpha_lo : (alpha_hi/2.0);
      status = gsl_root_test_interval(temp_lo, alpha_hi, 0, ALPHA_RELTOL);
    } while ((status == GSL_CONTINUE) && (iter < max_iter));
    if (iter >= max_iter) {
      printf("!!!ERROR: root finding iteration reached maximum iterations; QUITING!\n");
      return ITER_MAXED_OUT;
    }
    printf("%19.15e,\n",alpha);
    /* DEBUG: printf("%19.15e  (in orig bracket [%19.15e,%19.15e])\n",alpha,
              (double(k) * pi) / params.HZBsum, (double(k+1) * pi) / params.HZBsum); */
  }
  
  gsl_root_fsolver_free(solv);
  return status;
}
#endif /* COMPUTE_ALPHA */

