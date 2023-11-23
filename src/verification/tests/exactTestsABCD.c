/*
   Copyright (C) 2004-2006, 2014, 2016, 2023 Jed Brown and Ed Bueler
  
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

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>       /* M_PI */
#include "pism/verification/tests/exactTestsABCD.h"

#define SperA 31556926.0  /* seconds per year; 365.2422 days */

int exactA_old(const double r, double *H, double *M) {
  /* NOTE: t is in seconds */
  const double L = 750000.0;       /* m; distance of margin from center */
  const double M0 = 0.3 / SperA;   /* 30 cm year-1 constant accumulation */
  const double g = 9.81;           /* m/s^2; accel of gravity */
  const double rho = 910.0;        /* kg/m^3; density */
  const double n = 3.0;            /* Glen exponent */
  const double A = 1.0E-16/SperA;  /* = 3.17e-24  1/(Pa^3 s); */
                                   /* (EISMINT value) flow law parameter */
  const double Gamma = 2 * pow(rho * g,n) * A / (n+2);
  
  double       m, p, C;
  
  if (r < L) {
    m = 2.0 * n + 2.0;
    p = 1.0 + 1.0 / n;
    C = pow(pow(2.0, n - 1) * M0 / Gamma, 1.0 / m);
    *H = C * pow(pow(L, p) - pow(r, p), n / m); 
  } else {
    *H = 0.0;
  }
  
  *M = M0;
  
  return 0;
}

struct TestABCDParameters exactA(double r) {
  struct TestABCDParameters result;

  result.error_code = exactA_old(r, &result.H, &result.M);

  return result;
}

int exactB_old(const double t, const double r, double *H, double *M) {
  /* NOTE: t and t0 are in seconds */
  double alpha, beta, t0, Rmargin;
  const double n = 3.0, H0 = 3600.0, R0=750000.0;
  
  /* lambda=0.0 case of Bueler et al (2005) family of similarity solns;
     is Halfar (1983) soln */
  alpha=1.0/9.0;  /* alpha=(2-(n+1)*lambda)/(5*n+3)=1/9 */
  beta=1.0/18.0;  /* beta=(1+(2*n+1)*lambda)/(5*n+3)=1/18 */
  t0=422.45*SperA;  /* t0 = (beta/Gamma)
                             * pow((2n+1)/(n+1),n)*(pow(R0,n+1)/pow(H0,2n+1)) */

  Rmargin=R0*pow(t/t0,beta);
  if (r<Rmargin)
    *H = H0 * pow(t/t0,-alpha) 
          * pow(1.0-pow(pow(t/t0,-beta)*(r/R0), (n+1)/n),  n/(2*n+1));
  else
    *H = 0.0;

  *M=0.0;
  return 0;
}

struct TestABCDParameters exactB(const double t, const double r) {
  struct TestABCDParameters result;

  result.error_code = exactB_old(t, r, &result.H, &result.M);

  return result;
}

int exactC_old(const double t, const double r, double *H, double *M) {
  double lambda, alpha, beta, t0, Rmargin;
  const double n = 3.0, H0 = 3600.0, R0=750000.0;

  lambda=5.0;
  alpha=-1.0;  /* alpha=(2-(n+1)*lambda)/(5*n+3) */
  beta=2.0;  /* beta=(1+(2*n+1)*lambda)/(5*n+3) */
  t0=15208.0*SperA;  /* t0 = (beta/Gamma)
                             * pow((2n+1)/(n+1),n)*(pow(R0,n+1)/pow(H0,2n+1)) */

  Rmargin=R0*pow(t/t0,beta);
  if (r<Rmargin)
    *H = H0 * pow(t/t0,-alpha)
          * pow(1.0-pow(pow(t/t0,-beta)*(r/R0), (n+1)/n),  n/(2*n+1));
  else
    *H = 0.0;

  if (t>0.1*SperA)
    *M = (lambda/t)* (*H);
  else {  /* when less than 0.1 year, avoid division by time */
    Rmargin=R0*pow(0.1*SperA/t0,beta);
    if (r<Rmargin)
      *M=5*H0/t0;  /* constant value in disc of Rmargin radius */
    else
      *M=0.0;
  }
  return 0;
}

struct TestABCDParameters exactC(const double t, const double r) {
  struct TestABCDParameters result;

  result.error_code = exactC_old(t, r, &result.H, &result.M);

  return result;
}

int exactD_old(const double t, const double rin, double *H, double *M) {

  /* parameters describing extent of sheet: */
  const double H0=3600.0;          /* m */
  const double L=750000.0;         /* m */
  /* parameters for perturbation: */
  const double Tp=5000.0*SperA;    /* s */
  const double Cp=200.0;           /* m */
  /* fundamental physical constants */
  const double g=9.81;             /* m/s^2; accel of gravity */
  /* ice properties; parameters which appear in constitutive relation: */
  const double rho=910.0;          /* kg/m^3; density */
  const double n=3.0;              /* Glen exponent */
  const double A=1.0E-16/SperA;    /* = 3.17e-24  1/(Pa^3 s); */
                                   /* (EISMINT value) flow law parameter */

  double       r=rin;
  double       Gamma, power, s, lamhat, f, Hs, temp, C, Ms, df, ddf, chi, dchi,
               ddchi, c1, dHs, ddHs, dH, ddH, divterms, Mc;
 
  if (r < 0.0) r=-r;
 
  if (r >= L - 0.01) {
    *H = 0.0;
    *M = -0.1 / SperA;
  } else {
          if (r < 0.01) r = 0.01;  /* avoid r=0 singularity in formulas */
          
	  /* important derived quantities */
	  Gamma = 2 * pow(rho * g,n) * A / (n+2);
	  power = n / (2*n+2);
	  s = r/L;

	  /* compute H from analytical steady state Hs plus perturbation */
	  lamhat = (1+1/n)*s - (1/n) + pow(1-s,1+1/n) - pow(s,1+1/n);
          if ((r>0.3*L) && (r<0.9*L))    f = pow(cos(M_PI*(r-0.6*L)/(0.6*L)) ,2.0);
	  else                           f = 0.0;
	  Hs = (H0 / pow(1-1/n,power)) * pow(lamhat,power);
          *H = Hs + Cp * sin(2.0*M_PI*t/Tp) * f;

	  /* compute steady part of accumulation */
	  temp = pow(s,1/n) + pow(1-s,1/n) - 1;
	  C = Gamma * pow(H0,2*n+2) / pow(2.0 * L * (1-1/n),n);
	  Ms = (C/r) * pow(temp,n-1)
	                *  (2*pow(s,1/n) + pow(1-s,(1/n)-1) * (1 - 2*s) - 1);

	  /* derivs of H */
	  if ((r>0.3*L) && (r<0.9*L)) {
            df = -(M_PI/(0.6*L)) * sin(M_PI*(r-0.6*L)/(0.3*L));
            ddf = -(M_PI*M_PI/(0.18*L*L)) * cos(M_PI*(r-0.6*L)/(0.3*L));
	    chi = (4.0/3.0)*s - 1.0/3.0 + pow(1-s,4.0/3.0) - pow(s,4.0/3.0);
	    dchi = -(4.0/(3.0*L)) * (pow(s,1.0/3.0) + pow(1-s,1.0/3.0) - 1);
	    ddchi = -(4.0/(9.0*L*L)) * (pow(s,-2.0/3.0) - pow(1-s,-2.0/3.0));
	    c1 = (3.0*H0) / (8.0*pow(2.0/3.0,3.0/8.0));
	    dHs = c1 * pow(chi,-5.0/8.0) * dchi;
	    ddHs = c1 * ((-5.0/8.0) * pow(chi,-13.0/8.0) * dchi * dchi
	                   + pow(chi,-5.0/8.0) * ddchi);
            dH = dHs + Cp * sin(2.0*M_PI*t/Tp) * df;
            ddH = ddHs + Cp * sin(2.0*M_PI*t/Tp) * ddf;
	    divterms = Gamma * pow(*H,4.) * dH * dH
	         * ((1.0/r) * (*H) * dH + 5.0 * dH * dH + 3.0 * (*H) * ddH);
            Mc = (2.0*M_PI*Cp/Tp) * cos(2.0*M_PI*t/Tp) * f - Ms - divterms;
	  } else {
	    Mc = 0.0;
	  }

	  /* actual calculate total M */
	  *M = Ms + Mc;
  }
  return 0;
}

struct TestABCDParameters exactD(const double t, const double r) {
  struct TestABCDParameters result;

  result.error_code = exactD_old(t, r, &result.H, &result.M);

  return result;
}
