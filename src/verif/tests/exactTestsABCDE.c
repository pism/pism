/*
   Copyright (C) 2004-2006 Jed Brown and Ed Bueler
  
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
#include "exactTestsABCDE.h"

#define pi 3.14159265358979
#define SperA 31556926.0  /* seconds per year; 365.2422 days */

int exactA(const double r, double *H, double *M) {
  /* NOTE: t is in seconds */
  const double L = 750000.0;       /* m; distance of margin from center */
  const double M0 = 0.3 / SperA;   /* 30 cm/year constant accumulation */
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


int exactB(const double t, const double r, double *H, double *M) {
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
          * pow(  1.0-pow( pow(t/t0,-beta)*(r/R0), (n+1)/n ),  n/(2*n+1)  );
  else
    *H = 0.0;

  *M=0.0;
  return 0;
}


int exactC(const double t, const double r, double *H, double *M) {
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
          * pow(  1.0-pow( pow(t/t0,-beta)*(r/R0), (n+1)/n ),  n/(2*n+1)  );
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


int exactD(const double t, const double rin, double *H, double *M) {

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
	  if ((r>0.3*L) && (r<0.9*L))    f = pow( cos(pi*(r-0.6*L)/(0.6*L)) ,2.0);
	  else                           f = 0.0;
	  Hs = (H0 / pow(1-1/n,power)) * pow(lamhat,power);
	  *H = Hs + Cp * sin(2.0*pi*t/Tp) * f;

	  /* compute steady part of accumulation */
	  temp = pow(s,1/n) + pow(1-s,1/n) - 1;
	  C = Gamma * pow(H0,2*n+2) / pow(2.0 * L * (1-1/n),n);
	  Ms = (C/r) * pow(temp,n-1)
	                *  ( 2*pow(s,1/n) + pow(1-s,(1/n)-1) * (1 - 2*s) - 1 );

	  /* derivs of H */
	  if ((r>0.3*L) && (r<0.9*L)) {
	    df = -(pi/(0.6*L)) * sin(pi*(r-0.6*L)/(0.3*L));
	    ddf = -(pi*pi/(0.18*L*L)) * cos(pi*(r-0.6*L)/(0.3*L));
	    chi = (4.0/3.0)*s - 1.0/3.0 + pow(1-s,4.0/3.0) - pow(s,4.0/3.0);
	    dchi = -(4.0/(3.0*L)) * (pow(s,1.0/3.0) + pow(1-s,1.0/3.0) - 1);
	    ddchi = -(4.0/(9.0*L*L)) * ( pow(s,-2.0/3.0) - pow(1-s,-2.0/3.0) );
	    c1 = (3.0*H0) / (8.0*pow(2.0/3.0,3.0/8.0));
	    dHs = c1 * pow(chi,-5.0/8.0) * dchi;
	    ddHs = c1 * ( (-5.0/8.0) * pow(chi,-13.0/8.0) * dchi * dchi
	                   + pow(chi,-5.0/8.0) * ddchi );
	    dH = dHs + Cp * sin(2.0*pi*t/Tp) * df;
	    ddH = ddHs + Cp * sin(2.0*pi*t/Tp) * ddf;
	    divterms = Gamma * pow(*H,4.) * dH * dH
	         * ( (1.0/r) * (*H) * dH + 5.0 * dH * dH + 3.0 * (*H) * ddH );
	    Mc = (2.0*pi*Cp/Tp) * cos(2.0*pi*t/Tp) * f - Ms - divterms;
	  } else {
	    Mc = 0.0;
	  }

	  /* actual calculate total M */
	  *M = Ms + Mc;
  }
  return 0;
}


int exactE(const double xIN, const double yIN, 
           double *H, double *M, double *mu, double *ub, double *vb) {

  const double L = 750000.0;       /* m; distance of margin from center */
  const double M0 = 0.3 / SperA;   /* 30 cm/year constant accumulation */
  const double g = 9.81;           /* m/s^2; accel of gravity */
  const double rho = 910.0;        /* kg/m^3; density */
  const double n = 3.0;            /* Glen exponent */
  const double A = 1.0E-16/SperA;  /* = 3.17e-24  1/(Pa^3 s); */
                                   /* (EISMINT value) flow law parameter */
  const double Gamma = 2 * pow(rho * g,n) * A / (n+2);

  const double mu_max = 2.5e-11;         /* Pa^-1 m s^-1; max sliding coeff */
  const double r1 = 200e3, r2 = 700e3,   /* define region of sliding */
               theta1 = 10 * (pi/180), theta2 = 40 * (pi/180);
  const double rbot = (r2 - r1) * (r2 - r1),
               thetabot = (theta2 - theta1) * (theta2 - theta1);

  /* note all features are reflected across coordinate axes */
  double       x = fabs(xIN), y = fabs(yIN), 
               sgnx = 1.0, sgny = 1.0, r, theta;
  double       m, q, C, chi;
  double       mufactor, dchidr, P, dhdr, h_x, h_y, d2hdr2, dmudr, Mb;
  
  r = sqrt(x * x + y * y);
  if (xIN < 0)  sgnx = -1.0;
  if (yIN < 0)  sgny = -1.0;

  if (r < L) {
    m = 2.0 * n + 2.0;
    q = 1.0 + 1.0 / n;
    C = pow(pow(2.0, n - 1) * M0 / Gamma, 1.0 / m);
    chi = pow(L, q) - pow(r, q);
    *H = C * pow(chi, n / m); 

    if (x < 1.0)
      theta = pi / 2.0;
    else
      theta = atan(y / x);

    if ( (r <= r1) || (r >= r2) || (theta <= theta1) || (theta >= theta2) ) {
      /* if outside sliding region but within ice cap, return as in test A */
      *M = M0;
      *mu = 0.0;
      *ub = 0.0;
      *vb = 0.0;
    } else {
      /* if INSIDE sliding region */
      mufactor = mu_max * (4.0 * (theta - theta1) * (theta2 - theta) / thetabot);
      *mu = mufactor * (4.0 * (r - r1) * (r2 - r) / rbot);
 
      dchidr = -q * pow(r, 1.0/n);
      dhdr = ((n * C) / m) * pow(chi, (n - m) / m) * dchidr;
      /* also: dhdr = -(C / 2.0) * pow(r, 1.0 / n) * pow(chi, -(n + 2.0) / m) */
      P = rho * g * (*H);
      h_x = dhdr * cos(theta);
      h_y = dhdr * sin(theta);
      *ub = sgnx * ( -(*mu) * P * h_x );
      *vb = sgny * ( -(*mu) * P * h_y );
 
      d2hdr2 = dhdr * ( -((n + 2.0) / m) * dchidr / chi + 1.0 / (n * r) );
      dmudr = mufactor * 4.0 * (r1 + r2 - 2.0 * r) / rbot;
      Mb = -P * ( ((*mu) / r) * (*H) * dhdr + dmudr * (*H) * dhdr  
                  + 2.0 * (*mu) * dhdr * dhdr + (*mu) * (*H) * d2hdr2 );
      *M = M0 + Mb;
    }
  } else { /* outside of ice cap */
    *H = 0.0;
    *M = M0;
    *mu = 0.0;
    *ub = 0.0;
    *vb = 0.0;
  }
  return 0;
}

