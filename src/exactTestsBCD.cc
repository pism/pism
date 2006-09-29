// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
                         
#include <cstdio>
#include <cmath>
#include "exactTestsBCD.hh"


const double pi = 3.14159265358979;
const double SperA=31556926.0;  // seconds per year; 365.2422 days


int exactB(double t, double r, double &H, double &M) {
  // NOTE: t and t0 are in seconds
  double lambda, alpha, beta, t0, Rmargin;
  const double n = 3.0, H0 = 3600.0, R0=750000.0;
  
  lambda=0.0;
  alpha=1.0/9.0;  // alpha=(2-(n+1)*lambda)/(5*n+3)=1/9
  beta=1.0/18.0;  // beta=(1+(2*n+1)*lambda)/(5*n+3)=1/18
  t0=422.45*SperA;  // t0=(beta/Gamma)*pow((2n+1)/(n+1),n)*(pow(R0,n+1)/pow(H0,2n+1))

  Rmargin=R0*pow(t/t0,beta);
  if (r<Rmargin)
    H = H0 * pow(t/t0,-alpha) * pow(  1.0-pow( pow(t/t0,-beta)*(r/R0), (n+1)/n ),  n/(2*n+1)  );
  else
    H = 0.0;

  M=0.0;

  return 0;
}


int exactC(double t, double r, double &H, double &M) {
  double lambda, alpha, beta, t0, Rmargin;
  const double n = 3.0, H0 = 3600.0, R0=750000.0;

  lambda=5.0;
  alpha=-1.0;  // alpha=(2-(n+1)*lambda)/(5*n+3)
  beta=2.0;  // beta=(1+(2*n+1)*lambda)/(5*n+3)
  t0=15208.0*SperA;  // t0=(beta/Gamma)*pow((2n+1)/(n+1),n)*(pow(R0,n+1)/pow(H0,2n+1))

  Rmargin=R0*pow(t/t0,beta);
  if (r<Rmargin)
    H = H0 * pow(t/t0,-alpha) * pow(  1.0-pow( pow(t/t0,-beta)*(r/R0), (n+1)/n ),  n/(2*n+1)  );
  else
    H = 0.0;

  if (t>0.1*SperA)
    M = (lambda/t)*H;
  else {  // when less than 0.1 year, avoid division by time
    Rmargin=R0*pow(0.1*SperA/t0,beta);
    if (r<Rmargin)
      M=5*H0/t0;  // constant value in disc of Rmargin radius
    else
      M=0.0;
  }
  return 0;
}


int exactD(double t, double r, double &H, double &M) {

  // parameters describing extent of sheet:
  const double H0=3600.0;          // m
  const double L=750000.0;         // m
  // parameters for perturbation:
  const double Tp=5000.0*SperA;    // s
  const double Cp=200.0;           // m
  // fundamental physical constants
  const double g=9.81;             // m/s^2; accel of gravity
  // ice properties; parameters which appear in constitutive relation:
  const double rho=910.0;          // kg/m^3; density
  const double n=3.0;              // Glen exponent
  const double A=1.0E-16/SperA;    // = 3.17e-24  1/(Pa^3 s); (EISMINT value) flow law parameter


  if (r < 0.0) r=-r;
 
  if (r >= L - 0.01) {
    H = 0.0;
    M = -0.1 / SperA;
  } else {
          if (r < 0.01) r = 0.01;  // avoid r=0 singularity in formulas
          
	  // important derived quantities
	  double Gamma = 2 * pow(rho * g,n) * A / (n+2);
	  double power = n / (2*n+2);
	  double s = r/L;

	  // compute H from analytical steady state Hs plus perturbation
	  double lamhat = (1+1/n)*s - (1/n) + pow(1-s,1+1/n) - pow(s,1+1/n);
	  double f;
	  if ((r>0.3*L) && (r<0.9*L))      f = pow( cos(pi*(r-0.6*L)/(0.6*L)) ,2);
	  else                             f = 0.0;
	  double Hs = (H0 / pow(1-1/n,power)) * pow(lamhat,power);
	  H = Hs + Cp * sin(2.0*pi*t/Tp) * f;

	  // compute steady part of accumulation
	  double temp = pow(s,1/n) + pow(1-s,1/n) - 1;
	  double C = Gamma * pow(H0,2*n+2) / pow(2.0 * L * (1-1/n),n);
	  double Ms = (C/r) * pow(temp,n-1) *  ( 2*pow(s,1/n) + pow(1-s,(1/n)-1) * (1 - 2*s) - 1 );

	  // derivs of H
	  double df, ddf, chi, dchi, ddchi, c1, dHs, ddHs, dH, ddH, divterms, Mc;
	  if ((r>0.3*L) && (r<0.9*L)) {
	    df = -(pi/(0.6*L)) * sin(pi*(r-0.6*L)/(0.3*L));
	    ddf = -(pi*pi/(0.18*L*L)) * cos(pi*(r-0.6*L)/(0.3*L));
	    chi = (4.0/3.0)*s - 1.0/3.0 + pow(1-s,4.0/3.0) - pow(s,4.0/3.0);
	    dchi = -(4.0/(3.0*L)) * (pow(s,1.0/3.0) + pow(1-s,1.0/3.0) - 1);
	    ddchi = -(4.0/(9.0*L*L)) * ( pow(s,-2.0/3.0) - pow(1-s,-2.0/3.0) );
	    c1 = (3.0*H0) / (8.0*pow(2.0/3.0,3.0/8.0));
	    dHs = c1 * pow(chi,-5.0/8.0) * dchi;
	    ddHs = c1 * ( (-5.0/8.0) * pow(chi,-13.0/8.0) * dchi * dchi + pow(chi,-5.0/8.0) * ddchi );
	    dH = dHs + Cp * sin(2.0*pi*t/Tp) * df;
	    ddH = ddHs + Cp * sin(2.0*pi*t/Tp) * ddf;
	    divterms= Gamma * pow(H,4.) * dH * dH * ( (1.0/r) * H * dH + 5.0 * dH * dH + 3.0 * H * ddH );
	    Mc = (2.0*pi*Cp/Tp) * cos(2.0*pi*t/Tp) * f - Ms - divterms;
	  } else {
	    Mc = 0.0;
	  }

	  // actual calculate total M
	  M = Ms + Mc;
  }
  return 0;
}
