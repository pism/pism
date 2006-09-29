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
#include "exactTestsFG.hh"

double p3(double x) {
  // p_3=x^3-3*x^2+6*x-6, using Horner's
  return -6.0 + x*(6.0 + x*(-3.0 + x));
}

double p4(double x) {
  // p_4=x^4-4*x^3+12*x^2-24*x+24, using Horner's
  return 24.0 + x*(-24.0 + x*(12.0 + x*(-4.0 + x)));
}

int bothexact(const double t, const double r, const double z[], const int Mz,
              const double Cp, double &H, double &M, double TT[], double U[],
              double w[], double Sig[], double Sigc[]) {

  const double pi = 3.14159265358979;
  const double SperA=31556926.0;  // seconds per year; 365.2422 days

  // parameters describing extent of sheet:
  const double H0=3000.0;    // m
  const double L=750000.0;   // m
  // period of perturbation; inactive in Test F:
  const double Tp=2000.0*SperA;  // s

  // fundamental physical constants
  const double g=9.81;       // m/s^2; accel of gravity
  const double Rgas=8.314;   // J/(mol K)
  // ice properties; parameters which appear in constitutive relation:
  const double rho=910.0;    // kg/m^3; density
  const double k=2.1;        // J/m K s; thermal conductivity
  const double cpheat=2009.0;// J/kg K; specific heat capacity
  const double n=3.0;          // Glen exponent
  // next two are EISMINT II values; Paterson-Budd for T<263
  const double A=3.615E-13;  // Pa^-3 s^-1
  const double Q=6.0E4;      // J/mol
  // EISMINT II temperature boundary condition (Experiment F):
  const double Ggeo=0.042;   // J/m^2 s; geo. heat flux
  const double ST=1.67E-5;   // K m^-1
  const double Tmin=223.15;  // K
  const double Kcond=k/(rho*cpheat);  // constant in temp eqn

  // declare all temporary quantities; computed in blocks below
  double power, Hconst, s, lamhat, f, goft, Ts, nusqrt, nu;
  double lamhatr, fr, Hr, mu, surfArr, Uconst, omega;
  double Sigmu, lamhatrr, frr, Hrr, Tsr, nur, mur, phi, gamma;
  double I4H, divQ, Ht, nut;
  double I4,dTt,Tr,Tz,Tzz;
  double *I3;
  I3 = new double[Mz];  // only needed temporary array
  int i;

  if ( (r<=0) || (r>=L) ) {
    printf("\nERROR: code and derivation assume 0<r<L  !\n\n");
    return 1;
  }

  // compute H from analytical steady state Hs (Test D) plus perturbation
  power = n/(2*n+2);
  Hconst = H0/pow(1-1/n,power);
  s = r/L;
  lamhat = (1+1/n)*s - (1/n) + pow(1-s,1+1/n) - pow(s,1+1/n);
  if ((r>0.3*L) && (r<0.9*L))
    f = pow( cos(pi*(r-0.6*L)/(0.6*L)) ,2);
  else
    f = 0.0;
  goft = Cp*sin(2.0*pi*t/Tp);
  H = Hconst*pow(lamhat,power) + goft*f;

  // compute TT = temperature
  Ts = Tmin+ST*r;
  nusqrt = sqrt( 1 + (4.0*H*Ggeo)/(k*Ts) );
  nu = ( k*Ts/(2.0*Ggeo) )*( 1 + nusqrt );
  for (i=0; i<Mz; i++)
    TT[i] = Ts * (nu+H) / (nu+z[i]);

  // compute surface slope and horizontal velocity
  lamhatr = ((1+1/n)/L)*( 1 - pow(1-s,1/n) - pow(s,1/n) );
  if ( (r>0.3*L) && (r<0.9*L) )
    fr = -(pi/(0.6*L)) * sin(2.0*pi*(r-0.6*L)/(0.6*L));
  else
    fr = 0.0;
  Hr = Hconst * power * pow(lamhat,power-1) * lamhatr + goft*fr;   // chain rule
  if ( Hr>0 ) {
    printf("\nERROR: assumes H_r negative for all 0<r<L  !\n");
    return 1;
  }
  mu = Q/(Rgas*Ts*(nu+H));
  surfArr = exp(-Q/(Rgas*Ts));
  Uconst = 2.0 * pow(rho*g,n) * A;
  omega = Uconst * pow(-Hr,n) * surfArr * pow(mu,-n-1);
  for (i=0; i<Mz; i++) {
    I3[i] = p3(mu*H) * exp(mu*H) - p3(mu*(H-z[i])) * exp(mu*(H-z[i]));
    U[i] = omega * I3[i];
  }

  // compute strain heating
  for (i=0; i<Mz; i++) {
    Sigmu = -(Q*(nu+z[i])) / (Rgas*Ts*(nu+H));
    Sig[i] = (Uconst*g/cpheat) * exp(Sigmu) * pow( fabs(Hr)*( H -z[i]) ,n+1);
  }

  // compute vertical velocity
  lamhatrr = ((1+1/n) / (n*L*L)) * ( pow(1-s,(1/n)-1) - pow(s,(1/n)-1) );
  if ( (r>0.3*L) && (r<0.9*L) )
    frr = -(2.0*pi*pi/(0.36*L*L)) * cos(2.0*pi*(r-0.6*L)/(0.6*L));
  else
    frr = 0.0;
  Hrr = Hconst*power*(power-1)*pow(lamhat,power-2.0) * pow(lamhatr,2.0)  +
    Hconst*power*pow(lamhat,power-1)*lamhatrr + goft*frr;
  Tsr = ST;
  nur = (k*Tsr/(2.0*Ggeo)) * (1 + nusqrt) +
    (1/Ts) * (Hr*Ts-H*Tsr) / nusqrt;
  mur = ( -Q/(Rgas*Ts*Ts*pow(nu+H,2.0)) ) * ( Tsr*(nu+H)+Ts*(nur+Hr) );
  phi = 1/r + n*Hrr/Hr + Q*Tsr/(Rgas*Ts*Ts) - (n+1)*mur/mu;   // division by r
  gamma = pow(mu,n) * exp(mu*H) * (mur*H+mu*Hr) * pow(H,n);
  for (i=0; i<Mz; i++) {
    I4 = p4(mu*H) * exp(mu*H) - p4(mu*(H-z[i])) * exp(mu*(H-z[i]));
    w[i] = omega * ((mur/mu - phi)*I4/mu + (phi*(H-z[i])+Hr)*I3[i] - gamma*z[i]);
  }

  // compute compensatory accumulation M
  I4H = p4(mu*H) * exp(mu*H) - 24.0;
  divQ = - omega * (mur/mu - phi) * I4H / mu + omega * gamma * H;
  Ht = (Cp*2.0*pi/Tp) * cos(2.0*pi*t/Tp) * f;
  M = Ht + divQ;

  // compute compensatory heating
  nut = Ht/nusqrt;
  for (i=0; i<Mz; i++) {
    dTt = Ts * ((nut+Ht)*(nu+z[i])-(nu+H)*nut) * pow(nu+z[i],-2.0);
    Tr = Tsr*(nu+H)/(nu+z[i])
      + Ts * ((nur+Hr)*(nu+z[i])-(nu+H)*nur) * pow(nu+z[i],-2.0);
    Tz = -Ts * (nu+H) * pow(nu+z[i],-2.0);
    Tzz = 2.0 * Ts * (nu+H) * pow(nu+z[i],-3.0);
    Sigc[i] = dTt + U[i]*Tr + w[i]*Tz - Kcond*Tzz - Sig[i];
  }

  delete [] I3;
  return 0;
}
