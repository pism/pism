/*
   Copyright (C) 2004-2008, 2014, 2015, 2016, 2023 Ed Bueler and Jed Brown and Constantine Khroulev

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

#include "pism/verification/tests/exactTestsFG.hh"
#include <cmath>      // for pow, exp, cos, sin, M_PI, fabs, sqrt
#include <cstddef>    // for size_t
#include <stdexcept>  // for runtime_error

namespace pism {

static double p3(double x) {
  // p_3=x^3-3*x^2+6*x-6, using Horner's
  return -6.0 + x*(6.0 + x*(-3.0 + x));
}

static double p4(double x) {
  // p_4=x^4-4*x^3+12*x^2-24*x+24, using Horner's
  return 24.0 + x*(-24.0 + x*(12.0 + x*(-4.0 + x)));
}

TestFGParameters exactFG(double t, double r, const std::vector<double> &z, double Cp) {

  const double SperA = 31556926.0;  // seconds per year; 365.2422 days

  // parameters describing extent of sheet:
  const double H0 = 3000.0;     // m
  const double L  = 750000.0;   // m

  // period of perturbation; inactive in Test F:
  const double Tp = 2000.0*SperA; // s

  // fundamental physical constants
  const double g    = 9.81;     // m / s^2; accel of gravity
  const double Rgas = 8.314;    // J / (mol K)

  // ice properties; parameters which appear in constitutive relation:
  const double rho    = 910.0;  // kg / m^3; density
  const double k      = 2.1;    // J / m K s; thermal conductivity
  const double cpheat = 2009.0; // J / kg K; specific heat capacity
  const double n      = 3.0;    // Glen exponent

  // next two are EISMINT II values; Paterson - Budd for T < 263
  const double A = 3.615E-13;   // Pa^-3 s^-1
  const double Q = 6.0E4;       // J / mol

  // EISMINT II temperature boundary condition (Experiment F):
  const double Ggeo  = 0.042;   // J / m^2 s; geo. heat flux
  const double ST    = 1.67E-5; // K m^-1
  const double Tmin  = 223.15;  // K
  const double Kcond = k / (rho*cpheat); // constant in temp equation

  const size_t Mz = z.size();
  TestFGParameters result(Mz);
  double
    &H    = result.H,
    &M    = result.M;
  std::vector<double>
    &T    = result.T,
    &U    = result.U,
    &w    = result.w,
    &Sig  = result.Sig,
    &Sigc = result.Sigc;

  // temporary storage
  std::vector<double> I3(Mz);

  if ((r <= 0) or (r >= L)) {
    throw std::runtime_error("exactFG(): code and derivation assume 0 < r < L  !");
  }

  // compute H from analytical steady state Hs (Test D) plus perturbation
  const double power = n / (2*n + 2);
  const double Hconst = H0 / pow(1 - 1 / n, power);
  const double s = r / L;
  const double lamhat = (1 + 1 / n)*s - (1 / n) + pow(1 - s, 1 + 1 / n) - pow(s, 1 + 1 / n);

  double f = 0.0;
  if ((r > 0.3*L) and (r < 0.9*L)) {
    f = pow(cos(M_PI*(r - 0.6*L) / (0.6*L)) , 2.0);
  } else {
    f = 0.0;
  }

  const double goft = Cp*sin(2.0*M_PI*t / Tp);

  H = Hconst*pow(lamhat, power) + goft*f;

  // compute T = temperature
  const double Ts = Tmin + ST*r;
  const double nusqrt = sqrt(1 + (4.0*H*Ggeo) / (k*Ts));
  const double nu = (k*Ts / (2.0*Ggeo))*(1 + nusqrt);
  for (size_t i = 0; i < Mz; i++) {
    if (z[i] < H) {
      T[i] = Ts * (nu + H) / (nu + z[i]);
    } else { // surface value above ice surface; matches numerical way
      T[i] = Ts;
    }
    // old way: extend formula above surface: T[i] = Ts * (nu + H) / (nu + z[i]);
  }

  // compute surface slope and horizontal velocity
  const double lamhatr = ((1 + 1 / n) / L)*(1 - pow(1 - s, 1 / n) - pow(s, 1 / n));

  double fr = 0.0;
  if ((r > 0.3*L) and (r < 0.9*L)) {
    fr =  - (M_PI / (0.6*L)) * sin(2.0*M_PI*(r - 0.6*L) / (0.6*L));
  } else {
    fr = 0.0;
  }

  const double Hr = Hconst * power * pow(lamhat, power - 1) * lamhatr + goft*fr;   // chain rule
  if (Hr > 0) {
    throw std::runtime_error("exactFG(): assumes H_r negative for all 0 < r < L !");
  }

  const double mu      = Q / (Rgas*Ts*(nu + H));
  const double surfArr = exp(-Q / (Rgas*Ts));
  const double Uconst  = 2.0 * pow(rho*g, n) * A;
  const double omega   = Uconst * pow( - Hr, n) * surfArr * pow(mu, -n - 1);

  for (size_t i = 0; i < Mz; i++) {
    if (z[i] < H) {
      I3[i] = p3(mu*H) * exp(mu*H) - p3(mu*(H - z[i])) * exp(mu*(H - z[i]));
      U[i] = omega * I3[i];
    } else { // surface value above ice surface; matches numerical way
      I3[i] = p3(mu*H) * exp(mu*H) - p3(0.0);  // z[i] = H case in above
      U[i] = omega * I3[i];
    }
  }

  // compute strain heating
  for (size_t i = 0; i < Mz; i++) {
    if (z[i] < H) {
      const double Sigmu =  - (Q*(nu + z[i])) / (Rgas*Ts*(nu + H));
      Sig[i] = (Uconst*g / cpheat) * exp(Sigmu) * pow(fabs(Hr)*(H - z[i]) , n + 1);
    } else {
      Sig[i] = 0.0;
    }
  }

  // compute vertical velocity
  const double lamhatrr = ((1 + 1 / n) / (n*L*L)) * (pow(1 - s, (1 / n) - 1) - pow(s, (1 / n) - 1));

  double frr = 0.0;
  if ((r > 0.3*L) and (r < 0.9*L)) {
    frr =  - (2.0*M_PI*M_PI / (0.36*L*L)) * cos(2.0*M_PI*(r - 0.6*L) / (0.6*L));
  } else {
    frr = 0.0;
  }

  const double Hrr = Hconst*power*(power - 1)*pow(lamhat, power - 2.0) * pow(lamhatr, 2.0) +
    Hconst*power*pow(lamhat, power - 1)*lamhatrr + goft*frr;
  const double Tsr = ST;
  const double nur = (k*Tsr / (2.0*Ggeo)) * (1 + nusqrt) + (1 / Ts) * (Hr*Ts - H*Tsr) / nusqrt;
  const double mur = ( - Q / (Rgas*Ts*Ts*pow(nu + H, 2.0))) * (Tsr*(nu + H) + Ts*(nur + Hr));
  const double phi = 1 / r + n*Hrr / Hr + Q*Tsr / (Rgas*Ts*Ts) - (n + 1)*mur / mu;   // division by r
  const double gam = pow(mu, n) * exp(mu*H) * (mur*H + mu*Hr) * pow(H, n);
  for (size_t i = 0; i < Mz; i++) {
    const double I4 = p4(mu*H) * exp(mu*H) - p4(mu*(H - z[i])) * exp(mu*(H - z[i]));
    w[i] = omega * ((mur / mu - phi)*I4 / mu + (phi*(H - z[i]) + Hr)*I3[i] - gam*z[i]);
  }

  // compute compensatory accumulation M
  const double I4H = p4(mu*H) * exp(mu*H) - 24.0;
  const double divQ =  - omega * (mur / mu - phi) * I4H / mu + omega * gam * H;
  const double Ht = (Cp*2.0*M_PI / Tp) * cos(2.0*M_PI*t / Tp) * f;
  M = Ht + divQ;

  // compute compensatory heating
  const double nut = Ht / nusqrt;
  for (size_t i = 0; i < Mz; i++) {
    if (z[i] < H) {
      const double dTt = Ts * ((nut + Ht)*(nu + z[i]) - (nu + H)*nut) * pow(nu + z[i], - 2.0);
      const double Tr = Tsr*(nu + H) / (nu + z[i])
        + Ts * ((nur + Hr)*(nu + z[i]) - (nu + H)*nur) * pow(nu + z[i], - 2.0);
      const double Tz =  - Ts * (nu + H) * pow(nu + z[i], - 2.0);
      const double Tzz = 2.0 * Ts * (nu + H) * pow(nu + z[i], - 3.0);
      Sigc[i] = dTt + U[i]*Tr + w[i]*Tz - Kcond*Tzz - Sig[i];
    } else {
      Sigc[i] = 0.0;
    }
  }

  return result;
}

} // end of namespace pism
