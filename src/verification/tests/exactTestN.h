/*
   Copyright (C) 2010, 2014, 2016, 2023 Ed Bueler

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

#ifndef __exactTestN_h
#define __exactTestN_h 1

#ifdef __cplusplus
extern "C"
{
#endif

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exactTestN is a C implementation of the parabolic solution in
! Bodvardsson (1955), treated here as a manufactured exact solution to
! a steady-state SSA flow problem, including the mass continuity equation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

struct TestNConstants {
  double H0;    /* = dome thickness (m) */
  double L0;    /* = full flow-line length from dome to margin where H->0 (m) */
  double xc;    /* = in Bueler interpretation, the location of the calving front (m) */
  double a;     /* = surface mass balance lapse rate, with elevation (s-1) */
  double H_ela; /* = elevation of equilibrium line (m) */
  double k;     /* = coefficient for sliding */
  double H_xc;  /* = thickness at calving front */
  double T_xc;  /* = vertically-integrated longitudinal stress at calving front */
};

struct TestNConstants exactNConstants(void);

struct TestNParameters {
  int error_code;
  double H;               /* (m; ice thickness) */
  double h_x;             /* (; surface slope) */
  double u;               /* (m s-1; ice horizontal velocity) */
  double M;               /* (m s-1; surface mass balance) */
  double B;               /* (Pa s^(1/3); ice hardness) */
  double beta;            /* (Pa s m-1; linear sliding coefficient) */
};

struct TestNParameters exactN(double x);
   /* input    : x                   (m; 0.0 <= x <= L0)

      Assumes n = 3.

      In Bueler interpretation, M(x) and A(x) are constructed so that the
      solution in Bodvardsson (1955) can be thought of as solving mass continuity
      and SSA stress balance simultaneously:

         M(x) - (u H)_x = 0

         (2 H B(x) |u_x|^((1/n)-1) u_x)_x - beta(x) u = rho g H h_x

      Here H = H(x) is ice thickness and u = u(x) is ice velocity.  Also
      h(x) = H(x), up to a constant the user may choose, because the bed is flat.
      Following Bodvardsson, here is the equilibrium line altitude, surface
      mass balance, and the sliding coefficient:

         Hela = H0 / 1.5
         M(x) = a (h(x) - Hela)
         beta(x) = k rho g H(x)

      The boundary conditions are

         H(0) = H0

      and

         T(xc) = 0.5 (1 - rho/rhow) rho g H(xc)^2

      where T(x) is the vertically-integrated viscous stress,

         T(x) = 2 H(x) B(x) |u_x|^((1/n)-1) u_x.

      But B(x) is chosen so that this quantity is constant:  T(x) = T0.

      The boundary condition at x = xc implies that the calving front is
      exactly at the location where the ice sheet reaches flotation.

      return value =
         0 if successful
         1 if x < 0
         2 if x > L0

   */

#ifdef __cplusplus
}
#endif

#endif  /* __exactTestN_h */
