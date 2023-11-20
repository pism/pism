/*
   Copyright (C) 2008, 2016 Ed Bueler

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

#ifndef __exactTestM_h
#define __exactTestM_h 1

#ifdef __cplusplus
extern "C"
{
#endif

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exactTestM is a C implementation of an isothermal "exact" solution
! to the diagnostic SSA flow problem for a constant thickness annular
! ice shelf, with calving front, attached to a grounded sheet and with
! Dirichlet (prescribed velocity) condition at the grounding line.  A
! first order ODE in the radial coordinate is solved numerically, as with
! test L, for example, so the solution is not exactly exact.
!
! The velocity solution has
!    u(r) = alpha(r) \hat r + w(r,z) \hat z
! alpha(r) is found for R_g=300km < r < R_c=600km.  For r < R_g, a smoothly
! decreasing to zero value is returned for alpha, while for r > R_c
! alpha(r)=0 is returned.
!
! A supporting preprint is in preparation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

struct TestMParameters {
  int error_code;         /* GSL_SUCCESS = 0 if successful */
  double alpha;           /* (m s^-1;  always positive) */
  double Drr;             /* = alpha'(r) (s^-1; radial strain rate) */
};

struct TestMParameters exactM(double r,
                              double EPS_ABS, double EPS_REL, int ode_method);
   /* input    : r                             (m;   r >= 0)
      numerical: EPS_ABS                       (=1.0e-12 recommended)
                 EPS_REL                       (=0.0     recommended)
                 ode_method                    (=1       recommended; =Runge-Kutta-Cash-Karp)
   */


#ifdef __cplusplus
}
#endif

#endif  /* __exactTestM_h */
