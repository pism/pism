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
! test L, for example.
!
! The velocity solution has
!    u(r) = \alpha(r) \hat r + w(r,z) \hat z
!
! A supporting preprint is in preparation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

int exactM(double r,
           double *alpha,
           const double EPS_ABS, const double EPS_REL, const int ode_method);
   /* input    : r in m,  0 <= r <= 700000 */
   /* output   : alpha in m s^-1; always positive */
   /* numerical: EPS_ABS
                 EPS_REL
                 ode_method */

#ifdef __cplusplus
}
#endif


#endif  /* __exactTestM_h */

