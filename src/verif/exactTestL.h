/*
   Copyright (C) 2007 Ed Bueler
  
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

#ifndef __exactTestL_h
#define __exactTestL_h 1

#ifdef __cplusplus
extern "C"
{
#endif

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exactTestL is a C implementation of an isothermal "exact" solution on a
! no-flat bed described in section 2.3 of an incomplete preprint
!
!    Ed Bueler (March 2006) "Equilibrium ice sheets solve variational
!       inequalities"
!
! in this case the exact solution requires solving an ODE numerically
! (see src/exact/simpleL.c and src/exact/gridL.c)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

int exactL(double r, double *H, double *b, double *a, 
           const double EPS_ABS, const double EPS_REL, const int ode_method);
   /* r in meters,  0 <= r <= 750000 */
  
int exactL_list(double *r, int N, double *H, double *b, double *a);
   /* N values r[0] > r[1] > ... > r[N-1]  (_decreasing_)
      assumes r, H, b, a are _allocated_ length N arrays  */
   /* uses defaults EPS_ABS=1.0e-12, EPS_REL=0.0, ode_method=1 (RK Cash-Karp) */

#ifdef __cplusplus
}
#endif


#endif  /* __exactTestL_h */

