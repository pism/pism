/*
   Copyright (C) 2010 Ed Bueler
  
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

int geometry_exactN(double *H0, double *L0, double *xc);
   /* output:  H0 = dome thickness (m)
               L0 = full flow-line length from dome to margin where H->0 (m)
               xc = in Bueler interpretation, the location of the calving front (m) */

int exactN(double x, double *h, double *hx, double *u, double *M, double *A);
   /* input    : x                   (m; 0.0 <= x <= L0)

      output   : h = h(x)            (m; surface elevation)
                 hx = h_x(x)         (; surface slope)
                 u = u(x)            (m s-1; ice horizontal velocity)
                 M = M(x)            (m s-1; surface mass balance)
                 A = A(x)            (Pa-3 s-1; ice softness)

      Assumes n = 3.
      
      In Bueler interpretation, M(x) and A(x) are constructed so that the
      solution in Bodvardsson (1955) can be thought of as solving mass continuity
      and SSA stress balance simultaneously:

         M(x) - (u H)_x = 0

         ( 2 H A(x)^(-1/n) |u_x|^((1/n)-1) u_x )_x - beta(x) u = rho g H h_x
         
      Here H = H(x) is ice thickness and u = u(x) is ice velocity.  Also
      h(x) = H(x) because the bed is flat and, following Bodvardsson,
         
         M(x) = a (h(x) - Hela),   Hela = H0 / 1.5
      
         beta(x) = k rho g H(x).
      
      The boundary conditions are
      
         H(0) = H0

      and
      
         T(xc) = 0.5 (1 - rho/rhow) rho g H(xc)^2

      where T(x) is the vertically-integrated viscous stress,
      
         T(x) = 2 H(x) A(x)^(-1/n) |u_x|^((1/n)-1) u_x.

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

