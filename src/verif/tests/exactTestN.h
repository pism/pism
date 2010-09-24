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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exactTestN is a C implementation of a manufactured exact solution
! to the prognostic SSA flow problem, including the mass continuity equation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

int exactN(double t, double x,
           double *u, double *ux, double *h, double *b, double *hx,
           double *taud, double *taub, double *tau11,
           double *C);
   /* input    : t                   (s; t >= 0)
                 x                   (m; -400e3 <= x <= 400e3)

      output   : u = u(x)            (m s^-1; ice horizontal velocity)
                 ux = u_x(x)         (s^-1; longitudinal strain rate)
                 h = h(t,x)          (m; surface elevation)
                 b = b(x)            (m; bed elevation)
                 hx = h_x(t,x)       (; surface slope)
                 taud = taud(t,x)    (Pa; driving stress)
                 taub = taub(t,x)    (Pa; basal shear stress)
                 tau11 = tau11(t,x)  (Pa; longitudinal stress)
                 C = C(t,x)          (Pa s m-1; COMPENSATORY sliding coefficient
                                      to balance SSA)

      return value = 0 if successful
   */

#ifdef __cplusplus
}
#endif

#endif  /* __exactTestN_h */

