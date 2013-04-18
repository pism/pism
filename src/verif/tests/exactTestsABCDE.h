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

#ifndef __exactTestsABCDE_h
#define __exactTestsABCDE_h 1

#ifdef __cplusplus
extern "C"
{
#endif

/*
ELB 5/12/06; 10/14/06; 10/24/06
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exactTestsABCDE is a C implementation of the isothermal exact solutions 
! (Tests) A, B, C, D, and E from:
!
!    Ed Bueler, Craig S. Lingle, Jed A. Kallen-Brown, David N. Covey, and
!       Latrice N. Bowman (2005) "Exact solutions and numerical verification
!       for isothermal ice sheets," J. Glaciol. 51 (no. 173) 291--306.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

int exactA(const double r, double *H, double *M);

int exactB(const double t, const double r, double *H, double *M);

int exactC(const double t, const double r, double *H, double *M);

int exactD(const double t, const double r, double *H, double *M);

int exactE(const double x, const double y, 
           double *H, double *M, double *mu, double *ub, double *vb);


#ifdef __cplusplus
}
#endif


#endif  /* __exactTestsABCDE_h */

