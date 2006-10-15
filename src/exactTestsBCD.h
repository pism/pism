
/*
   Copyright (C) 2004-2006 Jed Brown and Ed Bueler
  
   This file is part of Pism.
  
   Pism is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2 of the License, or (at your option) any later
   version.
  
   Pism is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
   details.
  
   You should have received a copy of the GNU General Public License
   along with Pism; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#ifndef __exactTestsBCD_h
#define __exactTestsBCD_h 1

#ifdef __cplusplus
extern "C"
{
#endif

#include <math.h>

/*
ELB 5/12/06; 10/14/06
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exactTestsBCD is a C implementation of the isothermal exact solutions 
! (Tests) B, C, and D from:
!
!    Ed Bueler, Craig S. Lingle, Jed A. Kallen-Brown, David N. Covey, and
!       Latrice N. Bowman (2005) "Exact solutions and numerical verification
!       for isothermal ice sheets," J. Glaciol. 51 (no. 173) 291--306.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

int exactB(double t, double r, double *H, double *M);

int exactC(double t, double r, double *H, double *M);

int exactD(double t, double r, double *H, double *M);


#ifdef __cplusplus
}
#endif


#endif

