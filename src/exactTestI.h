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

#ifndef __exactTestI_h
#define __exactTestI_h 1

#ifdef __cplusplus
extern "C"
{
#endif

/*
12/8/06
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! exactTestI is a C implementation of the exact solution for an ice stream
! sliding over plastic till described on pages 237 and 238 of
! C. Schoof 2006 "A variational approach to ice streams" J Fluid Mech 556 pp 227--251
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

int exactI(const double m, const double x, const double y, 
           double *bed, double *tauc, double *u, double *v);

#ifdef __cplusplus
}
#endif


#endif  /* __exactTestI_h */
