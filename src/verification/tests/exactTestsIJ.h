/*
   Copyright (C) 2004-2006, 2015 Jed Brown, Ed Bueler, and Constantine Khroulev

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

#ifndef __exactTestIJ_h
#define __exactTestIJ_h 1

#ifdef __cplusplus
extern "C"
{
#endif

/*
12/8/06; 8/24/07
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! exactTestIJ contains C implementations of:
! 1.  an exact solution for an ice stream sliding over plastic till described
!     on pages 237 and 238 of C. Schoof 2006 "A variational approach to ice
!     streams" J Fluid Mech 556 pp 227--251
! 2.  an exact solution for a linearized ice shelf with periodic boundary
!     conditions [CREATED BY ELB; ONLY REFERENCE IS EARLY PREPRINT]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

struct TestIParameters {
  double bed, tauc, u, v;
};

struct TestIParameters exactI(const double m, const double x, const double y);

struct TestJParameters {
  double H, nu, u, v;
};

struct TestJParameters exactJ(const double x, const double y);

#ifdef __cplusplus
}
#endif


#endif  /* __exactTestIJ_h */
