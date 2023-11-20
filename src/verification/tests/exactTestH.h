/*
   Copyright (C) 2004-2006, 2016 Jed Brown and Ed Bueler

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

#ifndef __exactTestH_h
#define __exactTestH_h 1

#ifdef __cplusplus
extern "C"
{
#endif

/*
10/24/06
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exactTestH is a C implementation of a single isothermal exact solution
! which is a concatenation of Tests C and B from
!
!    Ed Bueler, Craig S. Lingle, Jed A. Kallen-Brown, David N. Covey, and
!       Latrice N. Bowman (2005) "Exact solutions and numerical verification
!       for isothermal ice sheets," J. Glaciol. 51 (no. 173), 291--306.
!
! Test H includes pointwise isostasy and was used in generating results in
!
!    Ed Bueler, Craig S. Lingle, and Jed Brown (2007) "Fast
!       computation of a deformable Earth model for ice-sheet simulations,"
!       Ann. Glaciol. 46, 97--105.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

struct TestHParameters {
  int error_code;
  double H;
  double M;
};

struct TestHParameters exactH(const double f, const double t, const double r);

#ifdef __cplusplus
}
#endif


#endif  /* __exactTestH_h */
