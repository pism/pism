// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef __exactTestsFG_hh
#define __exactTestsFG_hh

/*
ELB 9/12/05;  05/12/06
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exactTestsFG is a C++ implementation of the exact solutions Test F & G for a
! thermocoupled ice sheet.  References:
!
!    Ed Bueler, Jed A. Kallen-Brown, and Craig S. Lingle, "Accuracy
!       analysis of a numerical scheme for thermocoupled ice sheets," in
!       preparation 2006
!
!    Ed Bueler and Jed A. Kallen-Brown, "An exact solution to a model of a
!       cold, shallow, and thermocoupled ice sheet," in preparation 2006
!
!    Ed Bueler, , Jed A. Kallen-Brown, David N. Covey, and
!       Latrice N. Bowman (2005) "Exact solutions and numerical verification
!       for isothermal ice sheets," J. Glaciol. 51 (no. 173) 291--306.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

int bothexact(const double t, const double r, const double z[], const int Mz,
              const double Cp, double &H, double &M, double TT[], double U[],
              double w[], double Sig[], double Sigc[]);

#endif
