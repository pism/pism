/*
   Copyright (C) 2004-2006, 2014, 2016 Jed Brown and Ed Bueler and Constantine Khroulev

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

#ifndef EXACTTESTSFG_H
#define EXACTTESTSFG_H

#include <vector>

namespace pism {

/*
ELB 9/12/05;  05/12/06;  10/14/06;  5/30/08
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exactTestsFG is a C++ implementation of the exact solutions Test F & G for a
! thermocoupled ice sheet.  References:
!
!    Ed Bueler, Jed Brown, and Craig Lingle, "Exact solutions to the
!       thermomechanically coupled shallow ice approximation: effective
!       tools for verification,"  J. Glaciol. 53 (182), 499--516.
!
!    Ed Bueler and Jed Brown, "On exact solutions for cold, shallow, and
!       thermocoupled ice sheets," preprint arXiv:physics/0610106, 2006
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

struct TestFGParameters {
  TestFGParameters(int Mz)
    : T(Mz), U(Mz), w(Mz), Sig(Mz), Sigc(Mz) {
    // empty
  }

  double H, M;
  std::vector<double> T, U, w, Sig, Sigc;
};

TestFGParameters exactFG(double t, double r, const std::vector<double> &z, double Cp);

/*
 * NOTE:  Units returned for Sig and Sigc are K/s (i.e. temperature) not J/s.
 * This matches the published sources above but requires conversion in
 * PISM as of revision 311.
 */

} // end of namespace pism

#endif /* EXACTTESTSFG_H */
