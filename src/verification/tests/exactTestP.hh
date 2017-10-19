/*
   Copyright (C) 2012-2013, 2016 Ed Bueler

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

#ifndef __exactTestP_h
#define __exactTestP_h 1

#include <vector>
#include <string>

namespace pism {

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exactTestP is a C++ implementation of a nearly-exact solution to the 'distributed'
! subglacial hydrology model described in the draft manuscript
!
!    Ed Bueler & Ward van Pelt (2013) "A distributed model of subglacial
!    and englacial hydrology in tidewater glaciers and ice sheets"
!
! This nearly-exact solution requires solving an ODE numerically.
! Only the steady water thickness solution W(r) is computed here.  The
! pressure P can be computed by the formula P(W) which applies in steady state.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/

/* determines range on which W(r) is valid in Test P */
#define TESTP_R0       25000.0         /* m */
#define TESTP_L        0.9 * TESTP_R0  /* m */

/* exit status of exactP() could be one of these; return of zero indicates success */
#define TESTP_R_NEGATIVE          78463
#define TESTP_W_EXCEEDS_WR        78464
#define TESTP_W_BELOW_WCRIT       78465
#define TESTP_INVALID_METHOD      78466
#define TESTP_NOT_DONE            78467
#define TESTP_NO_LIST             78482
#define TESTP_LIST_NOT_DECREASING 78483
#define TESTP_OLD_GSL             78484

struct TestPParameters {
  TestPParameters(int N)
    : r(N), h(N), magvb(N), Wcrit(N), W(N), P(N) {
    error_code = 0;
  }

  int error_code;
  std::string error_message;
  std::vector<double> r, h, magvb, Wcrit, W, P;
};


TestPParameters exactP(const std::vector<double> &r,
                       double EPS_ABS, double EPS_REL, int ode_method);
/* Input r in meters, assumes that values in are decreasing (r[0] > r[1] > ... > r[N-1] >= 0).
   ode_method = 1  : rk8pd is Runge-Kutta Prince-Dormand (8,9) [default]
                2  : rk2   is Runge-Kutta (2,3)
                3  : rkf45 is Runge-Kutta-Felberg (4,5)
                4  : rkck  is Runge-Kutta Cash-Karp (4,5)
   None of these are implicit.  The Jacobian has not been implemented.
   Returns h (m), magvb (m s-1), W_c (m), W (m), P (Pa). */

} // end of namespace pism

#endif  /* __exactTestP_h */
