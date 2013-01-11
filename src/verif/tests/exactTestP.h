/*
   Copyright (C) 2012 Ed Bueler

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

#ifndef __exactTestP_h
#define __exactTestP_h 1

#ifdef __cplusplus
extern "C"
{
#endif

/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! exactTestP is a C implementation of a nearly-exact solution to the 'distributed'
! subglacial hydrology model described in the draft manuscript
!
!    Ward van Pelt & Ed Bueler (2012) "A diffusive-closure model of
!    subglacial hydrology"
!
! This nearly-exact solution requires solving an ODE numerically.
! Only the steady water thickness solution W(r) is computed here.  The
! pressure P can be computed by the formula P(W) which applies in steady state.
! See simpleP.c.
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

int exactP(double r, double *h, double *magvb, double *Wcrit, double *W, double *P,
           const double EPS_ABS, const double EPS_REL, const int ode_method);
   /* Input r in meters,  0 <= r.
      ode_method = 1  : rk8pd is Runge-Kutta Prince-Dormand (8,9) [default]
                   2  : rk2   is Runge-Kutta (2,3)
                   3  : rkf45 is Runge-Kutta-Felberg (4,5)
                   4  : rkck  is Runge-Kutta Cash-Karp (4,5)
      None of these are implicit.  The Jacobian has not been implemented.
      Returns h (m), magvb (m s-1), W_c (m), W (m), P (Pa). */


/* exit status of exactP_list() could be one of the above or one of these or zero for success */
#define TESTP_NO_LIST             78482
#define TESTP_LIST_NOT_DECREASING 78483

int exactP_list(double *r, int N, double *h, double *magvb, double *Wcrit, double *W, double *P,
                const double EPS_ABS, const double EPS_REL, const int ode_method);
   /* 1. assumes N values r[0] > r[1] > ... > r[N-1] >= 0  (_decreasing_)
      2. assumes r, h, magvb, Wcrit, W are _allocated_ length N arrays  */

int error_message_testP(int status);

#ifdef __cplusplus
}
#endif

#endif  /* __exactTestP_h */

