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

/* exit status could be one of these; return of zero indicates success */
#define TESTP_R_OUT_OF_RANGE 78465
#define TESTP_NOT_DONE       78466
#define TESTP_NOT_DECREASING 78467
#define TESTP_INVALID_METHOD 78468
#define TESTP_NO_LIST        78469

/* determines range on which W(r) is valid in Test P */
#define R0       25000.0       /* m */
#define L        0.9 * R0      /* m */

int exactP(double r, double *h, double *magvb, double *W,
           const double EPS_ABS, const double EPS_REL, const int ode_method);
   /* r in meters,  0 <= r <= L */
   /* returns h (m), magvb (m s-1), W (m) */

int exactP_list(double *r, int N, double *h, double *magvb, double *W,
           const double EPS_ABS, const double EPS_REL, const int ode_method);
   /* 1. assumes N values L >= r[0] > r[1] > ... > r[N-1] >= 0  (_decreasing_)
      2. assumes r, h, magvb, W are _allocated_ length N arrays  */

#ifdef __cplusplus
}
#endif

#endif  /* __exactTestP_h */

