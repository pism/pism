/*
   Copyright (C) 2012-2017 Ed Bueler and Constantine Khroulev

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

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

#include <gsl/gsl_version.h>
#if (defined GSL_MAJOR_VERSION) && (defined GSL_MINOR_VERSION) && \
  ((GSL_MAJOR_VERSION >= 1 && GSL_MINOR_VERSION >= 15) || (GSL_MAJOR_VERSION >= 2))
#define PISM_USE_ODEIV2 1
#include <gsl/gsl_odeiv2.h>
#endif

#include "exactTestP.hh"

namespace pism {

static const double SperA = 31556926.0; // seconds per year; 365.2422 days
static const double g     = 9.81;       // m s-2
static const double rhoi  = 910.0;      // kg m-3
static const double rhow  = 1000.0;     // kg m-3

// major model parameters:
static const double Aglen = 3.1689e-24;          // Pa-3 s-1
static const double k     = (0.01 / (rhow * g)); // FIXME:  this is extremely low but it matches what we were using
static const double Wr    = 1.0;                 // m
static const double c1    = 0.500;               // m-1
static const double c2    = 0.040;               // [pure]

// specific to exact solution
static const double m0 = ((0.20 / SperA) * rhow); // kg m-2 s-1; = 20 cm year-1
static const double h0 = 500.0;                   // m
static const double v0 = (100.0 / SperA);         // m s-1
static const double R1 = 5000.0;                  // m

int getsb(double r, double *sb, double *dsbdr) {
  double CC;
  if (r < R1) {
    *sb    = 0.0;
    *dsbdr = 0.0;
  } else {
    CC     = pow((c1 * v0) / (c2 * Aglen * pow((TESTP_L - R1), 5.0)), (1.0 / 3.0));
    *sb    = CC * pow(r - R1, (5.0 / 3.0));
    *dsbdr = (5.0 / 3.0) * CC * pow(r - R1, (2.0 / 3.0));
  }
  return 0;
}


double criticalW(double r) {
  double
    h  = h0 * (1.0 - (r / TESTP_R0) * (r / TESTP_R0)),
    Po = rhoi * g * h;
  double sb, dsb;
  getsb(r, &sb, &dsb);

  double sbcube = sb * sb * sb;
  double Pocube = Po * Po * Po;

  return (sbcube / (sbcube + Pocube)) * Wr;
}


int funcP(double r, const double W[], double f[], void *params) {
  /* Computes RHS f(r,W) for differential equation as given in dampnotes.pdf
  at https://github.com/bueler/hydrolakes:
      dW
      -- = f(r,W)
      dr
  Compare doublediff.m.  Assumes Glen power n=3.
  */

  double sb, dsb, dPo, tmp1, omega0, numer, denom;

  (void)params; /* quash warning "unused parameters" */

  if (r < 0.0) {
    f[0] = 0.0; /* place-holder */
    return TESTP_R_NEGATIVE;
  } else if (r > TESTP_L) {
    f[0] = 0.0;
    return GSL_SUCCESS;
  } else {
    getsb(r, &sb, &dsb);
    omega0 = m0 / (2.0 * rhow * k);
    dPo    = -(2.0 * rhoi * g * h0 / (TESTP_R0 * TESTP_R0)) * r;
    tmp1   = pow(W[0], 4.0 / 3.0) * pow(Wr - W[0], 2.0 / 3.0);
    numer  = dsb * W[0] * (Wr - W[0]);
    numer  = numer - (omega0 * r / W[0] + dPo) * tmp1;
    denom  = (1.0 / 3.0) * sb * Wr + rhow * g * tmp1;
    f[0]   = numer / denom;
    return GSL_SUCCESS;
  }
}


/* Computes initial condition W(r=L) = W_c(L^-). */
double initialconditionW() {
  double hL, PoL, sbL;
  hL  = h0 * (1.0 - (TESTP_L / TESTP_R0) * (TESTP_L / TESTP_R0));
  PoL = rhoi * g * hL;
  sbL = pow(c1 * v0 / (c2 * Aglen), 1.0 / 3.0);
  return (pow(sbL, 3.0) / (pow(sbL, 3.0) + pow(PoL, 3.0))) * Wr;
}


double psteady(double W, double magvb, double Po) {
  double sbcube, frac, P;
  sbcube = c1 * fabs(magvb) / (c2 * Aglen);
  frac   = (W < Wr) ? (Wr - W) / W : 0.0;
  P      = Po - pow(sbcube * frac, 1.0 / 3.0);
  if (P < 0.0) {
    P = 0.0;
  }
  return P;
}

#ifdef PISM_USE_ODEIV2

/* Solves ODE for W(r), the exact solution.  Input r[] and output W[] must be
allocated vectors of length N.  Input r[] must be decreasing.  The combination
EPS_ABS = 1e-12, EPS_REL=0.0, method = RK Dormand-Prince O(8)/O(9)
is believed for now to be predictable and accurate.  Note hstart is negative
so that the ODE solver does negative steps.  Assumes
   0 <= r[N-1] <= r[N-2] <= ... <= r[1] <= r[0] <= L.                            */
int getW(const double *r, int N, double *W, double EPS_ABS, double EPS_REL, int ode_method) {
  int count;
  int status = TESTP_NOT_DONE;
  double rr, hstart;
  const gsl_odeiv2_step_type *Tpossible[4];
  const gsl_odeiv2_step_type *T;
  gsl_odeiv2_system sys = { funcP, NULL, 1, NULL }; /* Jac-free method and no params */
  gsl_odeiv2_driver *d;

  /* setup for GSL ODE solver; these choices don't need Jacobian */
  Tpossible[0] = gsl_odeiv2_step_rk8pd;
  Tpossible[1] = gsl_odeiv2_step_rk2;
  Tpossible[2] = gsl_odeiv2_step_rkf45;
  Tpossible[3] = gsl_odeiv2_step_rkck;
  if ((ode_method > 0) && (ode_method < 5)) {
    T = Tpossible[ode_method - 1];
  } else {
    return TESTP_INVALID_METHOD;
  }

  hstart = -1000.0;
  d      = gsl_odeiv2_driver_alloc_y_new(&sys, T, hstart, EPS_ABS, EPS_REL);

  /* initial conditions: (r,W) = (L,W_c(L^-));  r decreases from L toward 0 */
  rr = TESTP_L;
  for (count = 0; count < N; count++) {
    /* except at start, use value at end of last interval as initial value for subinterval */
    W[count] = (count == 0) ? initialconditionW() : W[count - 1];
    while (rr > r[count]) {
      status = gsl_odeiv2_driver_apply(d, &rr, r[count], &(W[count]));
      if (status != GSL_SUCCESS) {
        break;
      }
      if (W[count] > Wr) {
        return TESTP_W_EXCEEDS_WR;
      } else if (W[count] < criticalW(r[count])) {
        return TESTP_W_BELOW_WCRIT;
      }
    }
  }

  gsl_odeiv2_driver_free(d);
  return status;
}

#else
int getW(const double *r, int N, double *W,
         double EPS_ABS, double EPS_REL, int ode_method) {
  (void) r;
  (void) EPS_ABS;
  (void) EPS_REL;
  (void) ode_method;

  for (int j = 0; j < N; ++j) {
    W[j] = 0;
  }
  return TESTP_OLD_GSL;
}
#endif

int exactP_list(const double *r, int N, double *h, double *magvb, double *Wcrit, double *W, double *P,
                double EPS_ABS, double EPS_REL, int ode_method) {

  int i, M, status;
  /* check first: we have a list, r is decreasing, r is in range [0,L) */
  if (N < 1) {
    return TESTP_NO_LIST;
  }

  for (i = 0; i < N; i++) {
    if ((i > 0) && (r[i] > r[i - 1])) {
      return TESTP_LIST_NOT_DECREASING;
    }
    if (r[i] < 0.0) {
      return TESTP_R_NEGATIVE;
    }
  }

  M = 0;
  while (r[M] > TESTP_L) {
    h[M]     = 0.0;
    magvb[M] = 0.0;
    Wcrit[M] = 0.0;
    W[M]     = 0.0;
    P[M]     = 0.0;
    M++;
  }

  for (i = M; i < N; i++) {
    h[i] = h0 * (1.0 - (r[i] / TESTP_R0) * (r[i] / TESTP_R0));

    if (r[i] > R1) {
      magvb[i] = v0 * pow((r[i] - R1) / (TESTP_L - R1), 5.0);
    } else {
      magvb[i] = 0.0;
    }

    Wcrit[i] = criticalW(r[i]);
  }

  status = getW(&(r[M]), N - M, &(W[M]), EPS_ABS, EPS_REL, ode_method);

  if (status) {
    for (i = M; i < N; i++) {
      P[i] = 0.0;
    }
    return status;
  } else {
    for (i = M; i < N; i++) {
      P[i] = psteady(W[i], magvb[i], rhoi * g * h[i]);
    }
    return 0;
  }
}

TestPParameters exactP(const std::vector<double> &r,
                       double EPS_ABS, double EPS_REL, int ode_method) {
  TestPParameters result(r.size());
  result.r = r;

  result.error_code = exactP_list(&r[0], r.size(),
                                  &result.h[0],
                                  &result.magvb[0],
                                  &result.Wcrit[0],
                                  &result.W[0],
                                  &result.P[0],
                                  EPS_ABS, EPS_REL, ode_method);

  switch (result.error_code) {
  case 0:
    result.error_message = "success";
    break;
  case TESTP_R_NEGATIVE:
    result.error_message = "error: r < 0";
    break;
  case TESTP_W_EXCEEDS_WR:
    result.error_message = "error: W > W_r";
    break;
  case TESTP_W_BELOW_WCRIT:
    result.error_message = "error: W < W_crit";
    break;
  case TESTP_INVALID_METHOD:
    result.error_message = "error: invalid choice for ODE method";
    break;
  case TESTP_NOT_DONE:
    result.error_message = "error: ODE integrator not done";
    break;
  case TESTP_NO_LIST:
    result.error_message = "error: no list of r values at input to exactP_list()";
    break;
  case TESTP_LIST_NOT_DECREASING:
    result.error_message = "error: input list of r values to exactP_list() is not decreasing";
    break;
  default:
    result.error_message = "unknown error code";
  }

  return result;
}

} // end of namespace pism
