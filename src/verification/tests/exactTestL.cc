/* Copyright (C) 2016, 2017, 2022, 2023 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "pism/verification/tests/exactTestL.hh"

#include <cstdlib>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>       /* M_PI */

#include <gsl/gsl_version.h>
#if (defined GSL_MAJOR_VERSION) && (defined GSL_MINOR_VERSION) && \
  ((GSL_MAJOR_VERSION >= 1 && GSL_MINOR_VERSION >= 15) || (GSL_MAJOR_VERSION >= 2))
#define PISM_USE_ODEIV2 1
#include <gsl/gsl_odeiv2.h>
#endif

#include "pism/util/error_handling.hh"
using pism::RuntimeError;

static const double SperA = 31556926.0;    /* seconds per year; 365.2422 days */

static const double L   = 750.0e3; /* m;    i.e. 750 km */
static const double b0  = 500.0;   /* m */
static const double z0  = 1.2;
static const double g   = 9.81;
static const double rho = 910.0;
static const double n   = 3.0;  /* Glen power */

int funcL(double r, const double u[], double f[], void *params) {
  /*
  RHS for differential equation:
      du                  5/8   / a_0  r  (L^2 - r^2) \ 1/3
      -- = - (8/3) b'(r) u    - |---------------------|
      dr                        \ 2 L^2 \tilde\Gamma  /
  */

  const double Lsqr = L * L;
  const double a0 = 0.3 / SperA;   /* m/s;  i.e. 0.3 m/a */
  const double A = 1.0e-16 / SperA;  /* = 3.17e-24  1/(Pa^3 s); EISMINT I flow law */
  const double Gamma = 2 * pow(rho * g,n) * A / (n+2);
  const double tilGamma = Gamma * pow(n,n) / (pow(2.0 * n + 2.0, n));
  const double C = a0 / (2.0 * Lsqr * tilGamma);

  if (params == NULL) {} /* quash warning "unused parameters" */

  if ((r >= 0.0) && (r <= L)) {
    const double freq = z0 * M_PI / L;
    const double bprime = b0 * freq * sin(freq * r);
    f[0] =  - (8.0/3.0) * bprime * pow(u[0], 5.0/8.0)
            - pow(C * r * (Lsqr - r * r), 1.0/3.0);
  } else {
    f[0] = 0.0;  /* no changes outside of defined interval */
  }
  return GSL_SUCCESS;
}

/* exit status could be one of these; returned zero indicates success */
#define TESTL_NOT_DONE       8966
#define TESTL_NOT_DECREASING 8967
#define TESTL_INVALID_METHOD 8968
#define TESTL_NO_LIST        8969

#ifdef PISM_USE_ODEIV2

/* combination EPS_ABS = 1e-12, EPS_REL=0.0, method = 1 = RK Cash-Karp
   is believed to be predictable and accurate */
int getU(const double *r, int N, double *u,
         const double EPS_ABS, const double EPS_REL, const int ode_method) {
   /* solves ODE for u(r)=H(r)^{8/3}, 0 <= r <= L, for test L
      r and u must be allocated vectors of length N; r[] must be decreasing */
   int k, count;
   int status = TESTL_NOT_DONE;
   double rr, hstart;
   const gsl_odeiv2_step_type* Tpossible[4];
   const gsl_odeiv2_step_type *T;
   gsl_odeiv2_system sys = {funcL, NULL, 1, NULL};  /* Jac-free method and no params */
   gsl_odeiv2_driver *d;

   /* setup for GSL ODE solver; these choices don't need Jacobian */
   Tpossible[0] = gsl_odeiv2_step_rk8pd;
   Tpossible[1] = gsl_odeiv2_step_rk2;
   Tpossible[2] = gsl_odeiv2_step_rkf45;
   Tpossible[3] = gsl_odeiv2_step_rkck;
   if ((ode_method > 0) && (ode_method < 5))
     T = Tpossible[ode_method-1];
   else {
     return TESTL_INVALID_METHOD;
   }

   /* check first: we have a list, and r is decreasing */
   if (N < 1) return TESTL_NO_LIST;
   for (k = 1; k<N; k++) {  if (r[k] > r[k-1]) return TESTL_NOT_DECREASING;  }

   /* outside of ice cap, u = 0 */
   k = 0;
   while (r[k] >= L) {
     u[k] = 0.0;
     k++;
     if (k == N) return GSL_SUCCESS;
   }

   /* initialize GSL ODE solver */
   hstart = -10000.0;
   d = gsl_odeiv2_driver_alloc_y_new(&sys, T, hstart, EPS_ABS, EPS_REL);

   /* initial conditions: (r,u) = (L,0);  r decreases from L */
   rr = L;
   for (count = k; count < N; count++) {
     /* except at start, use value at end of last interval as initial guess */
     u[count] = (count == 0) ? 0.0 : u[count-1];
     while (rr > r[count]) {
       status = gsl_odeiv2_driver_apply(d, &rr, r[count], &(u[count]));
       if (status != GSL_SUCCESS){
         break;
       }
     }
   }

   gsl_odeiv2_driver_free(d);
   return status;
}

#else  /* old GSL (< 1.15) */
int getU(const double *r, int N, double *u,
         const double EPS_ABS, const double EPS_REL, const int ode_method) {
  (void) r;
  (void) N;
  (void) u;
  (void) EPS_ABS;
  (void) EPS_REL;
  (void) ode_method;
  throw RuntimeError(PISM_ERROR_LOCATION, "Test L requires GSL version 1.15 or later.");
  return 0;
}
#endif

int exactL_list(const double *r, int N, double *H, double *b, double *a) {
  /* N values r[0] > r[1] > ... > r[N-1]  (decreasing)
     assumes r, H, b, a are allocated length N arrays  */

  const double Lsqr = L * L;
  const double a0 = 0.3 / SperA;   /* m/s;  i.e. 0.3 m/a */
  int stat, i;

  std::vector<double> u(N);

  /* combination EPS_ABS = 1e-12, EPS_REL=0.0, method = 1 = RK Cash-Karp
     believed to be predictable and accurate */
  stat = getU(r,N,u.data(),1.0e-12,0.0,1);
  if (stat != GSL_SUCCESS) {
    return stat;
  }

  for (i = 0; i < N; i++) {
    H[i] = pow(u[i],3.0/8.0);
    b[i] = - b0 * cos(z0 * M_PI * r[i] / L);
    a[i] = a0 * (1.0 - (2.0 * r[i] * r[i] / Lsqr));
  }

  return 0;
}

ExactLParameters::ExactLParameters(size_t size)
  : H(size), a(size), b(size) {
}

ExactLParameters exactL(const std::vector<double> &r) {
  ExactLParameters result(r.size());

  int error_code = exactL_list(&r[0], r.size(), &result.H[0], &result.b[0], &result.a[0]);

  switch (error_code) {
  case TESTL_NOT_DONE:
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Test L ERROR: exactL_list() returns 'NOT_DONE' (%d) ...",
                                  error_code);
    break;
  case TESTL_NOT_DECREASING:
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Test L ERROR: exactL_list() returns 'NOT_DECREASING' (%d) ...",
                                  error_code);
    break;
  case TESTL_INVALID_METHOD:
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Test L ERROR: exactL_list() returns 'INVALID_METHOD' (%d) ...",
                                  error_code);
    break;
  case TESTL_NO_LIST:
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Test L ERROR: exactL_list() returns 'NO_LIST' (%d) ...",
                                  error_code);
    break;
  default:
    break;
  }

  return result;
}
