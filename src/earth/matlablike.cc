/* Copyright (C) 2004-2009, 2014 Ed Bueler

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

#include <math.h>
#include <gsl/gsl_spline.h>
#include <petscvec.h>
#include "cubature.h"
#include "error_handling.hh"


PetscErrorCode conv2_same(Vec vA, int mA, int nA,  Vec vB, int mB, int nB,
                          Vec &vresult) {
  PetscErrorCode  ierr;
  double sum, **A, **B, **result;
  int i,j,r,s;

  ierr = VecGetArray2d(vA, mA, nA, 0, 0, &A);
  PISM_PETSC_CHK(ierr, "VecGetArray2d");
  ierr = VecGetArray2d(vB, mB, nB, 0, 0, &B);
  PISM_PETSC_CHK(ierr, "VecGetArray2d");
  ierr = VecGetArray2d(vresult, mA, nA, 0, 0, &result);
  PISM_PETSC_CHK(ierr, "VecGetArray2d");
  for (i=0; i < mA; i++) {
    for (j=0; j < nA; j++) {
      sum = 0.0;
      for (r = std::max(0, i - mB + 1); r < std::min(mA, i); r++) {
        for (s = std::max(0, j - nB + 1); s < std::min(nA, j); s++) {
          sum += A[r][s] * B[i - r][j - s];
        }
      }
      result[i][j] = sum;
    }
  }
  ierr = VecRestoreArray2d(vA, mA, nA, 0, 0, &A);
  PISM_PETSC_CHK(ierr, "VecRestoreArray2d");
  ierr = VecRestoreArray2d(vB, mB, nB, 0, 0, &B);
  PISM_PETSC_CHK(ierr, "VecRestoreArray2d");
  ierr = VecRestoreArray2d(vresult, mA, nA, 0, 0, &result);
  PISM_PETSC_CHK(ierr, "VecRestoreArray2d");
  return 0;
}


double interp1_linear(double* x, double* Y, int N, double xi) {
  
  gsl_interp_accel* acc    = gsl_interp_accel_alloc();
  gsl_spline*       spline = gsl_spline_alloc(gsl_interp_linear,N);
  gsl_spline_init(spline,x,Y,N);

  double result = gsl_spline_eval(spline,xi,acc);

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  return result;
}


double dblquad_cubature(integrand f, double ax, double bx, double ay, double by,
                        double reqRelError, void *fdata) {

  double   xmin[2] = {ax, ay};
  double   xmax[2] = {bx, by};
  unsigned maxEval = 5000;
  double   val, estimated_error;

  /* see cubature.h: */
  adapt_integrate(f, fdata, 2, xmin, xmax, 
                  maxEval, 0.0, reqRelError, &val, &estimated_error);
  return val;
}

