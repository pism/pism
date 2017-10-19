/* Copyright (C) 2004-2009, 2014, 2015, 2017 Ed Bueler and Constantine Khroulev

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
#include "pism/external/cubature/cubature.h"
#include "pism/util/petscwrappers/Vec.hh"

void conv2_same(Vec vA, int mA, int nA,  Vec vB, int mB, int nB,
                Vec vresult) {

  pism::petsc::VecArray2D
    A(vA, mA, nA),
    B(vB, mB, nB),
    result(vresult, mA, nA);

  for (int j=0; j < nA; j++) {
    for (int i=0; i < mA; i++) {
      double sum = 0.0;
      for (int r = std::max(0, i - mB + 1); r < std::min(mA, i); r++) {
        for (int s = std::max(0, j - nB + 1); s < std::min(nA, j); s++) {
          sum += A(r, s) * B(i - r, j - s);
        }
      }
      result(i, j) = sum;
    }
  }
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
  double   result = 0.0, estimated_error = 0.0;

  /* see cubature.h: */
  adapt_integrate(f, fdata, 2, xmin, xmax, 
                  maxEval, 0.0, reqRelError, &result, &estimated_error);
  return result;
}

