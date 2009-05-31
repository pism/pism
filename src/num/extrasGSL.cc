// Copyright (C) 2004-2009 Ed Bueler
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

#include <cmath>
#include <stdio.h>
#include <gsl/gsl_spline.h>
#include "cubature.h"
#include "extrasGSL.hh"


double interp1_linear(const double x[], const double Y[], int N, double xi) {
// no-input-checking version of Matlab's
//     yi=interp1(x,Y,xi,'linear')  *or*  yi=interp1(x,Y,xi,'linear','extrap')
// invokes GSL's interpolation
  double result;
  
  gsl_interp_accel* acc = gsl_interp_accel_alloc();
  gsl_spline* spline = gsl_spline_alloc(gsl_interp_linear,N);
  gsl_spline_init(spline,x,Y,N);

  result = gsl_spline_eval(spline,xi,acc);

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  return result;
}


double dblquad_cubature(integrand f, 
           const double ax, const double bx, const double ay, const double by,
           double reqRelError, void *fdata) {

  double *xmin, *xmax;
  
  xmin = new double[2];
  xmin[0] = ax; 
  xmin[1] = ay;
  xmax = new double[2];
  xmax[0] = bx; 
  xmax[1] = by;
  const unsigned maxEval = 5000;
  double val, estimated_error;

  // see cubature.h:
  adapt_integrate(f, fdata, 2, (double*) xmin, (double*) xmax, 
                  maxEval, 0.0, reqRelError, &val, &estimated_error);

  delete [] xmin;
  delete [] xmax;
  return val;
}

