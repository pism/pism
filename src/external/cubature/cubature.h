/*
 * Copyright (c) 2005 Steven G. Johnson
 *
 * Portions (see comments) based on HIntLib (also distributed under
 * the GNU GPL), copyright (c) 2002-2005 Rudolf Schuerer.
 *
 * Portions (see comments) based on GNU GSL (also distributed under
 * the GNU GPL), copyright (c) 1996-2000 Brian Gough.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef _cubature_h
#define _cubature_h 1

#ifdef __cplusplus
extern "C"
{
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>

/* Adaptive multidimensional integration on hypercubes (or, really,
   hyper-rectangles) using cubature rules.

   A cubature rule takes a function and a hypercube and evaluates
   the function at a small number of points, returning an estimate
   of the integral as well as an estimate of the error, and also
   a suggested dimension of the hypercube to subdivide.

   Given such a rule, the adaptive integration is simple:

   1) Evaluate the cubature rule on the hypercube(s).
      Stop if converged.

   2) Pick the hypercube with the largest estimated error,
      and divide it in two along the suggested dimension.

   3) Goto (1).

*/

typedef double (*integrand) (unsigned ndim, const double *x, void *);

/* Integrate the function f from xmin[dim] to xmax[dim], with at
   most maxEval function evaluations (0 for no limit),
   until the given absolute is achieved relative error.  val returns
   the integral, and estimated_error returns the estimate for the
   absolute error in val.  The return value of the function is 0
   on success and non-zero if there was an error. */
int adapt_integrate(integrand f, void *fdata,
            unsigned dim, const double *xmin, const double *xmax,
            unsigned maxEval,
            double reqAbsError, double reqRelError,
            double *val, double *estimated_error);

#ifdef __cplusplus
}
#endif

#endif /* ifndef _cubature_h */
