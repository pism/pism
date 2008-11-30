// Copyright (C) 2004-2008 Ed Bueler
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


/*
// original code from
//  http://en.wikibooks.org/wiki/Algorithm_implementation/Sorting/Heapsort
void heapsort(int arr[], unsigned int N)
{
   unsigned int n = N, i = n/2, parent, child;
   int t;

   for (;;) { // Loops until arr is sorted 
      if (i > 0) { // First stage - Sorting the heap 
            i--;           // Save its index to i 
            t = arr[i];    // Save parent value to t 
      } else {     // Second stage - Extracting elements in-place 
            n--;           // Make the new heap smaller 
            if (n == 0) return; // When the heap is empty, we are done 
            t = arr[n];    // Save last value (it will be overwritten) 
            arr[n] = arr[0]; // Save largest value at the end of arr 
      }

      parent = i; // We will start pushing down t from parent 
      child = i*2 + 1; // parent's left child 

      // Shift operation - pushing the value of t down the heap 
      while (child < n) {
            if (child + 1 < n  &&  arr[child + 1] > arr[child]) {
               child++; // Choose the largest child 
            }
            if (arr[child] > t) { // If any child is bigger than the parent 
               arr[parent] = arr[child]; // Move the largest child up 
               parent = child; // Move parent pointer to this child 
               child = parent*2 + 1; // Find the next child 
            } else {
               break; // t's place is found 
            }
      }
      arr[parent] = t; // We save t in the heap 
   }
}
*/

 
void heapsort_double_2indfollow(double arr[], int ia[], int ib[], unsigned int N) {
   unsigned int n = N, i = n/2, parent, child;
   double t;
   int tia, tib;

   for (;;) { /* Loops until arr is sorted */
      if (i > 0) { /* First stage - Sorting the heap */
            i--;           /* Save its index to i */
            t = arr[i];    /* Save parent value to t */
            tia = ia[i];
            tib = ib[i];
      } else {     /* Second stage - Extracting elements in-place */
            n--;           /* Make the new heap smaller */
            if (n == 0) return; /* When the heap is empty, we are done */
            t = arr[n];    /* Save last value (it will be overwritten) */
            tia = ia[n];
            tib = ib[n];
            arr[n] = arr[0]; /* Save largest value at the end of arr */
            ia[n] = ia[0];
            ib[n] = ib[0];
      }
      parent = i; /* We will start pushing down t from parent */
      child = i*2 + 1; /* parent's left child */
      while (child < n) { /* Shift operation - pushing the value of t down the heap */
            if (child + 1 < n  &&  arr[child + 1] > arr[child]) {
               child++; /* Choose the largest child */
            }
            if (arr[child] > t) { /* If any child is bigger than the parent */
               arr[parent] = arr[child]; /* Move the largest child up */
               ia[parent] = ia[child];
               ib[parent] = ib[child];
               parent = child; /* Move parent pointer to this child */
               child = parent*2 + 1; /* Find the next child */
            } else { break; /* t's place is found */ }
      }
      arr[parent] = t; /* We save t in the heap */
      ia[parent] = tia;
      ib[parent] = tib;
   }
}

