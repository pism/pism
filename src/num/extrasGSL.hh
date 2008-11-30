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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

/* 
   extrasGSL.[hh|cc] includes procedures which use the GSL library for numerical 
   computations.  These are useful to the ice model but are not dependent on 
   any members of the class IceModel.  These procedures mimic Matlab functionality: 
      interp1_linear, dblquad_cubature.
   
   The cubature program by Steven Johnson is used for double integration, rather 
   than the iterated use of a single-variable adaptive rule as in Matlab's 
   dblquad.  cubature.[h|c] by Johnson is based on HIntLib by Rudolf Schuerer 
   as well on GSL itself.  Note cubature is distributed with GPL.
   
   The heapsort code below is not GSL.
*/

double interp1_linear(const double x[], const double Y[], int N, double xi);

// get 'integrand' type
#include "cubature.h"
double dblquad_cubature(integrand f, 
          const double ax, const double bx, const double ay, const double by,
          double reqRelError, void *fdata);

// heapsort from 
//   http://en.wikibooks.org/wiki/Algorithm_implementation/Sorting/Heapsort
// (found 7/7/07 by ELB;  see also http://en.wikipedia.org/wiki/Heapsort)
// modified to have two index arrays "follow along" and get rearranged the same way
//
// note it is possible I should have used a struct and qsort() from stdlib.h, but
// but what I have here seems to work well; also there are routines in gsl_sort.h
// and gsl_sort_vector.h
void heapsort_double_2indfollow(double arr[], int ia[], int ib[], unsigned int N);


