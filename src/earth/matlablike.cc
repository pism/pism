/* Copyright (C) 2004-2009, 2014, 2015, 2017, 2018, 2023 Ed Bueler and Constantine Khroulev

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

#include "pism/earth/matlablike.hh"

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

