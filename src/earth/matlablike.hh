/* Copyright (C) 2004-2009, 2015, 2017, 2018 Ed Bueler

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

#ifndef __matlablike_hh
#define __matlablike_hh

#include "pism/external/cubature/cubature.h"  /* for 'integrand' type */

/*! Functions like Matlab's dblquad() for 2D integrals.  Wrapper for cubature by Steven Johnson.

The cubature code is used for double integration.  It replaces iterated use of
a single-variable adaptive rule as in Matlab's dblquad.  cubature.{h|c} by Johnson
is based on HIntLib by Rudolf Schuerer as well on GSL itself.  cubature is GPLed.
 */
double dblquad_cubature(integrand f, double ax, double bx, double ay, double by,
                        double reqRelError, void *fdata);


#endif // ifndef __matlablike_hh
