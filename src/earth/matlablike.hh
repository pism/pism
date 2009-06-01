/* Copyright (C) 2004-2009 Ed Bueler

 This file is part of PISM.

 PISM is free software; you can redistribute it and/or modify it under the
 terms of the GNU General Public License as published by the Free Software
 Foundation; either version 2 of the License, or (at your option) any later
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

#include <petscvec.h>
#include "cubature.h"  /* for 'integrand' type */

/*! Functions like Matlab's conv2(A,B,'same').

Only works on sequential PETSc Vecs.

Works by naively (sans FFT) convolving two sequential Vecs A,B.  Size of A is mA x nA.
Size of B is mB x nB.  Returns a result Vec which is the same size as A.
Effectively pads by zero, in that the result matches the discrete but
infinite case .  (The infinite case is where A(i,j) and B(i,j) are defined for
-inf < i,j < inf but A(i,j)=0 if i<0 or i>mA-1 or j<0 or j>nA-1
and B(i,j)=0 if i<0 or i>mB-1 or j<0 or j>nB-1.)

This operation is O(mA^2 nA^2), but an alternate FFT implementation would
be O(mA nA log(mA nA)), presumably.
 */
PetscErrorCode conv2_same(Vec vA, int mA, int nA,  Vec vB, int mB, int nB,
                          Vec &vresult);


/*! Functions like Matlab's interp1.  Wrapper for GSL 1D interpolation.

Compare Matlab's \code  yi=interp1(x,Y,xi,'linear') \endcode
Also \code yi=interp1(x,Y,xi,'linear','extrap') \endcode.

Does no input checking.
 */
double interp1_linear(double* x, double* Y, int N, double xi);


/*! Functions like Matlab's dblquad() for 2D integrals.  Wrapper for cubature by Steven Johnson.

The cubature code is used for double integration.  It replaces iterated use of
a single-variable adaptive rule as in Matlab's dblquad.  cubature.{h|c} by Johnson
is based on HIntLib by Rudolf Schuerer as well on GSL itself.  cubature is GPLed.
 */
double dblquad_cubature(integrand f, double ax, double bx, double ay, double by,
                        double reqRelError, void *fdata);


#endif // ifndef __matlablike_hh

