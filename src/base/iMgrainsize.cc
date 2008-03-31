// Copyright (C) 2004-2008 Jed Brown and Ed Bueler
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

#include <petscda.h>
#include "iceModel.hh"


//! Compute a grain size from a pseudo-age, which is determined only by the vertical velocity component.
/*!
PISM allows the choice of the Goldsby-Kohlstedt flow law with option <tt>-gk</tt>.  
That flow law requires a grain size to compute the softness/viscosity.  To determine 
the grain size we use the Vostok core \lo\cite{VostokCore}\elo as a source for a 
universal relation between the age of the ice and its grain size; see grainSizeVostok().

By default we do not use the full model age, which takes a very long time
to equilibriate.  (If you want to use the full model age for this purpose add 
option <tt>-real_age_grainsize</tt>.)  Instead we compute a pseudo age which uses only 
the vertical (scalar) component of the velocity field and makes a steady state assumption.

In fact we solve this equation in each ice column:
     \f[w\frac{\partial a}{\partial z} \stackrel{\ast}{=} 1,\f]
where \f$a\f$ is the pseudo-age.  This represents a major simplification of 
the actual age equation
    \f[ \frac{\partial \tau}{\partial t} + u \frac{\partial \tau}{\partial x}
        + v \frac{\partial \tau}{\partial y} + w \frac{\partial \tau}{\partial z} = 1\f]
which is solved in ageStep().  There are two simplifications, first that we 
assume steady state for the age field, and secondly that the horizontal velocity 
is assumed to be zero.  (Alternately we could explain dropping the
horizontal advection terms by assuming the age field does not vary in horizontal 
directions; each column is assumed to have the same age profile as its neighbor.)

The boundary value for the first order hyperbolic equation \f$\ast\f$ is 
\f$a(z=H)=0\f$, where \f$H\f$ is the elevation of the surface of the ice.  We 
work down the column from the top, using the downward velocity to add
to the age.  If at any point the vertical velocity is positive then we assume 
the ice is old below that level; in the Vostok time scale, for producing a 
grain size, ``old'' means \f$10^6\f$ years; see grainSizeVostok().

The numerical method is to approximate \f$\ast\f$ by
    \f[\left(\frac{w_k + w_{k+1}}{2}\right)\,
         \left(\frac{a_{k+1} - a_k}{z_{k+1}-z_k}\right) = 1.\f]
or
    \f[a_k = a_{k+1} - \frac{2(z_{k+1}-z_k)}{w_k + w_{k+1}}.\f]
This has second order truncation error (at \f$z_{k+1/2}\f$) whether or not 
vertical grid is equally-spaced.
 */
PetscErrorCode  IceModel::computeGrainSize_PseudoAge(
                     const PetscScalar H, const PetscInt Mz, 
                     PetscScalar *w, PetscScalar *age_wspace,
                     PetscScalar **gs) {
  // don't call this method when realAgeForGrainSize == PETSC_TRUE
  PetscScalar *age = age_wspace;
  
  const PetscScalar  old_ice = 1.0e6 * secpera; // A million years
  const PetscInt     ks = grid.kBelowHeight(H);
  bool               downward = true;
  age[Mz-1] = 0.0;  // top always new
  for (PetscInt k = Mz-2; k >= 0; k--) {
    if (k+1 >= ks) {          // if either z_k or z_k+1 are at top of ice
      age[k] = 0.0;
    } else if (!downward) {   // once upward vel has been set, rest of column is old
      age[k] = old_ice;
    } else if (w[k] >= 0.0) { // upward (non-downward) velocity found; from now on will be old
      age[k] = old_ice;
      downward = false;
    } else { // at this point we know:   z_k < H,   z_k+1 < H,   w_k < 0,   w_k+1 < 0
      // implement a_k = a_{k+1} - \frac{2(z_{k+1}-z_k)}{w_k + w_{k+1}}
      const PetscScalar dz = grid.zlevels[k+1] - grid.zlevels[k];
      age[k] = age[k+1] - (2.0 * dz) / (w[k] + w[k+1]);
    }
  }
      
  // convert age or pseudo-age to grainsize and put in gs3
  for (PetscInt k = 0; k < Mz; k++) {
    (*gs)[k] = grainSizeVostok(age[k]);
  }

  return 0;
}


//! Use the Vostok core as a source of a relationship between the age of the ice and the grain size.
/*! 
A data set is interpolated.  The data is from \lo\cite{VostokCore}\elo.
 */
PetscScalar IceModel::grainSizeVostok(PetscScalar age) const {
  const PetscInt numPoints = 22;
  const PetscScalar ageAt[numPoints] = {  // ages in ka
    0.0000e+00, 5.0000e+01, 1.0000e+02, 1.2500e+02, 1.5000e+02,
    1.5800e+02, 1.6500e+02, 1.7000e+02, 1.8000e+02, 1.8800e+02,
    2.0000e+02, 2.2500e+02, 2.4500e+02, 2.6000e+02, 3.0000e+02,
    3.2000e+02, 3.5000e+02, 4.0000e+02, 5.0000e+02, 6.0000e+02,
    8.0000e+02, 1.0000e+04 };
  const PetscScalar gsAt[numPoints] = {   // grain sizes in m
    1.8000e-03, 2.2000e-03, 3.0000e-03, 4.0000e-03, 4.3000e-03,
    3.0000e-03, 3.0000e-03, 4.6000e-03, 3.4000e-03, 3.3000e-03,
    5.9000e-03, 6.2000e-03, 5.4000e-03, 6.8000e-03, 3.5000e-03,
    6.0000e-03, 8.0000e-03, 8.3000e-03, 3.6000e-03, 3.8000e-03,
    9.5000e-03, 1.0000e-02 };
  const PetscScalar a = age * 1.0e-3 / secpera; // Age in ka
  PetscInt l = 0;               // Left end of the binary search
  PetscInt r = numPoints - 1;   // Right end

  // If we are out of range
  if (a < ageAt[l]) {
    return gsAt[l];
  } else if (a > ageAt[r]) {
    return gsAt[r];
  }
  // Binary search for the interval
  while (r > l + 1) {
    const PetscInt j = (r + l) / 2;
    if (a < ageAt[j]) {
      r = j;
    } else {
      l = j;
    }
  }
  if ((r == l) || (PetscAbsReal(r - l) > 1)) {
    PetscPrintf(grid.com, "binary search in grainSizeVostok: oops.\n");
  }
  // Linear interpolation on the interval
  return gsAt[l] + (a - ageAt[l]) * (gsAt[r] - gsAt[l]) / (ageAt[r] - ageAt[l]);
}

