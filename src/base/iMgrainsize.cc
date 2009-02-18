// Copyright (C) 2004-2009 Jed Brown and Craig Lingle and Ed Bueler
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


//! Use the Vostok core as a source of a relationship between the age of the ice and the grain size.
/*! 
A data set is interpolated here.  The intention is that the softness of the ice has
nontrivial dependence on its age, through its grainsize, because of variable dustiness
of the global climate.  The grainsize is partly determined by at which point in 
the glacial cycle the given ice fell as snow.

The data is from \lo\cite{DeLaChapelleEtAl98} and 
\cite{LipenkovEtAl89}\elo.  In particular, Figure A2 in the former reference was
hand-sampled with an attempt to include the ``wiggles'' in that figure.  Ages of
the oldest ice (>= 300 ka) were estimated in a necessarily ad hoc way.  The 
age value of 10000 ka was added simply to give interpolation for very old ice;
ages beyond that get constant extrapolation.  Linear interpolation is done between
the samples.
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

