// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
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

// ELB moved grain size to separate source file 7/16/06: it is a physical model

PetscErrorCode  IceModel::updateGrainSizeIfNeeded() {
  // This is a front end to the grain size update system.  It initializes the
  // first time it is called and then performs an update no more often than
  // gsIntervalYears.
  PetscErrorCode ierr;
  static PetscScalar lastGSUpdateYear = grid.p->year;  // only happens when first called

  // If the current grain sizes are not expired, exit cleanly.
  if (grid.p->year - lastGSUpdateYear >= gsIntervalYears) {
    ierr = updateGrainSizeNow(); CHKERRQ(ierr);
    lastGSUpdateYear = grid.p->year;
    ierr = verbPrintf(2,grid.com, "g"); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com, "$"); CHKERRQ(ierr);
  }
  return 0;
}


void  IceModel::setConstantGrainSize(PetscScalar d) {
  gs3.setToConstant(d);
}


PetscErrorCode  IceModel::updateGrainSizeNow() {
  PetscErrorCode ierr;
  PetscScalar **H;
  PetscScalar *age, *gs, *w, *zz;

  // "age" here is a pseudo-age computed using vertical velocity; compare vtau which 
  //   is updated by temperatureStep()
  age = new PetscScalar[grid.p->Mz];
  gs = new PetscScalar[grid.p->Mz];
  w = new PetscScalar[grid.p->Mz];
  
  zz = new PetscScalar[grid.p->Mz];
  for (PetscInt k = 0; k < grid.p->Mz; k++) {
    zz[k] = k * grid.p->dz;
  }
  
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = gs3.needAccessToVals(); CHKERRQ(ierr);
  ierr = w3.needAccessToVals(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i < grid.xs + grid.xm; i++) {
    for (PetscInt j = grid.ys; j < grid.ys + grid.ym; j++) {
      // get column of w vals
      ierr = w3.getValColumn(i,j,grid.p->Mz,zz,w); CHKERRQ(ierr); 
   
      // compute age in column from vertical velocity
      const PetscInt ks = static_cast<PetscInt>(floor(H[i][j] / grid.p->dz));
      for (PetscInt k = grid.p->Mz-1; k >= 0; k--) {
        if (k >= ks-1) {        // At the top of the ice
          age[k] = 0.0;
        } else if (w[k] >= 0.0) { // Upward velocity
          age[k] = 1.0e6 * secpera; // A million years
        } else {
          // Solve this equation in each vertical interval:
          //       a_z = (1/w)     with boundary value a_{surface} = 0
          // a(z) = a(0) + \int_0^z \frac{dz}{w_0 + w' z}
          // a(z) = a(0) + (1/w') \log \frac{w_0 + w'z}{w_0}
          // a(z) = a(0) + (z/w_0) \frac{\log (1 + x)}{x}  where x = \frac{w' z}{w_0}
          // log(1 + x)/x = 1 - x/2 + x^2/3 - x^3/4 + ...
          const PetscScalar w_prime = (w[k] - w[k+1]) / grid.p->dz;
          const PetscScalar x = w_prime * grid.p->dz / w[k];
          // This is second order approximation since computing \frac{\log(1 + x)}{x}
          // has problems as x \to 0
          age[k] = age[k+1] - (grid.p->dz / w[k+1]) * (1.0 - x / 2.0);
        }
      }

      // convert to grainsize
      for (PetscInt k = 0; k < grid.p->Mz; k++) {
        gs[k] = grainSizeVostok(age[k]);
      }

      // put in gs3
      ierr = gs3.setValColumn(i,j,grid.p->Mz,zz,gs); CHKERRQ(ierr); 
      
    }
  }

  delete [] age;  delete [] gs;  delete [] w;  delete [] zz;

  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = gs3.doneAccessToVals(); CHKERRQ(ierr);
  ierr = w3.doneAccessToVals(); CHKERRQ(ierr);

  // velocitySIAStaggered() uses neighbor values for gs
  ierr = gs3.beginGhostComm(); CHKERRQ(ierr);
  ierr = gs3.endGhostComm(); CHKERRQ(ierr);

  return 0;
}


PetscScalar IceModel::grainSizeVostok(PetscScalar age) const {
  // This age/grain size relationship comes from a Vostok core
  const PetscInt numPoints = 22;
  const PetscScalar ageAt[numPoints] = {
    0.0000e+00, 5.0000e+01, 1.0000e+02, 1.2500e+02, 1.5000e+02,
    1.5800e+02, 1.6500e+02, 1.7000e+02, 1.8000e+02, 1.8800e+02,
    2.0000e+02, 2.2500e+02, 2.4500e+02, 2.6000e+02, 3.0000e+02,
    3.2000e+02, 3.5000e+02, 4.0000e+02, 5.0000e+02, 6.0000e+02,
    8.0000e+02, 1.0000e+04 };
  const PetscScalar gsAt[numPoints] = {
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

