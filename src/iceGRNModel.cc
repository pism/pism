// Copyright (C) 2004-2007 Jed Brown and Ed Bueler
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <cstring>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"
#include "iceGRNModel.hh"

IceGRNModel::IceGRNModel(IceGrid &g, IceType &i) : IceModel(g, i) {
  // only call parent's constructor
}

void IceGRNModel::setflowlawNumber(PetscInt law) {
  flowlawNumber = law;
}

PetscInt IceGRNModel::getflowlawNumber() {
  return flowlawNumber;
}

PetscErrorCode IceGRNModel::setFromOptions() {
  PetscErrorCode ierr;
  
  ierr = IceModel::setFromOptions(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IceGRNModel::initFromOptions() {
  PetscErrorCode ierr;

  ierr = IceModel::initFromOptions(); CHKERRQ(ierr);

  char inFile[PETSC_MAX_PATH_LEN];
  PetscTruth inFileSet, bootFileSet;
  
  ierr = PetscOptionsGetString(PETSC_NULL, "-if", inFile,
                               PETSC_MAX_PATH_LEN, &inFileSet); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-bif", inFile,
                               PETSC_MAX_PATH_LEN, &bootFileSet); CHKERRQ(ierr);

  if (inFileSet == PETSC_TRUE) {
    //ierr = initFromFile(inFile); CHKERRQ(ierr);
  } else if (bootFileSet == PETSC_TRUE) {
    //ierr = bootstrapFromFile_netCDF(inFile);
    
    // after we set the new temperatures, we need
    // to set the 3D temps again
    ierr = initTs(); CHKERRQ(ierr);
    ierr = putTempAtDepth(); CHKERRQ(ierr);
    ierr = cleanExtraLand(); CHKERRQ(ierr);
  } else {
    SETERRQ(1, "Error: need input file.\n");
  }
                        
  if (!isInitialized()) {
    SETERRQ(1, "Model has not been initialized.\n");
  }       

  return 0;
}

PetscErrorCode IceGRNModel::additionalAtStartTimestep() {
    PetscErrorCode ierr;
 
    // at the beginning of each time step
    // we need to redo the model of the surface
    // temperatures.
    ierr = initTs(); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode IceGRNModel::initTs() {
  PetscErrorCode ierr;
  PetscScalar val, Z;
  PetscScalar **Ts, **lat, **h, **H, **b;
  
  ierr = DAVecGetArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      //EISMINT surface temperature model
      Z = PetscMax(h[i][j], 20 * (lat[i][j] - 65));
      val = 49.13 - 0.007992 * Z - 0.7576 * (lat[i][j]);
      Ts[i][j] = val + ice.meltingTemp;
    }
  }
  
  ierr = DAVecRestoreArray(grid.da2, vTs, &Ts); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vh, &h); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vH, &H); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(grid.da2, vbed, &b); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IceGRNModel::ellePiecewiseFunc(PetscScalar lon, PetscScalar *lat) {
  // first function
  float l1_x1 = -68.18, l1_y1 = 80.1;
  float l1_x2 = -62, l1_y2 = 82.24;

  // piecewise boundaries
  float m, b;

  m = (l1_y1 - l1_y2) / (l1_x1 - l1_x2);
  b = (l1_y2) - m * (l1_x2);
  *lat = m * lon + b;

  return 0;
}

PetscErrorCode IceGRNModel::cleanExtraLand(){
  PetscErrorCode ierr;
  PetscScalar lat_line;
  // remove mask SE of the following point
  float ice_lon = 30, ice_lat = 67;
  PetscScalar **lat, **lon, **mask;

  ierr = DAVecGetArray(grid.da2, vLatitude, &lat); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vLongitude, &lon); CHKERRQ(ierr);
  ierr = DAVecGetArray(grid.da2, vMask, &mask); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {

      ellePiecewiseFunc(lon[i][j], &lat_line);
      if (lat[i][j]>lat_line) {
          mask[i][j] = MASK_FLOATING_OCEAN0;
      } else if (lat[i][j] < ice_lat && lon[i][j] > -ice_lon) {
        mask[i][j] = MASK_FLOATING_OCEAN0;
      }
    }
  }

  return 0;
}
