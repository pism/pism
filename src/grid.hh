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

#ifndef __grid_hh
#define __grid_hh

#include <petscda.h>
#include "pism_const.hh"

class IceParam {
public:
  char history[HISTORY_STRING_LENGTH]; // history of commands used to generate this file
  PetscScalar Lx, Ly;  // half width of the ice model grid in x-direction, y-direction (m)
  PetscScalar Lz, Lbz; // extent of the ice, bedrock in z-direction (m)
  PetscInt    Mx, My; // number of grid points in x-direction, y-direction
  PetscInt    Mz, Mbz; // number of grid points in z-direction (ice), z-direction (bedrock).
  PetscScalar dx, dy, dz; // spacing of grid
  PetscScalar year;       // current time; units of years
};

int initIceParam(MPI_Comm com, IceParam **param);

class IceGrid {
public:
  MPI_Comm    com;
  PetscMPIInt rank, size;
  IceParam    *p;
  DA          da2, da3, da3b;
  PetscInt    xs, xm, ys, ym;

  IceGrid(MPI_Comm c, PetscMPIInt r, PetscMPIInt s);
  IceGrid(MPI_Comm c, PetscMPIInt r, PetscMPIInt s, IceParam *p);
  ~IceGrid();
  PetscErrorCode createDA();
  PetscErrorCode setCoordinatesDA();
  PetscErrorCode destroyDA();
  PetscErrorCode viewDA();
  PetscErrorCode rescale(PetscScalar lx, PetscScalar ly, PetscScalar lz);
private:
  PetscTruth createDA_done;
};

#endif	/* __grid_hh */
