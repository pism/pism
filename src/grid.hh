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

#include <petscbag.h>
#include <petscda.h>
#include "pism_const.hh"

class IceParam {
public:
  PetscErrorCode setFromOptions();
  char history[HISTORY_STRING_LENGTH];
  PetscScalar Lx, Ly, Lz, Lbz;
  PetscInt    Mx, My, Mz, Mbz;
  PetscScalar dx, dy, dz;
  PetscScalar year;              // years
};

int initIceParam(MPI_Comm com, IceParam **param, PetscBag *bag);

class IceGrid {
public:
  MPI_Comm    com;
  PetscMPIInt rank, size;
  IceParam    *p;
  PetscBag    bag;
  DA          da2, da3, da3b;
  PetscInt    xs, xm, ys, ym;

  IceGrid(MPI_Comm c, PetscMPIInt r, PetscMPIInt s);
  IceGrid(MPI_Comm c, PetscMPIInt r, PetscMPIInt s, IceParam *p, PetscBag b);
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
