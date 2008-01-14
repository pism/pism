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

#ifndef __grid_hh
#define __grid_hh

#include <petscda.h>
#include "pism_const.hh"

//! Collects the parameters for the grid and computational domain.
/*! 
Note that the default choices made when constructing an instance of IceParam will be overridden 
by PISM runtime options, by the input file (if an input file is used), and frequently 
by derived classes (of IceModel).
 */
class IceParam {
public:
  IceParam();
  char history[HISTORY_STRING_LENGTH]; // history of commands used to generate this file

  PetscScalar Lx, Ly;  // half width of the ice model grid in x-direction, y-direction (m)
  PetscInt    Mx, My; // number of grid points in x-direction, y-direction
  PetscScalar dx, dy; // spacing of grid

  PetscScalar Lz, Lbz; // extent of the ice, bedrock in z-direction (m)
  PetscInt    Mz, Mbz; // number of grid points in z-direction (ice), z-direction (bedrock).
//  PetscScalar dx, dy, dz; // spacing of grid [PREVIOUSLY WAS equal spaced in vert]

  PetscScalar year;       // current time; units of years
  
protected:
  static const PetscScalar DEFAULT_ICEPARAM_Lx, DEFAULT_ICEPARAM_Ly, DEFAULT_ICEPARAM_Lz, 
                           DEFAULT_ICEPARAM_year;
  static const PetscInt    DEFAULT_ICEPARAM_Mx, DEFAULT_ICEPARAM_My, DEFAULT_ICEPARAM_Mz,
                           DEFAULT_ICEPARAM_Mbz;
};

//! Describes the PISM grid and the distribution of data across processors.
/*! This class creates and destroys the PETSc DAs (distributed arrays) which control
    how all PISM variables are distributed across multiple processors.  This class also
    contains the parameters which describe the PISM computational box and the 
    three-dimensional PISM finite difference grid.
 */
class IceGrid {
public:
  IceGrid(MPI_Comm c, PetscMPIInt r, PetscMPIInt s);
  IceGrid(MPI_Comm c, PetscMPIInt r, PetscMPIInt s, IceParam *p);
  ~IceGrid();

  PetscErrorCode createDA();
  PetscErrorCode destroyDA();
//  PetscErrorCode viewDA();

  PetscErrorCode rescale(const PetscScalar lx, const PetscScalar ly, 
                         const PetscScalar lz);
  PetscErrorCode rescale(const PetscScalar lx, const PetscScalar ly, 
                         const PetscScalar lz, const PetscTruth truelyPeriodic);
  PetscErrorCode rescale_using_zlevels(const PetscScalar lx, const PetscScalar ly, 
                                       const PetscTruth truelyPeriodic);

  bool        equalVertSpacing();
  PetscInt    kBelowHeight(const PetscScalar height);
  PetscInt    kBelowHeightEQ(const PetscScalar height);
  
  MPI_Comm    com;
  PetscMPIInt rank, size;
  IceParam    *p;
  DA          da2;
  PetscInt    xs, xm, ys, ym;

  PetscScalar *zlevels, *zblevels; // z levels, in ice & bedrock, internally stored and represented in 3d Vecs

  PetscScalar dzEQ, *zlevelsEQ;    // for equal spacing on [0,p->Lz]
  PetscScalar dzbEQ, *zblevelsEQ;  //  "    "      "    "  [-p->Lbz,0]

private:
  PetscTruth      createDA_done, equally_spaced;
  PetscErrorCode  setLevelsFromLsMs();
  bool            isIncreasing(const PetscInt len, PetscScalar *vals);
};

#endif	/* __grid_hh */
