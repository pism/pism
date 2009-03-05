// Copyright (C) 2009 Jed Brown
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

#ifndef _iceShelfExtension_hh
#define _iceShelfExtension_hh

#include <petsc.h>
class IceType;

//! Jed's shelfExtension object: a strength extension that factors the nu*H coefficient of the SSA equations so that it can use your IceType.
class IceShelfExtension {
public:
  IceShelfExtension(MPI_Comm,const char *pre);
  ~IceShelfExtension();
  PetscErrorCode setIce(const IceType *);
  PetscErrorCode setTemperature(PetscScalar);
  PetscErrorCode setThickness(PetscScalar);
  PetscErrorCode setStrainRate(PetscScalar);
  PetscErrorCode forceNuH(PetscScalar);
  PetscErrorCode setFromOptions();
  PetscErrorCode printInfo(PetscInt);
  PetscErrorCode view(PetscViewer);
  // Should this take a strain rate and just ignore it in the typical (linear extension) case.  The extension is
  // somewhat nonphysical anyway and using a power-law relation will effectively make it stronger than necessary in the
  // low-strain limit.
  PetscScalar viscosity();
  PetscScalar thickness() const;
private:
  MPI_Comm comm;
  char prefix[256];
  PetscScalar T,H,Du,cached_viscosity,force_nuH;
  const IceType *ice;
  IceType *private_ice;
};

#endif /* _iceShelfExtension_hh */
