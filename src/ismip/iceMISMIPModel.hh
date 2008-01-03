// Copyright (C) 2008 Ed Bueler
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

#ifndef __iceMISMIPModel_hh
#define __iceMISMIPModel_hh

#include <petscvec.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"

class MISMIPIce : public GlenIce {
public:
  MISMIPIce();
  virtual PetscScalar flow(const PetscScalar stress, const PetscScalar temp,
                           const PetscScalar pressure) const;
  // this one returns nu; ignors temp & pressure
  virtual PetscScalar effectiveViscosity(const PetscScalar regularization,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar temp, const PetscScalar pressure) const;
  // this one returns nu * H; calls effectiveViscosity() for nu
  virtual PetscScalar effectiveViscosityColumn(const PetscScalar regularization,
                           const PetscScalar H, const PetscScalar dz,
                           const PetscScalar u_x, const PetscScalar u_y,
                           const PetscScalar v_x, const PetscScalar v_y,
                           const PetscScalar *T1, const PetscScalar *T2) const;
  PetscErrorCode setA(const PetscScalar myA);
  PetscScalar    getA();
  PetscErrorCode printInfo(const int thresh, MPI_Comm com);
  
protected:
  PetscScalar  A_MISMIP, // softness
               B_MISMIP; // hardness; B = A^{-1/n}
};


class MISMIPBasalType : public ViscousBasalType {
public:
  MISMIPBasalType(const PetscScalar m, const PetscScalar C, const PetscScalar regularize);
  virtual PetscScalar    drag(PetscScalar coeff, PetscScalar tauc, PetscScalar vx, PetscScalar vy);
  virtual PetscErrorCode printInfo(const int thresh, MPI_Comm com);

protected:
  PetscScalar m_MISMIP, C_MISMIP;
  PetscScalar regularize_MISMIP;
};


class IceMISMIPModel : public IceModel {

public:
  IceMISMIPModel(IceGrid &g, IceType &i, MISMIPIce &mismip_i);
  virtual PetscErrorCode setFromOptions();
  virtual PetscErrorCode initFromOptions();
  PetscErrorCode         setMISMIPBed();
  PetscErrorCode         additionalAtEndTimestep();
  virtual PetscErrorCode summaryPrintLine(const PetscTruth printPrototype, const PetscTruth tempAndAge,
                           const PetscScalar year, const PetscScalar dt, 
                           const PetscScalar volume_kmcube, const PetscScalar area_kmsquare,
                           const PetscScalar meltfrac, const PetscScalar H0, const PetscScalar T0);

protected:
  MISMIPIce   &mismip_ice;
  PetscInt    exper, gridmode, runindex;
  char        sliding;
  PetscScalar runtimeyears;
  char        initials[PETSC_MAX_PATH_LEN];

};

#endif  // __iceMISMIPModel_hh

