// Copyright (C) 2004-2009 Nathan Shemonski and Ed Bueler
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

#ifndef __iceGRNModel_hh
#define __iceGRNModel_hh

#include <petscvec.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"
#include "../base/forcing.hh"
#include "../coupler/pccoupler.hh"


class PISMEISGREENPDDCoupler : public PISMPDDCoupler {

protected:
  virtual PetscScalar getSummerWarming(
       const PetscScalar elevation, const PetscScalar latitude, const PetscScalar Tma);
};


class IceGRNModel : public IceModel {

public:
  IceGRNModel(IceGrid &g, IceType *i);
  virtual PetscErrorCode setFromOptions();
  using IceModel::initFromOptions;
  PetscErrorCode attachEISGREENPDDPCC(PISMEISGREENPDDCoupler &p);
  virtual PetscErrorCode initFromOptions(PetscTruth doHook = PETSC_TRUE);

protected:
  PISMEISGREENPDDCoupler *pddPCC; // points to same PISMAtmosCoupler as IceModel::atmosPCC,
                                  //   but we access PDD parameters through this pointer

  int expernum;  // SSL2 is 1, CCL3 is 3, GWL3 is 4
  virtual PetscErrorCode additionalAtStartTimestep();

private:
  PetscTruth     inFileSet, 
                 noEllesmereIcelandDelete,
                 haveSurfaceTemp,
                 haveGeothermalFlux;
  PetscScalar    calculateMeanAnnual(PetscScalar h, PetscScalar lat);
  PetscErrorCode updateTs();
  PetscErrorCode ellePiecewiseFunc(PetscScalar lon, PetscScalar *lat);
  PetscErrorCode cleanExtraLand();
};
#endif

