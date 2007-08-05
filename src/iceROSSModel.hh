// Copyright (C) 2006-2007 Ed Bueler
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

#ifndef __iceROSSModel_hh
#define __iceROSSModel_hh

#include <petscvec.h>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"

class IceROSSModel : public IceModel {
public:
    IceROSSModel(IceGrid &g, IceType &i);
    virtual PetscErrorCode initFromOptions();
    virtual PetscErrorCode  diagnosticRun();
    PetscErrorCode         finishROSS();

private:
    Vec             obsAzimuth, obsMagnitude, obsAccurate;    
    PetscErrorCode  createROSSVecs();
    PetscErrorCode  destroyROSSVecs();
    PetscErrorCode  fillinTemps();
    PetscErrorCode  readObservedVels(const char *fname);
    PetscErrorCode  putObservedVelsCartesian();
    PetscErrorCode  computeErrorsInAccurateRegion();
//    PetscErrorCode  readRIGGSandCompare();
};

#endif /* __iceROSSModel_hh */
