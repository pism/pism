// Copyright (C) 2004-2006 Jed Brown and Ed Bueler
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

#ifndef __iceROSSModel_hh
#define __iceROSSModel_hh

#include <petscda.h>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"

class IceROSSModel : public IceModel {
public:
    IceROSSModel(IceGrid &g, IceType &i);
    PetscInt               getflowlawNumber();
    void                   setflowlawNumber(PetscInt);
    virtual PetscErrorCode initFromOptions();
    PetscErrorCode         readROSSfiles();
    virtual PetscErrorCode run();
    PetscErrorCode         writeROSSfiles();
    PetscErrorCode         readRIGGSandCompare();

private:
    char            prefixROSS[PETSC_MAX_PATH_LEN];
    PetscInt        flowlawNumber, xsROSS, xmROSS;
    Vec             obsAzimuth, obsMagnitude, obsAccurate,
                    ubarBC, vbarBC;
    PetscScalar     gridLat[111], gridLon[147];
    PetscInt        kbcGridLoc[2][77], inletGridLoc[2][22];
    PetscScalar     vecErrAcc;
    
    PetscErrorCode  createROSSVecs();
    PetscErrorCode  destroyROSSVecs();
    PetscErrorCode  setBoundaryVels();
    PetscErrorCode  fillinTemps();
    PetscErrorCode  makeSurfaceFloating();
    PetscErrorCode  showObservedVels();
    PetscErrorCode  computeErrorsInAccurate();
    PetscErrorCode  runTune();
};

#endif /* __iceROSSModel_hh */
