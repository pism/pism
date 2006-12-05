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
    void                   setflowlawNumber(PetscInt);
    PetscInt               getflowlawNumber();
    virtual PetscErrorCode initFromOptions();
    virtual PetscErrorCode run();

private:
    PetscInt   flowlawNumber;
    Vec        obsAzimuth, obsMagnitude, obsAccurate;
    
    PetscErrorCode         readROSSfile();
    PetscErrorCode         createROSSVecs();
    PetscErrorCode         destroyROSSVecs();
};

#endif /* __iceROSSModel_hh */
