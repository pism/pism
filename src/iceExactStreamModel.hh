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

#ifndef __iceExactSteamModel_hh
#define __iceExactSteamModel_hh

#include <petscda.h>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"

class IceExactSteamModel : public IceModel {
public:
    IceExactSteamModel(IceGrid &g, IceType &i);
    PetscInt               getflowlawNumber();
    void                   setflowlawNumber(PetscInt);
    virtual PetscErrorCode initFromOptions();
    virtual PetscErrorCode run();

private:
    PetscInt        flowlawNumber;
    static const PetscScalar   
                    m_schoof, L_schoof, aspect_schoof, h0_schoof,
                    B_schoof, p_schoof;
    
    PetscScalar     taucGet(PetscInt i, PetscInt j) const;
    PetscScalar     basalDragx(PetscScalar **u, PetscScalar **v,
                                   PetscInt i, PetscInt j) const;
    PetscScalar     basalDragy(PetscScalar **u, PetscScalar **v,
                                   PetscInt i, PetscInt j) const;
    PetscErrorCode  fillinTemps();
    PetscErrorCode  setBoundaryVels();
};

#endif /* __iceExactSteamModel_hh */
