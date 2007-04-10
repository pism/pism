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

#ifndef __iceExactStreamModel_hh
#define __iceExactStreamModel_hh

#include <petscda.h>
#include "grid.hh"
#include "materials.hh"
#include "iceModel.hh"

class IceExactStreamModel : public IceModel {
public:
    IceExactStreamModel(IceGrid &g, IceType &i);
    PetscInt               getflowlawNumber();
    void                   setflowlawNumber(PetscInt);
    virtual PetscErrorCode initFromOptions();
    virtual PetscErrorCode run();
    PetscErrorCode         reportErrors();

private:
    PetscInt        flowlawNumber;
    PetscTruth      exactOnly;
    static const PetscScalar   
                    m_schoof, L_schoof, aspect_schoof, H0_schoof,
                    B_schoof, p_schoof, DEFAULT_PLASTIC_REGULARIZE;
    
//     PetscScalar     basalDragx(PetscScalar **beta, PetscScalar **tauc,
//                                PetscScalar **u, PetscScalar **v,
//                                PetscInt i, PetscInt j) const;
//     PetscScalar     basalDragy(PetscScalar **beta, PetscScalar **tauc,
//                                PetscScalar **u, PetscScalar **v,
//                                PetscInt i, PetscInt j) const;
    PetscErrorCode  taucSet();
    PetscErrorCode  fillinTemps();
    PetscErrorCode  setInitStateAndBoundaryVels();
};

#endif /* __iceExactStreamModel_hh */
