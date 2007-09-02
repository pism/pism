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

#ifndef __iceExactSSAModel_hh
#define __iceExactSSAModel_hh

#include <petscda.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"

class IceExactSSAModel : public IceModel {
public:
    IceExactSSAModel(IceGrid &g, IceType &i, const char mytest);
    virtual PetscErrorCode initFromOptions();
    virtual PetscErrorCode diagnosticRun();
    PetscErrorCode         reportErrors();

protected:
    char       test;       // only 'I', 'J' supported
    PetscTruth exactOnly;
    Vec*       vNuForJ;
          
    PetscErrorCode  fillFromExactSolution();
    PetscErrorCode  taucSetI();
    PetscErrorCode  setInitStateAndBoundaryVelsI();
    PetscErrorCode  setInitStateJ();

private:
    // constants for I
    static const PetscScalar   
               m_schoof, L_schoof, aspect_schoof, H0_schoof,
               B_schoof, p_schoof, DEFAULT_PLASTIC_REGULARIZE;
    // constants for J
    static const PetscScalar LforJ;
};

#endif /* __iceExactSSAModel_hh */
