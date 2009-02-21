// Copyright (C) 2004-2008 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __iceEISModel_hh
#define __iceEISModel_hh

#include <petscda.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"

//! This derived class does EISMINT II simplified geometry experiments.  
/*!
These experiments only involve the thermomechanically coupled shallow ice approximation.  See 
A. J. Payne and ten others, 200. <em> Results from the EISMINT model intercomparison:
the effects of thermomechanical coupling</em>.  J. Glaciol. 46(153), 227--238.
 */
class IceEISModel : public IceModel {
public:
    IceEISModel(IceGrid &g, IceType *i);
    virtual PetscErrorCode setFromOptions();
    using IceModel::initFromOptions;
    virtual PetscErrorCode initFromOptions(PetscTruth doHook = PETSC_TRUE);
    
protected:
    char        expername;
    bool        infileused;
    PetscScalar M_max, R_el, T_min, S_b, S_T;
 
    PetscErrorCode initAccumTs();
    PetscErrorCode fillintemps();
    virtual PetscScalar basalVelocity(const PetscScalar x, const PetscScalar y, 
         const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
         const PetscScalar mu);

    // for experiments I,J and K,L, respectively:
    PetscErrorCode generateTroughTopography();
    PetscErrorCode generateMoundTopography();
};

#endif /* __iceEISModel_hh */
