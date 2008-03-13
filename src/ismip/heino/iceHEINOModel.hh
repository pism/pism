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

#ifndef __iceHEINOModel_hh
#define __iceHEINOModel_hh

#include <petscda.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"

class IceHEINOModel : public IceModel {
public:
    IceHEINOModel(IceGrid &g, IceType &i);
    virtual PetscErrorCode initFromOptions();
    void           setExperName(char);
    char           getExperName();
    void           setflowlawNumber(PetscInt);
    PetscInt       getflowlawNumber();
    PetscErrorCode simpFinalize();
    
private:
    char       expername;
    PetscTruth inFileSet;
    PetscInt   flowlawNumber;
 
    char         ismipRunName[3], heinodatprefix[20];
    PetscTruth   ismipNoDeliver, ismipForceDT;
    PetscInt     ismipHeinoRun;
    PetscViewer  ts[2], tss[3], tsp[7][3], // viewers (ASCII .dat files) for HEINO
                 pf; // viewer (ASCII .dat files) for 2D plan form
    PetscScalar  time1yr, time2yr, time3yr, time4yr;
    PetscTruth   time1Set, time2Set, time3Set, time4Set;
    bool         someTimeSet;
    PetscScalar  C_S;  // soft sediment sliding parameter
    
    PetscErrorCode setExperNameFromOptions();
    PetscErrorCode applyDefaultsForExperiment();
    PetscErrorCode initAccumTs();
    PetscErrorCode fillintemps();
    virtual PetscScalar basalVelocity(const PetscScalar x, const PetscScalar y, 
         const PetscScalar H, const PetscScalar T, const PetscScalar alpha,
         const PetscScalar mu);
    
    PetscErrorCode heinoCreateDat();
    PetscErrorCode heinoCloseDat();
    bool inSoftSediment(const PetscScalar x, const PetscScalar y);
    bool nearHeino(const PetscScalar x1, const PetscScalar y1,
                   const PetscScalar x2, const PetscScalar y2);
    virtual PetscErrorCode additionalAtStartTimestep();
    virtual PetscErrorCode additionalAtEndTimestep();
    PetscErrorCode planFormWrite(int nn);
};

#endif /* __iceHEINOModel_hh */
