// Copyright (C) 2009 Ed Bueler
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

static char help[] = 
  "Draft driver for PISM simulations using PISMClimateCoupler.\n";

#include <petsc.h>
#include "../base/grid.hh"
#include "../base/materials.hh"
#include "../base/iceModel.hh"
#include "pccoupler.hh"


//! A draft derived class which shows how to use atmosphere and ocean instances of PISMClimateCoupler, to duplicate functionality of 'pismr'.
/*!
This class could be moved to a pair of files as usual (e.g. icePCCModel.hh
and icePCCModel.cc).

As is stands, IcePCCModel has two new members "atmosPCC" and "oceanPCC", 
instances of derived classes of PISMClimateCoupler.

IcePCCModel could re-implement most IceModel procedures.  These
re-implementations would include climate data by calling the atmosPCC
and oceanPCC members.  Four procedures inherited from IceModel are listed 
as candidates for re-implementation.  They are just suggestions; re-implemention
of either just one or many IceModel procedures can be imagined.

For example, a re-implementation of IceModel::massContExplicitStep() could 
call atmosPCC.updateSurfMassFluxAndProvide() with the current ice flow 
model state and time-step information and get back a net surface mass balance.  
It could call atmosPCC.updateSurfTempAndProvide() with the current ice flow 
model state and time-step information and get back a surface temperature.

IcePCCModel could be changed to be a derived class \e of a derived
class of IceModel, such as IceEISModel or IceGRNModel.

Re-implementations could call the IceModel version and then do additional 
computation.  Or they could completely replace the existing computation.

The intention of this draft implementation is that the executable 
\c pcc, if built, has identical function to \c pismr.
 */
class IcePCCModel : public IceModel {

public:
  IcePCCModel(IceGrid &g, IceType *i);

  //PISMAtmosphereCoupler atmosPCC;
  PISMConstAtmosCoupler atmosPCC;
  PISMOceanCoupler oceanPCC;

  // calls initialization for atmosPCC, oceanPCC
  virtual PetscErrorCode initFromOptions(PetscTruth doHook = PETSC_TRUE);

  // could re-implement these procedures
  PetscErrorCode write_model_state(const char filename[]);
  PetscErrorCode massContExplicitStep();
  PetscErrorCode temperatureStep(PetscScalar* vertSacrCount);
};


//! Does nothing except call to parent's constructor.
IcePCCModel::IcePCCModel(IceGrid &g, IceType *i) : IceModel(g, i) {
}


PetscErrorCode IcePCCModel::initFromOptions(PetscTruth doHook) {
  PetscErrorCode ierr;
  ierr = IceModel::initFromOptions(doHook); CHKERRQ(ierr);
  ierr = atmosPCC.initFromOptions(&grid); CHKERRQ(ierr);
  ierr = oceanPCC.initFromOptions(&grid); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IcePCCModel::write_model_state(const char filename[]) {
  PetscErrorCode ierr;
  ierr = IceModel::write_model_state(filename); CHKERRQ(ierr);
  ierr = atmosPCC.writeCouplingFieldsToFile(filename); CHKERRQ(ierr);
  return 0;
}

// see procedure in src/base/iMgeometry.cc
PetscErrorCode IcePCCModel::massContExplicitStep() {
  PetscErrorCode ierr;
  IceModelVec2* my_vsmf;

  ierr = atmosPCC.updateSurfMassFluxAndProvide(
             grid.year, dt / secpera, vMask, vh,
             my_vsmf); CHKERRQ(ierr);
  ierr = vAccum.copy_from(*my_vsmf); CHKERRQ(ierr);
  ierr = IceModel::massContExplicitStep(); CHKERRQ(ierr);
  return 0;
}


// see procedure in src/base/iMtemp.cc
PetscErrorCode IcePCCModel::temperatureStep(PetscScalar* vertSacrCount) {
  PetscErrorCode ierr;
  IceModelVec2* my_vst;

  ierr = atmosPCC.updateSurfTempAndProvide(
             grid.year, dt / secpera, vMask, vh,
             my_vst); CHKERRQ(ierr);
  ierr = vTs.copy_from(*my_vst); CHKERRQ(ierr);
  ierr = IceModel::temperatureStep(vertSacrCount); CHKERRQ(ierr);
  return 0;
}


// The executable 'pccr' is this procedure.
int main(int argc, char *argv[]) {
  PetscErrorCode  ierr;

  MPI_Comm    com;
  PetscMPIInt rank, size;

  ierr = PetscInitialize(&argc, &argv, PETSC_NULL, help); CHKERRQ(ierr);

  com = PETSC_COMM_WORLD;
  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(com, &size); CHKERRQ(ierr);

  /* This explicit scoping forces destructors to be called before PetscFinalize() */
  {
    IceGrid g(com, rank, size);
    IceType*   ice = PETSC_NULL;

    ierr = verbosityLevelFromOptions(); CHKERRQ(ierr);
    ierr = verbPrintf(1,com, "PCCR (demo of climate coupling mode)\n"); CHKERRQ(ierr);
    
    ierr = userChoosesIceType(com, ice); CHKERRQ(ierr); // allocates ice
    IcePCCModel m(g, ice);
    ierr = m.setExecName("pccr"); CHKERRQ(ierr);
    ierr = m.setFromOptions(); CHKERRQ(ierr);
    ierr = m.initFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "running ...\n"); CHKERRQ(ierr);
    ierr = m.run(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "... done with run\n"); CHKERRQ(ierr);

    // We provide a default base name if no -o option.
    ierr = m.writeFiles("unnamed_pccr.nc"); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "\n"); CHKERRQ(ierr);

    delete ice;
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
