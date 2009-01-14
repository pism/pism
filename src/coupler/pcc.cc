// Copyright (C) 2009 Ed Bueler and Ricarda Winkelmann
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


//! A draft derived class which shows how to use atmosphere and ocean instances of PISMClimateCoupler.
/*!
This class could be moved to a pair of files as usual (e.g. icePCCModel.hh
and icePCCModel.cc).

As is stands, IcePCCModel has two new members "atmosPCC" and "oceanPCC", instances of 
derived classes of PISMClimateCoupler.  Perhaps only one
of these is needed, and perhaps more. 

IcePCCModel could re-implement most IceModel procedures.  These
re-implementations would include climate data by calling the atmosPCC
and oceanPCC members.  Four procedures inherited from IceModel are listed 
as candidates for re-implementation.  They are just suggestions; re-implemention
of either just one or many IceModel procedures can be imagined.

For example, a re-implementation of IceModel::massContExplicitStep() could 
call atmosPCC.???() with the current ice flow model state, including the time, 
and get back a net surface mass balance.  Likewise, a re-implementation of
IceModel::temperatureStep() for ice shelves could call atmosPCC.???() to 
get a surface boundary temperature and also call oceanPCC.???() to get a
basal boundary temperature.

IcePCCModel could be changed to be a derived class \e of a derived
class of IceModel, such as IceEISModel or IceGRNModel.

Re-implementations could call the IceModel version and then do additional 
computation.  Or they could completely replace the existing computation.

The current empty implementation merely calls the IceModel version.  The
PISMClimateCoupler methods are not actually used.  Thus the executable 
\c pcc, if built, has identical function to \c pismr.
 */
class IcePCCModel : public IceModel {

public:
  IcePCCModel(IceGrid &g, IceType *i);

  PISMAtmosphereCoupler atmosPCC;
  PISMOceanCoupler oceanPCC;

  // could re-implement these procedures
  PetscErrorCode massContExplicitStep();
  PetscErrorCode temperatureStep(PetscScalar* vertSacrCount);
  PetscErrorCode additionalAtStartTimestep();
  PetscErrorCode additionalAtEndTimestep();
};


//! Does nothing except call to parent's constructor.  Could initialize stuff.
IcePCCModel::IcePCCModel(IceGrid &g, IceType *i) : IceModel(g, i) {
  atmosPCC.setGrid(&g);
  oceanPCC.setGrid(&g);
}


// see procedure in src/base/iMgeometry.cc
PetscErrorCode IcePCCModel::massContExplicitStep() {
  PetscErrorCode ierr;
  
  ierr = IceModel::massContExplicitStep(); CHKERRQ(ierr);
  return 0;
}

// see procedure in src/base/iMtemp.cc
PetscErrorCode IcePCCModel::temperatureStep(PetscScalar* vertSacrCount) {
  PetscErrorCode ierr;
  
  ierr = IceModel::temperatureStep(vertSacrCount); CHKERRQ(ierr);
  return 0;
}


// empty implementation found in src/base/iMutil.cc
PetscErrorCode IcePCCModel::additionalAtStartTimestep() {
  PetscErrorCode ierr;
  
  ierr = IceModel::additionalAtStartTimestep(); CHKERRQ(ierr);
  return 0;
}

// empty implementation found in src/base/iMutil.cc
PetscErrorCode IcePCCModel::additionalAtEndTimestep() {
  PetscErrorCode ierr;
  
  ierr = IceModel::additionalAtEndTimestep(); CHKERRQ(ierr);
  return 0;
}



// The executable 'pcc' is this procedure.
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
    ierr = verbPrintf(1,com, "PCC (demo of climate coupling mode)\n"); CHKERRQ(ierr);
    
    ierr = userChoosesIceType(com, ice); CHKERRQ(ierr); // allocates ice
    IcePCCModel m(g, ice);
    ierr = m.setExecName("pcc"); CHKERRQ(ierr);
    ierr = m.setFromOptions(); CHKERRQ(ierr);
    ierr = m.initFromOptions(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "running ...\n"); CHKERRQ(ierr);
    ierr = m.run(); CHKERRQ(ierr);

    ierr = verbPrintf(2,com, "... done with run\n"); CHKERRQ(ierr);

    // We provide a default base name if no -o option.
    ierr = m.writeFiles("unnamed_pcc.nc"); CHKERRQ(ierr);
    ierr = verbPrintf(2,com, "\n"); CHKERRQ(ierr);

    delete ice;
  }

  ierr = PetscFinalize(); CHKERRQ(ierr);
  return 0;
}
