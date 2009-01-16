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

#ifndef __pPDDcoupler_hh
#define __pPDDcoupler_hh

#include <petsc.h>
#include "../base/grid.hh"
#include "pccoupler.hh"

//! A derived class of PISMConstAtmosCoupler which provides a PDD to PISM.
/*!
The PDD here is the one already implemented in PISM.  That is, it is the one
from EISMINT-Greenland.  Thus it has various constants parameterizing the 
melt and refreeze processes.
 */
class PISMPDDCoupler : public PISMAtmosphereCoupler {

public:
  PISMPDDCoupler();
  ~PISMPDDCoupler();  // destroys PDD

  virtual PetscErrorCode initFromOptions(IceGrid* g);
  virtual PetscErrorCode writeCouplingFieldsToFile(const char *filename);

protected:
  IceModelVec2 vsurfaccum;
};


#endif

