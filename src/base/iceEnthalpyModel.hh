// Copyright (C) 2009-2010 Andreas Aschwanden, Ed Bueler and Constantine Khroulev
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

#ifndef __iceEnthalpyModel_hh
#define __iceEnthalpyModel_hh

#include <petsc.h>
#include "iceModelVec.hh"
#include "NCVariable.hh"
#include "iceModel.hh"


//! Temporary class for development of enthalpy-based polythermal PISM.
class IceEnthalpyModel : public IceModel {

public:
  IceEnthalpyModel(IceGrid &g, NCConfigVariable &conf, NCConfigVariable &conf_overrides);

protected:
  using IceModel::temperatureStep;
  virtual PetscErrorCode temperatureStep(PetscScalar* vertSacrCount, 
                                         PetscScalar* bulgeCount);

  virtual PetscErrorCode enthalpyAndDrainageStep(PetscScalar* vertSacrCount,
                                                 PetscScalar* liquifiedVol);

  virtual PetscErrorCode drainageToBaseModelEnth(
                PetscScalar omega_max, PetscScalar thickness,
                PetscScalar z, PetscScalar dz,
                PetscScalar &enthalpy, PetscScalar &Hmelt);
};

#endif  // __iceEnthalpyModel_hh

