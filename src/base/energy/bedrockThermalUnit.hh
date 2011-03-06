// Copyright (C) 2011 Ed Bueler
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

#ifndef _BEDROCKTHERMALUNIT_H_
#define _BEDROCKTHERMALUNIT_H_

#include "PISMComponent.hh"
#include "iceModelVec.hh"
#include "PISMVars.hh"
#include "materials.hh"
#include "enthalpyConverter.hh"
#include "PISMDiagnostic.hh"


//! Given ice/bedrock interface temperature over a time-step, provides upward geothermal flux at the interface.
class BedrockThermalUnit : public PISMComponent_TS {

public:
  BedrockThermalUnit(IceGrid &g, EnthalpyConverter &e, const NCConfigVariable &conf);

  virtual ~BedrockThermalUnit() { }

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode max_timestep(PetscReal /*t_years*/, PetscReal &dt_years);

  virtual PetscErrorCode update(PetscReal /*t_years*/, PetscReal /*dt_years*/);

  virtual PetscErrorCode write_model_state(PetscReal /*t_years*/, PetscReal /*dt_years*/,
					   string filename);

protected:
  virtual PetscErrorCode allocate();

  EnthalpyConverter &EC;

  IceModelVec2S     basal_temp;

  PetscScalar       rho, c, k, D;

  // pointers into IceModel space, generally:
  IceModelVec2Mask  *mask;
  IceModelVec3      *enthalpy;
  // FIXME: need thickness too, to compute basal temp from enthalpy and pressure
};

#endif /* _BEDROCKTHERMALUNIT_H_ */

