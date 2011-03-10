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

#ifndef _PISMBEDTHERMALUNIT_H_
#define _PISMBEDTHERMALUNIT_H_

#include "PISMComponent.hh"
#include "iceModelVec.hh"
#include "PISMVars.hh"
#include "materials.hh"
#include "enthalpyConverter.hh"
#include "PISMDiagnostic.hh"


//! Given ice/bedrock interface temperature over a time-step, provides upward geothermal flux at the interface.
class PISMBedThermalUnit : public PISMComponent_TS {

public:
  PISMBedThermalUnit(IceGrid &g, EnthalpyConverter &e, const NCConfigVariable &conf);

  virtual ~PISMBedThermalUnit() { }

  virtual PetscErrorCode init(PISMVars &vars);

  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc, nc_type nctype);  
  virtual PetscErrorCode write_variables(set<string> vars, string filename);

  virtual PetscErrorCode max_timestep(PetscReal /*t_years*/, PetscReal &dt_years);

  virtual PetscErrorCode get_upward_geothermal_flux(PetscReal t_years, PetscReal dt_years,
                                                    IceModelVec2S &result);


protected:
  virtual PetscErrorCode allocate();

  IceModelVec3Bedrock temp;     //!< storage for bedrock thermal layer temperature;
                                //!    part of state; units K

  IceModelVec2S     ghf,        //!< storage for geothermal heat flux at base of
                                //!    bedrock thermal layer; part of state; units W m-2
                    ice_base_temp;

  // parameters of the heat equation:  T_t = D T_xx  where D = k / (rho c)
  PetscScalar       bed_rho, bed_c, bed_k, bed_D;

  EnthalpyConverter &EC; //!< needed to extract base temperature from ice enthalpy

  // pointers into IceModel space, generally:
  IceModelVec2Mask  *mask;     //!< is ice floating so BTU sees ocean temp?
  IceModelVec2S     *thk;      //!< needed to get ice base pressure
  IceModelVec3      *enthalpy; //!< needed to get ice base temperature
};

#endif /* _PISMBEDTHERMALUNIT_H_ */

