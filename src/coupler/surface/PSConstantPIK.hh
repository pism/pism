// Copyright (C) 2011, 2012, 2013 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#ifndef _PSCONSTANTPIK_H_
#define _PSCONSTANTPIK_H_

#include "PISMSurface.hh"
#include "iceModelVec.hh"
#include "PISMAtmosphere.hh"

//! \brief A class implementing a constant-in-time surface model for the surface mass balance.
//!
//! Reads data from a PISM input file.
//!
//! Ice surface temperature is parameterized as in PISM-PIK, using a latitude
//! and surface elevation-dependent formula.

class PSConstantPIK : public PISMSurfaceModel {
public:
  PSConstantPIK(IceGrid &g, const NCConfigVariable &conf);

  virtual PetscErrorCode init(PISMVars &vars);

  virtual void attach_atmosphere_model(PISMAtmosphereModel *input);

  virtual void get_diagnostics(map<string, PISMDiagnostic*> &dict,
                               map<string, PISMTSDiagnostic*> &ts_dict);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);
  virtual PetscErrorCode ice_surface_mass_flux(IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_temperature(IceModelVec2S &result);
  virtual PetscErrorCode define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, const PIO &nc);
  virtual void add_vars_to_output(string keyword, set<string> &result);
protected:
  string input_file;
  IceModelVec2S climatic_mass_balance, ice_surface_temp;
  IceModelVec2S *lat, *usurf;
private:
  PetscErrorCode allocate_PSConstantPIK();
};

#endif /* _PSCONSTANTPIK_H_ */
