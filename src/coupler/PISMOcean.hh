// Copyright (C) 2008-2011 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
// Gudfinna Adalgeirsdottir and Andy Aschwanden
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

#ifndef __PISMOceanModel_hh
#define __PISMOceanModel_hh

#include "PISMComponent.hh"
#include "iceModelVec.hh"
#include "Timeseries.hh"

//! A very rudimentary PISM ocean model.
class PISMOceanModel : public PISMComponent_TS {
public:
  PISMOceanModel(IceGrid &g, const NCConfigVariable &conf)
    : PISMComponent_TS(g, conf)
  {
    sea_level = 0;
  }
  virtual ~PISMOceanModel() {};
  virtual PetscErrorCode sea_level_elevation(PetscReal t_years, PetscReal dt_years,
					     PetscReal &result) = 0;
  virtual PetscErrorCode shelf_base_temperature(PetscReal t_years, PetscReal dt_years,
						IceModelVec2S &result) = 0;
  virtual PetscErrorCode shelf_base_mass_flux(PetscReal t_years, PetscReal dt_years,
					      IceModelVec2S &result) = 0;
protected:
  PetscReal sea_level;
};

//! \brief A class implementing a constant (in terms of the ocean inputs) ocean
//! model. Uses configuration parameters for the sea level elevation and
//! sub-shelf heat flux.
class POConstant : public PISMOceanModel {
public:
  POConstant(IceGrid &g, const NCConfigVariable &conf);
  virtual ~POConstant() {}
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode sea_level_elevation(PetscReal t_years, PetscReal dt_years,
					     PetscReal &result);
  virtual PetscErrorCode shelf_base_temperature(PetscReal t_years, PetscReal dt_years,
						IceModelVec2S &result);
  virtual PetscErrorCode shelf_base_mass_flux(PetscReal t_years, PetscReal dt_years,
					      IceModelVec2S &result);

  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc,
                                          nc_type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
protected:
  IceModelVec2S *ice_thickness;	// is not owned by this class
  NCSpatialVariable shelfbmassflux, shelfbtemp;
};

//! A class defining the interface of a PISM ocean model modifier.
class POModifier : public PISMOceanModel {
public:
  POModifier(IceGrid &g, const NCConfigVariable &conf)
    : PISMOceanModel(g, conf)
  { input_model = NULL; }

  virtual ~POModifier()
  { delete input_model; }

  virtual void attach_input(PISMOceanModel *input);
  virtual void add_vars_to_output(string key, set<string> &result) {
    if (input_model != NULL)
      input_model->add_vars_to_output(key, result);
  }

  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc,
                                          nc_type nctype) {
    PetscErrorCode ierr;
    if (input_model != NULL) {
      ierr = input_model->define_variables(vars, nc, nctype); CHKERRQ(ierr);
    }
    return 0;
  }

  virtual PetscErrorCode write_variables(set<string> vars, string filename) {
    PetscErrorCode ierr;
    if (input_model != NULL) {
      ierr = input_model->write_variables(vars, filename); CHKERRQ(ierr);
    }
    return 0;
  }

  virtual void get_diagnostics(map<string, PISMDiagnostic*> &dict) {
    if (input_model)
      input_model->get_diagnostics(dict);
  }

protected:
  PISMOceanModel *input_model;
};

//! A class implementing sea level forcing.
class POForcing : public POModifier {
public:
  POForcing(IceGrid &g, const NCConfigVariable &conf)
    : POModifier(g, conf)
  {
    dSLforcing = NULL;
    delta_sea_level = NULL;
  }

  virtual ~POForcing()
  {
    delete dSLforcing;
    delete delta_sea_level;
  }

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode sea_level_elevation(PetscReal t_years, PetscReal dt_years,
					     PetscReal &result);
  virtual PetscErrorCode shelf_base_temperature(PetscReal t_years, PetscReal dt_years,
						IceModelVec2S &result);
  virtual PetscErrorCode shelf_base_mass_flux(PetscReal t_years, PetscReal dt_years,
					      IceModelVec2S &result);
protected:
  string forcing_file;
  Timeseries *dSLforcing;	//!< sea level forcing time-series
  DiagnosticTimeseries *delta_sea_level;
};




//! \brief A class implementing an ocean model.
//! Parameterization of sub-shelf melting with respect to sub-shelf heat flux like in Beckmann_Goosse 2003
class POConstantPIK : public PISMOceanModel {
public:
  POConstantPIK(IceGrid &g, const NCConfigVariable &conf);
  virtual ~POConstantPIK() {}
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode sea_level_elevation(PetscReal t_years, PetscReal dt_years,
					     PetscReal &result);
  virtual PetscErrorCode shelf_base_temperature(PetscReal t_years, PetscReal dt_years,
						IceModelVec2S &result);
  virtual PetscErrorCode shelf_base_mass_flux(PetscReal t_years, PetscReal dt_years,
					      IceModelVec2S &result);

  virtual void add_vars_to_output(string keyword, set<string> &result);
  virtual PetscErrorCode define_variables(set<string> vars, const NCTool &nc,
                                          nc_type nctype);
  virtual PetscErrorCode write_variables(set<string> vars, string filename);
protected:
  IceModelVec2S *ice_thickness;	// is not owned by this class
  NCSpatialVariable shelfbmassflux, shelfbtemp;
};



#endif	// __PISMOceanModel_hh
