// Copyright (C) 2008-2010 Ed Bueler, Constantine Khroulev, Ricarda Winkelmann,
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

#ifndef __PISMSurfaceModel_hh
#define __PISMSurfaceModel_hh

#include "../base/PISMComponent.hh"
#include "../base/iceModelVec.hh"
#include "PISMAtmosphere.hh"
#include "localMassBalance.hh"

class PISMSurfaceModel : public PISMComponent {
public:
  PISMSurfaceModel(IceGrid &g, const NCConfigVariable &conf)
    : PISMComponent(g, conf)
  { atmosphere = NULL; };

  virtual ~PISMSurfaceModel()
  { delete atmosphere; };

  virtual PetscErrorCode init(PISMVars &vars);
  virtual void attach_atmosphere_model(PISMAtmosphereModel *input);
  virtual PetscErrorCode ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
					       IceModelVec2 &result) = 0;
  virtual PetscErrorCode ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
						 IceModelVec2 &result) = 0;
  virtual PetscErrorCode write_model_state(PetscReal t_years, PetscReal dt_years,
					    string filename);
  virtual PetscErrorCode write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
						 string filename);
  virtual PetscErrorCode write_fields(set<string> vars, PetscReal t_years,
				      PetscReal dt_years, string filename);
protected:
  PISMAtmosphereModel *atmosphere;
};

//! A do-nothing (dummy) surface model. <b> Please avoid using it! </b>
class PSDummy : public PISMSurfaceModel {
public:
  PSDummy(IceGrid &g, const NCConfigVariable &conf)
    : PISMSurfaceModel(g, conf)
  {};

  virtual void attach_atmosphere_model(PISMAtmosphereModel *input)
  { delete input; }

  virtual PetscErrorCode init(PISMVars &) { return 0; };
  virtual PetscErrorCode ice_surface_mass_flux(PetscReal, PetscReal, IceModelVec2&)
  { return 0; }

  virtual PetscErrorCode ice_surface_temperature(PetscReal, PetscReal, IceModelVec2 &)
  { return 0; }
  virtual PetscErrorCode write_model_state(PetscReal, PetscReal, string)
  { return 0; }
  virtual PetscErrorCode write_diagnostic_fields(PetscReal, PetscReal, string)
  { return 0; }
  virtual PetscErrorCode write_fields(set<string>, PetscReal, PetscReal, string)
  { return 0; }
};

//! \brief A class implementing a primitive surface model.
/*! 
  This is an "invisible" model; it implements two modeling choices:
  \li accumulation obtained from an atmosphere model is interpreted as surface
  mass flux;

  \li mean-annual near-surface air temperature is interpreted as instantaneous
  temperature of the ice at the ice surface (i.e. the boundary condition of the
  conservation of energy scheme).
*/
class PSSimple : public PISMSurfaceModel {
public:
  PSSimple(IceGrid &g, const NCConfigVariable &conf)
    : PISMSurfaceModel(g, conf) {};
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
					       IceModelVec2 &result);
  virtual PetscErrorCode ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
						 IceModelVec2 &result);
};

//! \brief A class implementing a constant-in-time surface model. Reads data
//! from a PISM input file.
class PSConstant : public PISMSurfaceModel {
public:
  PSConstant(IceGrid &g, const NCConfigVariable &conf)
    : PISMSurfaceModel(g, conf)
  {};

  virtual PetscErrorCode init(PISMVars &vars);
  //! This surface model does not use an atmosphere model.
  virtual void attach_atmosphere_model(PISMAtmosphereModel *input)
  { delete input; }

  virtual PetscErrorCode ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
					       IceModelVec2 &result);
  virtual PetscErrorCode ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
						 IceModelVec2 &result);
  virtual PetscErrorCode write_model_state(PetscReal t_years, PetscReal dt_years,
					    string filename);
  virtual PetscErrorCode write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
						 string filename);
protected:
  string input_file;
  IceModelVec2 acab, artm;
};

class PSModifier : public PISMSurfaceModel {
public:
  PSModifier(IceGrid &g, const NCConfigVariable &conf)
    : PISMSurfaceModel(g, conf)
  { input_model = NULL; }

  virtual ~PSModifier()
  { delete input_model; }

  virtual void attach_input(PISMSurfaceModel *input);
protected:
  PISMSurfaceModel *input_model;
};

//! A class implementing a mechanism modifying surface mass balance to force
//! ice thickness to a given target at the end of the run.
class PSForceThickness : public PSModifier {
public:
  PSForceThickness(IceGrid &g, const NCConfigVariable &conf)
    : PSModifier(g, conf)
  {
    ice_thickness = NULL;
    alpha = config.get("force_to_thickness_alpha");
  }

  virtual ~PSForceThickness() {}
  PetscErrorCode init(PISMVars &vars);
  virtual void attach_atmosphere_model(PISMAtmosphereModel *input);
  virtual PetscErrorCode ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
					       IceModelVec2 &result);
  virtual PetscErrorCode ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
						 IceModelVec2 &result);
  virtual PetscErrorCode max_timestep(PetscReal t_years, PetscReal &dt_years);
  virtual PetscErrorCode write_model_state(PetscReal t_years, PetscReal dt_years,
					   string filename);
  virtual PetscErrorCode write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
						 string filename);
protected:
  string input_file;
  PetscReal alpha;
  IceModelVec2 *ice_thickness;	//!< current ice thickness produced by IceModel.
  IceModelVec2 target_thickness;
};

class PSLocalMassBalance : public PISMSurfaceModel {
public:
  PSLocalMassBalance(IceGrid &g, const NCConfigVariable &conf);
  virtual ~PSLocalMassBalance();
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
					       IceModelVec2 &result);
  virtual PetscErrorCode ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
						 IceModelVec2 &result);
protected:
  LocalMassBalance *mbscheme;	//!< mass balance scheme to use
  bool use_fausto_pdd_parameters;
  IceModelVec2 temp_mj,	//!< for the mean July temperature needed to set PDD parameters as in [\ref Faustoetal2009].
    *lat;		//!< latitude needed to set PDD parameters as in [\ref Faustoetal2009].
};

#endif	// __PISMSurfaceModel_hh
