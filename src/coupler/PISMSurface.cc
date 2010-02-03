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

#include "PISMSurface.hh"

///// PISMSurfaceModel base class:

void PISMSurfaceModel::attach_atmosphere_model(PISMAtmosphereModel *input) {
  if (atmosphere != NULL) {
    delete atmosphere;
  }
  atmosphere = input;
}

PetscErrorCode PISMSurfaceModel::init() {
  PetscErrorCode ierr;

  if (atmosphere == NULL)
    SETERRQ(1, "PISMSurfaceModel::init(): atmosphere == NULL");

  ierr = atmosphere->init(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMSurfaceModel::write_input_fields(PetscReal t_years, PetscReal dt_years,
						    string filename) {
  PetscErrorCode ierr;

  ierr = atmosphere->write_input_fields(t_years, dt_years, filename); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMSurfaceModel::write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
							 string filename) {
  PetscErrorCode ierr;

  ierr = atmosphere->write_diagnostic_fields(t_years, dt_years, filename); CHKERRQ(ierr);

  return 0;
}

///// Simple PISM surface model.

PetscErrorCode PSSimple::init() {
  PetscErrorCode ierr;

  if (atmosphere == NULL)
    SETERRQ(1, "PISMSurfaceModel::init(): atmosphere == NULL");

  ierr = atmosphere->init(); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
		    "* Initializing the simplest PISM surface model\n"
		    "  (precipitation == mass balance, 2m air temperature == ice surface temperature)...\n"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSSimple::ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
					       IceModelVec2 &result) {
  PetscErrorCode ierr;
  ierr = atmosphere->mean_precip(t_years, dt_years, result); CHKERRQ(ierr);

  string history = result.string_attr("history");
  history = "re-interpreted mean precipitation as surface mass balance\n" + history;
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSSimple::ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
						 IceModelVec2 &result) {
  PetscErrorCode ierr;
  ierr = atmosphere->mean_annual_temp(t_years, dt_years, result); CHKERRQ(ierr);

  string history = result.string_attr("history");
  history = "re-interpreted mean annual 2 m air temperature as instantaneous ice temperature at the ice surface\n" + history;
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

///// Constant-in-time surface model.

PetscErrorCode PSConstant::init() {
  PetscErrorCode ierr;
  LocalInterpCtx *lic = NULL;
  bool regrid = false;
  int start = -1;

  ierr = verbPrintf(2, grid.com, "  Initializing the constant-in-time surface process model...\n"); CHKERRQ(ierr);


  // allocate IceModelVecs for storing temperature and surface mass balance fields

  // create mean annual ice equivalent snow precipitation rate (before melt, and not including rain)
  ierr = acab.create(grid, "acab", false); CHKERRQ(ierr);
  ierr = acab.set_attrs("climate_state", 
			"constant-in-time ice-equivalent accumulation/ablation rate",
			"m s-1", 
			""); CHKERRQ(ierr); // no CF standard_name ??
  ierr = acab.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  acab.write_in_glaciological_units = true;

  ierr = artm.create(grid, "artm", false); CHKERRQ(ierr);
  ierr = artm.set_attrs("climate_state",
			"constant-in-time ice temperature (at the ice surface)",
			"K",
			""); CHKERRQ(ierr);
  
  // find PISM input file to read data from:

  ierr = find_pism_input(input_file, lic, regrid, start); CHKERRQ(ierr);

  // read snow precipitation rate and temperatures from file
  ierr = verbPrintf(2, grid.com, 
		    "    reading time-independent ice-equivalent accumulation/ablation rate 'acab'\n"
		    "    and time-independent ice temperature (at the ice surface) 'artm' from %s ... \n",
		    input_file.c_str()); CHKERRQ(ierr); 
  if (regrid) {
    ierr = acab.regrid(input_file.c_str(), *lic, true); CHKERRQ(ierr); // fails if not found!
    ierr = artm.regrid(input_file.c_str(), *lic, true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = acab.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
    ierr = artm.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }

  delete lic;

  t = grid.year;
  dt = 0;
	    
  return 0;
}

PetscErrorCode PSConstant::ice_surface_mass_flux(PetscReal /*t_years*/, PetscReal /*dt_years*/,
						 IceModelVec2 &result) {
  PetscErrorCode ierr;
  string history  = "read from " + input_file + "\n";

  ierr = acab.copy_to(result); CHKERRQ(ierr);
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSConstant::ice_surface_temperature(PetscReal /*t_years*/, PetscReal /*dt_years*/,
						   IceModelVec2 &result) {
  PetscErrorCode ierr;
  string history  = "read from " + input_file + "\n";

  ierr = artm.copy_to(result); CHKERRQ(ierr);
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

//! Does not ask the atmosphere model because it does not use one.
PetscErrorCode PSConstant::write_input_fields(PetscReal /*t_years*/, PetscReal /*dt_years*/,
						    string filename) {
  PetscErrorCode ierr;

  ierr = artm.write(filename.c_str(), NC_FLOAT); CHKERRQ(ierr);
  ierr = acab.write(filename.c_str(), NC_FLOAT); CHKERRQ(ierr);

  return 0;
}

//! Does not ask the atmosphere model because it does not use one.
PetscErrorCode PSConstant::write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
						   string filename) {
  PetscErrorCode ierr;

  ierr = write_input_fields(t_years, dt_years, filename); CHKERRQ(ierr);

  return 0;
}

///// PISM Surface model modifier.

void PSModifier::attach_input(PISMSurfaceModel *input) {
  if (input_model != NULL) {
    delete input_model;
  }
  input_model = input;
}

///// PISM surface model implementing a PDD scheme.

PSLocalMassBalance::PSLocalMassBalance(IceGrid &g, const NCConfigVariable &conf, PISMVars &vars)
  : PISMSurfaceModel(g, conf, vars) {
  mbscheme = NULL;
  use_fausto_pdd_parameters = false;
  lat = NULL;
}

PSLocalMassBalance::~PSLocalMassBalance() {
  delete mbscheme;
}

PetscErrorCode PSLocalMassBalance::init() {
  PetscErrorCode ierr;
  PetscTruth pdd_rand, pdd_rand_repeatable, fausto_params;

  ierr = PISMSurfaceModel::init(); CHKERRQ(ierr);

  ierr = check_option("-pdd_rand", pdd_rand); CHKERRQ(ierr);
  ierr = check_option("-pdd_rand_repeatable", pdd_rand_repeatable); CHKERRQ(ierr);
  ierr = check_option("-pdd_greenland", fausto_params); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com, "* Initializing the PDD-based surface process model...\n"); CHKERRQ(ierr);

  if (pdd_rand_repeatable) {
    ierr = verbPrintf(2, grid.com, "  Using a PDD implementation based on simulating a random process.\n"); CHKERRQ(ierr);
    mbscheme = new PDDrandMassBalance(config, true);
  } else if (pdd_rand) {
    ierr = verbPrintf(2, grid.com, "  Using a PDD implementation based on simulating a repeatable random process.\n\n"); CHKERRQ(ierr);
    mbscheme = new PDDrandMassBalance(config, false);
  } else {
    ierr = verbPrintf(2, grid.com, "  Using a PDD implementation based on an expectation integral.\n"); CHKERRQ(ierr);
    mbscheme = new PDDMassBalance(config);
  }

  if (fausto_params) {
    ierr = verbPrintf(2, grid.com, "  Setting PDD parameters using formulas (6) and (7) in [Faustoetal2009]...\n"); CHKERRQ(ierr);
    use_fausto_pdd_parameters = true;

    // allocate an IceModelVec2 to store mean July temperatures:
    ierr = temp_mj.create(grid, "temp_mj", false); CHKERRQ(ierr);
    ierr = temp_mj.set_attrs("internal",
			     "mean July temperature from the [\ref Faustoetal2009] parameterization",
			     "K",
			     ""); CHKERRQ(ierr);
    lat = dynamic_cast<IceModelVec2*>(variables.get("latitude"));
    if (!lat) SETERRQ(1, "ERROR: latitude is not available");
  }

  return 0;
}

PetscErrorCode PSLocalMassBalance::ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
							 IceModelVec2 &result) {
  PetscErrorCode ierr;
  PetscScalar **lat_degN;

  // to ensure that temperature time series are correct:
  ierr = atmosphere->update(t_years, dt_years); CHKERRQ(ierr);

  // This is a point-wise (local) computation, so we can use "result" to store
  // precipitation:
  ierr = atmosphere->mean_precip(t_years, dt_years, result); CHKERRQ(ierr);

  // set up air temperature time series
  PetscInt Nseries;
  ierr = mbscheme->getNForTemperatureSeries(t_years * secpera,
					    dt_years * secpera, Nseries); CHKERRQ(ierr);

  const PetscScalar tseries = (t_years - floor(t_years)) * secpera,
    dtseries = (dt_years * secpera) / ((PetscScalar) Nseries);

  // times for the air temperature time-series, in years:
  vector<PetscScalar> ts(Nseries), T(Nseries);
  for (PetscInt k = 0; k < Nseries; ++k)
    ts[k] = t_years + k * dt_years / Nseries;

  if (use_fausto_pdd_parameters) {
    // this is a nasty hack: time T is computed so that the call to
    // temp_snapshot produces a July temperature field. This is bad, because
    // this snapshot might have temperature offsets/anomalies applied to it
    // (and this time is the wrong one).
    const PetscReal
      sperd = 8.64e4, // exact number of seconds per day
      julydaysec = sperd * config.get("snow_temp_july_day"),
      T = floor(t_years) + julydaysec / secpera;

    ierr = atmosphere->temp_snapshot(T, 0.0, temp_mj); CHKERRQ(ierr);
    ierr = temp_mj.begin_access(); CHKERRQ(ierr);
    ierr = lat->get_array(lat_degN); CHKERRQ(ierr);
  }

  ierr = atmosphere->begin_pointwise_access(); CHKERRQ(ierr);
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = atmosphere->temp_time_series(i, j, Nseries, &ts[0], &T[0]); CHKERRQ(ierr);

      if (use_fausto_pdd_parameters) {
	// if we can (and if we are asked to), set mass balance parameters
	// according to formula (6) in [\ref Faustoetal2009]
	PDDMassBalance* pddscheme = dynamic_cast<PDDMassBalance*>(mbscheme);
	if (pddscheme != NULL) {
	  ierr = pddscheme->setDegreeDayFactorsFromSpecialInfo(lat_degN[i][j], temp_mj(i,j)); CHKERRQ(ierr);
	}
      }

      result(i,j) = mbscheme->getMassFluxFromTemperatureTimeSeries(tseries, dtseries, &T[0],
								   Nseries, result(i,j));
    }
  }
  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = atmosphere->end_pointwise_access(); CHKERRQ(ierr);

  if (use_fausto_pdd_parameters) {
  }


  return 0;
}

PetscErrorCode PSLocalMassBalance::ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
							   IceModelVec2 &result) {
  PetscErrorCode ierr;
  ierr = atmosphere->mean_annual_temp(t_years, dt_years, result); CHKERRQ(ierr);

  string history = result.string_attr("history");
  history = "re-interpreted mean annual near-surface air temperature as instantaneous ice temperature at the ice surface\n" + history;
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}
