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
#include "../base/PISMIO.hh"

///// PISMSurfaceModel base class:

void PISMSurfaceModel::attach_atmosphere_model(PISMAtmosphereModel *input) {
  if (atmosphere != NULL) {
    delete atmosphere;
  }
  atmosphere = input;
}

PetscErrorCode PISMSurfaceModel::init(PISMVars &vars) {
  PetscErrorCode ierr;

  if (atmosphere == NULL)
    SETERRQ(1, "PISMSurfaceModel::init(PISMVars &vars): atmosphere == NULL");

  ierr = atmosphere->init(vars); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PISMSurfaceModel::write_model_state(PetscReal t_years, PetscReal dt_years,
						    string filename) {
  PetscErrorCode ierr;

  if (atmosphere != NULL) {
    ierr = atmosphere->write_model_state(t_years, dt_years, filename); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PISMSurfaceModel::write_fields(set<string> vars, PetscReal t_years,
					      PetscReal dt_years, string filename) {
  PetscErrorCode ierr;

  if (atmosphere != NULL) {
    ierr = atmosphere->write_fields(vars, t_years, dt_years, filename); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PISMSurfaceModel::write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
							 string filename) {
  PetscErrorCode ierr;

  if (atmosphere != NULL) {
    ierr = atmosphere->write_diagnostic_fields(t_years, dt_years, filename); CHKERRQ(ierr);
  }

  return 0;
}

///// Simple PISM surface model.

PetscErrorCode PSSimple::init(PISMVars &vars) {
  PetscErrorCode ierr;

  if (atmosphere == NULL)
    SETERRQ(1, "PISMSurfaceModel::init(PISMVars &vars): atmosphere == NULL");

  ierr = atmosphere->init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
     "* Initializing the simplest PISM surface (snow) processes model PSSimple.\n"
     "  It passes atmospheric state directly to upper ice fluid surface:\n"
     "    surface mass balance          := precipitation,\n"
     "    ice upper surface temperature := 2m air temperature.\n"
     "  Any choice of atmosphere coupler (option '-atmosphere') is ignored.\n");
     CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PSSimple::ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
					       IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = atmosphere->mean_precip(t_years, dt_years, result); CHKERRQ(ierr);

  string history = result.string_attr("history");
  history = "re-interpreted precipitation as surface mass balance (PSSimple)\n" + history;
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSSimple::ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
						 IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = atmosphere->mean_annual_temp(t_years, dt_years, result); CHKERRQ(ierr);

  string history = result.string_attr("history");
  history = "re-interpreted mean annual 2 m air temperature as instantaneous ice temperature at the ice surface (PSSimple)\n" + history;
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

///// Constant-in-time surface model.

PetscErrorCode PSConstant::init(PISMVars &/*vars*/) {
  PetscErrorCode ierr;
  LocalInterpCtx *lic = NULL;
  bool regrid = false;
  int start = -1;

  ierr = verbPrintf(2, grid.com, 
     "* Initializing the constant-in-time surface processes model PSConstant.\n"
     "  It reads surface mass balance and ice upper-surface temperature\n"
     "  directly from the file and holds them constant.\n"); CHKERRQ(ierr);

  // allocate IceModelVecs for storing temperature and surface mass balance fields

  // create mean annual ice equivalent snow precipitation rate (before melt, and not including rain)
  ierr = acab.create(grid, "acab", false); CHKERRQ(ierr);
  ierr = acab.set_attrs("climate_state", 
			"constant-in-time ice-equivalent surface mass balance (accumulation/ablation) rate",
			"m s-1", 
			"land_ice_surface_specific_mass_balance"); // CF standard_name
			CHKERRQ(ierr);
  ierr = acab.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  acab.write_in_glaciological_units = true;

  ierr = artm.create(grid, "artm", false); CHKERRQ(ierr);
  ierr = artm.set_attrs("climate_state",
			"constant-in-time ice temperature at the ice surface",
			"K",
			""); CHKERRQ(ierr);
  
  // find PISM input file to read data from:
  ierr = find_pism_input(input_file, lic, regrid, start); CHKERRQ(ierr);

  // read snow precipitation rate and temperatures from file
  ierr = verbPrintf(2, grid.com, 
    "    reading ice-equivalent surface mass balance (accumulation/ablation) rate 'acab'\n"
    "    and ice surface temperature  'artm' from %s ... \n",
    input_file.c_str()); CHKERRQ(ierr); 
  if (regrid) {
    ierr = acab.regrid(input_file.c_str(), *lic, true); CHKERRQ(ierr); // fails if not found!
    ierr = artm.regrid(input_file.c_str(), *lic, true); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = acab.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
    ierr = artm.read(input_file.c_str(), start); CHKERRQ(ierr); // fails if not found!
  }

  delete lic;
  return 0;
}

PetscErrorCode PSConstant::ice_surface_mass_flux(PetscReal /*t_years*/, PetscReal /*dt_years*/,
						 IceModelVec2S &result) {
  PetscErrorCode ierr;
  string history  = "read from " + input_file + "\n";

  ierr = acab.copy_to(result); CHKERRQ(ierr);
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSConstant::ice_surface_temperature(PetscReal /*t_years*/, PetscReal /*dt_years*/,
						   IceModelVec2S &result) {
  PetscErrorCode ierr;
  string history  = "read from " + input_file + "\n";

  ierr = artm.copy_to(result); CHKERRQ(ierr);
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

//! Does not ask the atmosphere model because it does not use one.
PetscErrorCode PSConstant::write_model_state(PetscReal /*t_years*/, PetscReal /*dt_years*/,
						    string filename) {
  PetscErrorCode ierr;

  ierr = artm.write(filename.c_str()); CHKERRQ(ierr);
  ierr = acab.write(filename.c_str()); CHKERRQ(ierr);

  return 0;
}

//! Does not ask the atmosphere model because it does not use one.
PetscErrorCode PSConstant::write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
						   string filename) {
  PetscErrorCode ierr;

  ierr = write_model_state(t_years, dt_years, filename); CHKERRQ(ierr);

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

PSTemperatureIndex::PSTemperatureIndex(IceGrid &g, const NCConfigVariable &conf)
  : PISMSurfaceModel(g, conf) {
  mbscheme = NULL;
  faustogreve = NULL;
  base_ddf.snow = config.get("pdd_factor_snow");
  base_ddf.ice  = config.get("pdd_factor_ice");
  base_ddf.refreezeFrac = config.get("pdd_refreeze");
  base_pddStdDev = config.get("pdd_std_dev");
  base_pddThresholdTemp = config.get("pdd_positive_threshold_temp");
}

PSTemperatureIndex::~PSTemperatureIndex() {
  delete mbscheme;
  delete faustogreve;
}

PetscErrorCode PSTemperatureIndex::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool           pdd_rand, pdd_rand_repeatable, fausto_params, pSet;
  PetscScalar    factor;

  ierr = PISMSurfaceModel::init(vars); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", 
                           "Temperature-index (PDD) scheme for surface (snow) processes", "");
                           CHKERRQ(ierr);
  {
    ierr = PISMOptionsIsSet("-pdd_rand",
                            "Use a PDD implementation based on simulating a random process",
			    pdd_rand); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-pdd_rand_repeatable",
                            "Use a PDD implementation based on simulating a repeatable random process",
			    pdd_rand_repeatable); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-pdd_fausto",
                            "Set PDD parameters using formulas (6) and (7) in [Faustoetal2009]",
			    fausto_params); CHKERRQ(ierr);

    ierr = PISMOptionsReal("-pdd_factor_snow", "PDD snow factor",
                           base_ddf.snow, pSet); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-pdd_factor_ice", "PDD ice factor",
                           base_ddf.ice, pSet); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-pdd_refreeze", "PDD refreeze fraction",
                           base_ddf.refreezeFrac, pSet); CHKERRQ(ierr);

    ierr = PISMOptionsReal("-pdd_std_dev", "PDD standard deviation",
                           base_pddStdDev, pSet); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-pdd_positive_threshold_temp", 
                           "PDD uses this temp in K to determine 'positive' temperatures",
                           base_pddThresholdTemp, pSet); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);


  ierr = verbPrintf(2, grid.com,
    "* Initializing the default temperature-index, PDD-based surface (snow) processes\n"
    "  scheme.  Precipitation and 2m air temperature provided by atmosphere above\n"
    "  (PISMAtmosphere) are inputs.  Surface mass balance and ice fluid upper surface\n"
    "  temperature are outputs.  See PISM User's Manual for control of PDD parameters\n"
    "  (degree-day factors).\n");
    CHKERRQ(ierr);

  if (pdd_rand_repeatable) {
    ierr = verbPrintf(2, grid.com,
      "  Using a PDD implementation based on simulating a random process.\n");
      CHKERRQ(ierr);
    mbscheme = new PDDrandMassBalance(config, true);
  } else if (pdd_rand) {
    ierr = verbPrintf(2, grid.com,
      "  Using a PDD implementation based on simulating a repeatable random process.\n\n");
      CHKERRQ(ierr);
    mbscheme = new PDDrandMassBalance(config, false);
  } else {
    ierr = verbPrintf(2, grid.com,
      "  Using a PDD implementation based on an expectation integral.\n");
      CHKERRQ(ierr);
    mbscheme = new PDDMassBalance(config);
  }

  if (config.get_flag("pdd_limit_timestep")) {
    ierr = verbPrintf(2, grid.com, "  NOTE: Limiting time-steps to 1 year.\n"); CHKERRQ(ierr);
  }

  ierr = acab.create(grid, "acab", false); CHKERRQ(ierr);
  ierr = acab.set_attrs("diagnostic",
			"instantaneous ice-equivalent surface mass balance (accumulation/ablation) rate",
			"m s-1",  // m *ice-equivalent* per second
			"land_ice_surface_specific_mass_balance");  // CF standard_name
                        CHKERRQ(ierr);
  ierr = acab.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  acab.write_in_glaciological_units = true;
  ierr = acab.set_attr("comment", "positive values correspond to ice gain"); CHKERRQ(ierr); 

  // diagnostic fields:

  ierr = accumulation_rate.create(grid, "saccum", false); CHKERRQ(ierr);
  ierr = accumulation_rate.set_attrs("diagnostic",
                                     "instantaneous ice-equivalent surface accumulation rate (precip minus rain)",
                                     "m s-1",
                                     ""); CHKERRQ(ierr);
  ierr = accumulation_rate.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  accumulation_rate.write_in_glaciological_units = true;

  ierr = melt_rate.create(grid, "smelt", false); CHKERRQ(ierr);
  ierr = melt_rate.set_attrs("diagnostic",
                             "instantaneous ice-equivalent surface melt rate",
                             "m s-1",
                             ""); CHKERRQ(ierr);
  ierr = melt_rate.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  melt_rate.write_in_glaciological_units = true;

  ierr = runoff_rate.create(grid, "srunoff", false); CHKERRQ(ierr);
  ierr = runoff_rate.set_attrs("diagnostic",
                               "instantaneous ice-equivalent surface meltwater runoff rate",
                               "m s-1",
                               ""); CHKERRQ(ierr);
  ierr = runoff_rate.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  runoff_rate.write_in_glaciological_units = true;

  if (fausto_params) {
    ierr = verbPrintf(2, grid.com, 
       "  Setting PDD parameters using formulas (6) and (7) in [Faustoetal2009]...\n");
       CHKERRQ(ierr);

    lat = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
    if (!lat) SETERRQ(1, "ERROR: 'latitude' is not available or is of wrong type in vars dictionary");
    lon = dynamic_cast<IceModelVec2S*>(vars.get("longitude"));
    if (!lon) SETERRQ(1, "ERROR: 'longitude' is not available or is of wrong type in vars dictionary");
    usurf = dynamic_cast<IceModelVec2S*>(vars.get("usurf"));
    if (!usurf) SETERRQ(1, "ERROR: 'usurf' is not available or is of wrong type in vars dictionary");
   
    faustogreve = new FaustoGrevePDDObject(grid, config);
  } else {
    // generally, this is the case in which degree day factors do not depend
    //   on location; we use base_ddf
    lat = NULL;
    lon = NULL;
    usurf = NULL;
  }

  return 0;
}

PetscErrorCode PSTemperatureIndex::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;
  PetscScalar **lat_degN;

  if ((fabs(t_years - t) < 1e-12) &&
      (fabs(dt_years - dt) < 1e-12))
    return 0;

  t  = t_years;
  dt = dt_years;

  // to ensure that temperature time series are correct:
  ierr = atmosphere->update(t_years, dt_years); CHKERRQ(ierr);

  // This is a point-wise (local) computation, so we can use "acab" to store
  // precipitation:
  ierr = atmosphere->mean_precip(t_years, dt_years, acab); CHKERRQ(ierr);

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

  if (faustogreve != NULL) {
    if (lat == NULL) { SETERRQ(1,"faustogreve object is allocated BUT lat==NULL"); }
    if (lon == NULL) { SETERRQ(2,"faustogreve object is allocated BUT lon==NULL"); }
    if (usurf == NULL) { SETERRQ(3,"faustogreve object is allocated BUT usurf==NULL"); }
    ierr = lat->begin_access(); CHKERRQ(ierr)
    ierr = lon->begin_access(); CHKERRQ(ierr)
    ierr = usurf->begin_access(); CHKERRQ(ierr)
    ierr = faustogreve->update_temp_mj(usurf, lat, lon); CHKERRQ(ierr);
  }

  DegreeDayFactors  ddf = base_ddf;

  ierr = atmosphere->begin_pointwise_access(); CHKERRQ(ierr);
  ierr = acab.begin_access(); CHKERRQ(ierr);

  ierr = accumulation_rate.begin_access(); CHKERRQ(ierr);
  ierr = melt_rate.begin_access(); CHKERRQ(ierr);
  ierr = runoff_rate.begin_access(); CHKERRQ(ierr);
  
  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = atmosphere->temp_time_series(i, j, Nseries, &ts[0], &T[0]); CHKERRQ(ierr);

      if (faustogreve != NULL) {
	// we have been asked to set mass balance parameters according to
	//   formula (6) in [\ref Faustoetal2009]; they overwrite ddf set above
	ierr = faustogreve->setDegreeDayFactors(i,j,(*usurf)(i,j),(*lat)(i,j),(*lon)(i,j),ddf);
	           CHKERRQ(ierr);
      }

      PetscScalar pddsum = mbscheme->getPDDSumFromTemperatureTimeSeries(
                                  base_pddStdDev, base_pddThresholdTemp,
                                  tseries, dtseries, &T[0], Nseries);

      PetscScalar snow_amount = mbscheme->getSnowFromPrecipAndTemperatureTimeSeries(
                                  acab(i,j), // precipitation rate (input)
                                  tseries, dtseries, &T[0], Nseries);

      ierr = mbscheme->getMassFluxesFromPDDs(ddf,
                                             dt_years * secpera, pddsum, snow_amount,
                                             accumulation_rate(i,j), // output
                                             melt_rate(i,j), // output
                                             runoff_rate(i,j), // output
                                             acab(i,j)); // acab = smb (output)
                                             CHKERRQ(ierr); 
    }
  }

  ierr = accumulation_rate.end_access(); CHKERRQ(ierr);
  ierr = melt_rate.end_access(); CHKERRQ(ierr);
  ierr = runoff_rate.end_access(); CHKERRQ(ierr);

  ierr = acab.end_access(); CHKERRQ(ierr);
  ierr = atmosphere->end_pointwise_access(); CHKERRQ(ierr);

  if (faustogreve != NULL) {
    ierr = lat->end_access(); CHKERRQ(ierr)
    ierr = lon->end_access(); CHKERRQ(ierr)
    ierr = usurf->end_access(); CHKERRQ(ierr)
  }

  return 0;
}


PetscErrorCode PSTemperatureIndex::ice_surface_mass_flux(PetscReal t_years, PetscReal dt_years,
							 IceModelVec2S &result) {
  PetscErrorCode ierr;

  // This flag is set in pclimate to allow testing the model. It does not
  // affect the normal operation of PISM.
  if (config.get_flag("pdd_limit_timestep")) {
    dt_years = PetscMin(dt_years, 1.0);
  }

  ierr = update(t_years, dt_years); CHKERRQ(ierr);

  ierr = acab.copy_to(result); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PSTemperatureIndex::ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
							   IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = atmosphere->mean_annual_temp(t_years, dt_years, result); CHKERRQ(ierr);

  string history = result.string_attr("history");
  history = "re-interpreted mean annual near-surface air temperature as instantaneous ice temperature at the ice surface\n" + history;
  ierr = result.set_attr("history", history); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSTemperatureIndex::write_fields(set<string> vars, PetscReal t_years,
					      PetscReal dt_years, string filename) {
  PetscErrorCode ierr;

  ierr = PISMSurfaceModel::write_fields(vars, t_years, dt_years, filename); CHKERRQ(ierr);

  if (vars.find("saccum") != vars.end()) {
    ierr = accumulation_rate.write(filename.c_str()); CHKERRQ(ierr);
  }  

  if (vars.find("smelt") != vars.end()) {
    ierr = melt_rate.write(filename.c_str()); CHKERRQ(ierr);
  }  

  if (vars.find("srunoff") != vars.end()) {
    ierr = runoff_rate.write(filename.c_str()); CHKERRQ(ierr); 
  }  

  return 0;
}

PetscErrorCode PSTemperatureIndex::write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
                                                           string filename) {
  PetscErrorCode ierr;

  ierr = PISMSurfaceModel::write_diagnostic_fields(t_years, dt_years, filename); CHKERRQ(ierr);

  ierr = accumulation_rate.write(filename.c_str()); CHKERRQ(ierr);
  ierr = melt_rate.write(filename.c_str()); CHKERRQ(ierr);
  ierr = runoff_rate.write(filename.c_str()); CHKERRQ(ierr); 

  return 0;
}


///// "Force-to-thickness" mechanism

void PSForceThickness::attach_atmosphere_model(PISMAtmosphereModel *input) {
  input_model->attach_atmosphere_model(input);
}

PetscErrorCode PSForceThickness::init(PISMVars &vars) {
  PetscErrorCode ierr;
  char fttfile[PETSC_MAX_PATH_LEN] = "";
  PetscTruth opt_set;
  PetscScalar fttalpha;
  PetscTruth  fttalphaSet;

  ierr = input_model->init(vars); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Surface model forcing", ""); CHKERRQ(ierr);

  ierr = PetscOptionsString("-force_to_thk",
			    "Specifies the target thickness file",
			    "", "",
			    fttfile, PETSC_MAX_PATH_LEN, &opt_set); CHKERRQ(ierr);

  if (!opt_set) {
    ierr = PetscPrintf(grid.com,
      "ERROR: surface model forcing requires the -force_to_thk option.\n"); CHKERRQ(ierr);
    PetscEnd();
  }
    
  ierr = PetscOptionsReal("-force_to_thk_alpha",
			  "Specifies the force-to-thickness alpha value in per-year units",
			  "", alpha * secpera,
			  &fttalpha, &fttalphaSet); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
		    "* Initializing force-to-thickness mass-balance modifier...\n"); CHKERRQ(ierr);

  ice_thickness = dynamic_cast<IceModelVec2S*>(vars.get("land_ice_thickness"));
  if (!ice_thickness) SETERRQ(1, "ERROR: land_ice_thickness is not available");

  ierr = target_thickness.create(grid, "thk", false); CHKERRQ(ierr); // name to read by
  ierr = target_thickness.set_attrs(
     "climate_state", 
     "target ice thickness (to be reached at the end of the run",
     "m", 
     "land_ice_thickness"); CHKERRQ(ierr); // standard_name to read by

  ierr = ftt_mask.create(grid, "ftt_mask", false); CHKERRQ(ierr);
  ierr = ftt_mask.set_attrs(
     "climate_state",
     "mask specifying where to apply the force-to-thickness mechanism",
     "", ""); CHKERRQ(ierr); // no units and no standard name
  ierr = ftt_mask.set(1.0); CHKERRQ(ierr); // default: applied in whole domain

  ierr = ftt_modified_acab.create(grid, "ftt_modified_acab", false); CHKERRQ(ierr);
  ierr = ftt_modified_acab.set_attrs(
     "diagnostic",
     "modified ice-equivalent surface mass balance (accumulation/ablation) rate;"
       " result from force-to-thickness PSModifier",
     "m s-1", 
     ""); CHKERRQ(ierr); // no standard name
  ierr = ftt_modified_acab.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  ftt_modified_acab.write_in_glaciological_units = true;
  ierr = ftt_modified_acab.set_attr("comment", "positive values correspond to ice gain"); CHKERRQ(ierr); 
  ierr = ftt_modified_acab.set(0.0); CHKERRQ(ierr); // no useful default

  input_file = fttfile;

  // determine exponential rate alpha from user option or from factor; option
  // is given in a^{-1}
  if (fttalphaSet == PETSC_TRUE) {
    ierr = verbPrintf(3, grid.com, "    option -force_to_thk_alpha seen\n");
       CHKERRQ(ierr);
    alpha = fttalpha / secpera;
  }
    
  ierr = verbPrintf(2, grid.com,
		    "    alpha = %.6f a-1 for %.3f a run, for -force_to_thk mechanism\n",
		    alpha * secpera, grid.end_year - grid.start_year); CHKERRQ(ierr);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  // fttfile now contains name of -force_to_thk file; now check
  // it is really there; if so, read the dimensions of computational grid so
  // that we can set up a LocalInterpCtx for actual reading of target thickness
  PISMIO nc(&grid);
  grid_info gi;
  bool mask_exists = false;
  ierr = nc.open_for_reading(fttfile); CHKERRQ(ierr);
  ierr = nc.get_grid_info_2d(gi); CHKERRQ(ierr);
  ierr = nc.find_variable("ftt_mask", NULL, mask_exists); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  LocalInterpCtx* lic;
  lic = new LocalInterpCtx(gi, NULL, NULL, grid); // 2D only
  ierr = verbPrintf(2, grid.com, 
		    "    reading target thickness 'thk' from %s ...\n", fttfile); CHKERRQ(ierr); 
  ierr = target_thickness.regrid(fttfile, *lic, true); CHKERRQ(ierr);

  if (mask_exists) {
    ierr = verbPrintf(2, grid.com, 
                      "    reading force-to-thickness mask 'ftt_mask' from %s ...\n", fttfile); CHKERRQ(ierr); 
    ierr = ftt_mask.regrid(fttfile, *lic, true); CHKERRQ(ierr);
    write_ftt_mask = true;
  }

  delete lic;

  // reset name to avoid confusion; attributes again because lost by set_name()?
  ierr = target_thickness.set_name("target_thickness"); CHKERRQ(ierr);
  ierr = target_thickness.set_attrs(
    "",  // pism_intent unknown
    "target thickness for force-to-thickness-spinup mechanism (hit this at end of run)",
    "m",
    "");  // no CF standard_name, to put it mildly
  CHKERRQ(ierr);

  return 0;
}

/*!
If \c -force_to_thk \c foo.nc is in use then vthktarget will have a target ice thickness
map.  Let \f$H_{\text{tar}}\f$ be this target thickness,
and let \f$H\f$ be the current model thickness.  Recall that the mass continuity 
equation solved by IceModel::massContExplicitStep() is
  \f[ \frac{\partial H}{\partial t} = M - S - \nabla\cdot \mathbf{q} \f]
and that this procedure is supposed to produce \f$M\f$.
In this context, the semantics of \c -force_to_thk are that \f$M\f$ is modified
by a multiple of the difference between the target thickness and the current thickness.
In particular, the \f$\Delta M\f$ that is produced here is 
  \f[\Delta M = \alpha (H_{\text{tar}} - H)\f]
where \f$\alpha>0\f$ is determined below.  Note \f$\Delta M\f$ is positive in
areas where \f$H_{\text{tar}} > H\f$, so we are adding mass there, and we are ablating
in the other case.

Let \f$t_s\f$ be the start time and \f$t_e\f$ the end time for the run.
Without flow or basal mass balance, or any surface mass balance other than the
\f$\Delta M\f$ computed here, we are solving
  \f[ \frac{\partial H}{\partial t} = \alpha (H_{\text{tar}} - H) \f]
Let's assume \f$H(t_s)=H_0\f$.  This initial value problem has solution
\f$H(t) = H_{\text{tar}} + (H_0 - H_{\text{tar}}) e^{-\alpha (t-t_s)}\f$
and so
  \f[ H(t_e) = H_{\text{tar}} + (H_0 - H_{\text{tar}}) e^{-\alpha (t_e-t_s)} \f]
The constant \f$\alpha\f$ has a default value \c pism_config:force_to_thickness_alpha
of \f$0.002\,\text{a}^{-1}\f$.

The final feature is that we turn on this mechanism so it is harshest near the end
of the run.  In particular,
  \f[\Delta M = \lambda(t) \alpha (H_{\text{tar}} - H)\f]
where
  \f[\lambda(t) = \frac{t-t_s}{t_e-t_s}\f]

The next exacmple uses files generated from the EISMINT-Greenland experiment;
see the corresponding chapter of the User's Manual.

Suppose we regard the SSL2 run as a spin-up to reach a better temperature field.
It is a spinup in which the surface was allowed to evolve.  Assume the 
early file \c green20km_y1.nc has the target thickness, because it essentially
has the input thickness.  This script adds a 500 a run, to finalize the spinup,
in which the ice sheet geometry goes from the the thickness values in 
\c green_ssl2_110ka.nc to values very close to those in \c green20km_y1.nc:
\code
#!/bin/bash

NN=8

mpiexec -n $NN pismr -ys -500.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
  -surface pdd,forcing -pdd_fausto -force_to_thk green20km_y1.nc \
  -o with_force.nc -ts_file ts_with_force.nc -ts_times -500:10:0
  
mpiexec -n $NN pismr -ys -500.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
  -surface pdd -pdd_fausto \
  -o no_force.nc -ts_file ts_no_force.nc -ts_times -500:10:0

mpiexec -n $NN pismr -ys -500.0 -ye 0 -skip 5 -i green_ssl2_110ka.nc -atmosphere searise_greenland \
  -surface pdd,forcing -pdd_fausto -force_to_thk green20km_y1.nc -force_to_thk_alpha 0.0002 \
  -o weak_force.nc -ts_file ts_weak_force.nc -ts_times -500:10:0
\endcode
The script also has a run with no forcing, and one with forcing at a lower alpha value,
a factor of ten smaller than the default.

As shown below, the time series for \c ivol in the
above time series files show that the force-to-thickness mechanism is forcing
a system with negative feedback.  We see decaying oscillations toward the
intended volume.

\image html ivol_force_to_thk.png "\b Volume results from the -force_to_thk mechanism."
\anchor ivol_force_to_thk
 */
PetscErrorCode PSForceThickness::ice_surface_mass_flux(
       PetscReal t_years, PetscReal dt_years, IceModelVec2S &result) {
  PetscErrorCode ierr;

  // get the surface mass balance result from the next level up
  ierr = input_model->ice_surface_mass_flux(t_years, dt_years, result); CHKERRQ(ierr);

  ierr = verbPrintf(5, grid.com,
     "    updating surface mass balance using -force_to_thk mechanism ...");
     CHKERRQ(ierr);
    
  // force-to-thickness mechanism is only full-strength at end of run
  const PetscScalar lambda = (t_years - grid.start_year) / (grid.end_year - grid.start_year);
  ierr = verbPrintf(5, grid.com,
		    " (t_years = %.3f a, start_year = %.3f a, end_year = %.3f a, alpha = %.5f, lambda = %.3f)\n",
		    t_years, grid.start_year , grid.end_year, alpha, lambda); CHKERRQ(ierr);
  if ((lambda < 0.0) || (lambda > 1.0)) {
    SETERRQ(4,"computed lambda (for -force_to_thk) out of range; in updateSurfMassFluxAndProvide()");
  }

  PetscScalar **H;
  ierr = ice_thickness->get_array(H);   CHKERRQ(ierr);
  ierr = target_thickness.begin_access(); CHKERRQ(ierr);
  ierr = ftt_mask.begin_access(); CHKERRQ(ierr); 
  ierr = ftt_modified_acab.begin_access(); CHKERRQ(ierr); 
  ierr = result.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (ftt_mask(i,j) > 0.5) {
        result(i,j) += lambda * alpha * (target_thickness(i,j) - H[i][j]);
      }
      ftt_modified_acab(i,j) = result(i,j);
    }
  }
  ierr = ice_thickness->end_access(); CHKERRQ(ierr);
  ierr = target_thickness.end_access(); CHKERRQ(ierr);
  ierr = ftt_mask.end_access(); CHKERRQ(ierr); 
  ierr = ftt_modified_acab.end_access(); CHKERRQ(ierr); 
  ierr = result.end_access(); CHKERRQ(ierr);
  // no communication needed

  return 0;
}

//! Does not modify ice surface temperature.
PetscErrorCode PSForceThickness::ice_surface_temperature(PetscReal t_years, PetscReal dt_years,
							 IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = input_model->ice_surface_temperature(t_years, dt_years, result); CHKERRQ(ierr);

  return 0;
}

/*!
The timestep restriction is, by direct analogy, the same as for
   \f[\frac{dy}{dt} = - \alpha y\f]
with explicit (forward) Euler.  If \f$\Delta t\f$ is the time step then Euler is
\f$y_{n+1} = (1-\alpha \Delta t) y_n\f$.  We require for stability that
\f$|y_{n+1}|\le |y_n|\f$, which is to say \f$|1-\alpha \Delta t|\le 1\f$.
Equivalently (since \f$\alpha \Delta t>0\f$),
   \f[\alpha \Delta t\le 2\f]
Therefore we set here
   \f[\Delta t = \frac{2}{\alpha}.\f]
 */
PetscErrorCode PSForceThickness::max_timestep(PetscReal t_years, PetscReal &dt_years) {
  PetscErrorCode ierr;
  PetscReal max_dt = 2.0 / (alpha * secpera);
  
  ierr = input_model->max_timestep(t_years, dt_years); CHKERRQ(ierr);

  if (dt_years > 0) {
    if (max_dt > 0)
      dt_years = PetscMin(max_dt, dt_years);
  }
  else dt_years = max_dt;

  return 0;
}

PetscErrorCode PSForceThickness::write_model_state(PetscReal t_years, PetscReal dt_years,
						    string filename) {
  PetscErrorCode ierr;

  ierr = input_model->write_model_state(t_years, dt_years, filename); CHKERRQ(ierr);

  if (write_ftt_mask) {
    ierr = ftt_mask.write(filename.c_str()); CHKERRQ(ierr);
  }  

  return 0;
}

//! Adds ftt_modified_acab to "big" output files.
void PSForceThickness::add_vars_to_output(string key, set<string> &result) {
  if (key == "big")
    result.insert("ftt_modified_acab");
}

PetscErrorCode PSForceThickness::write_fields(set<string> vars, PetscReal t_years,
					      PetscReal dt_years, string filename) {
  PetscErrorCode ierr;

  ierr = input_model->write_fields(vars, t_years, dt_years, filename); CHKERRQ(ierr);

  if (vars.find("ftt_modified_acab") != vars.end()) {
    ierr = ftt_modified_acab.write(filename.c_str()); CHKERRQ(ierr); 
  }  

  if (vars.find("ftt_mask") != vars.end()) {
    ierr = ftt_mask.write(filename.c_str()); CHKERRQ(ierr); 
  }  

  return 0;
}


PetscErrorCode PSForceThickness::write_diagnostic_fields(PetscReal t_years, PetscReal dt_years,
							 string filename) {
  PetscErrorCode ierr;

  ierr = input_model->write_diagnostic_fields(t_years, dt_years, filename); CHKERRQ(ierr);

  if (write_ftt_mask) {
    ierr = ftt_mask.write(filename.c_str()); CHKERRQ(ierr);
  }

  ierr = ftt_modified_acab.write(filename.c_str()); CHKERRQ(ierr);

  return 0;
}
