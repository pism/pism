// Copyright (C) 2011, 2012, 2013 PISM Authors
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

#include "PSTemperatureIndex.hh"
#include "localMassBalance.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include "PISMVars.hh"
#include "PISMTime.hh"
#include "PISMAtmosphere.hh"

///// PISM surface model implementing a PDD scheme.

PSTemperatureIndex::PSTemperatureIndex(IceGrid &g, const NCConfigVariable &conf)
  : PISMSurfaceModel(g, conf) {
  PetscErrorCode ierr = allocate_PSTemperatureIndex(); CHKERRCONTINUE(ierr);
  if (ierr != 0)
    PISMEnd();

}

PSTemperatureIndex::~PSTemperatureIndex() {
  delete mbscheme;
  delete faustogreve;
}

PetscErrorCode PSTemperatureIndex::allocate_PSTemperatureIndex() {
  PetscErrorCode ierr;
  bool flag;

  mbscheme		= NULL;
  faustogreve		= NULL;
  base_ddf.snow		= config.get("pdd_factor_snow");
  base_ddf.ice		= config.get("pdd_factor_ice");
  base_ddf.refreezeFrac = config.get("pdd_refreeze");
  base_pddStdDev	= config.get("pdd_std_dev");
  base_pddThresholdTemp = config.get("pdd_positive_threshold_temp");

  ierr = PetscOptionsBegin(grid.com, "",
                           "Temperature-index (PDD) scheme for surface (snow) processes", "");
                           CHKERRQ(ierr);
  {
    ierr = PISMOptionsIsSet("-pdd_rand",
                            "Use a PDD implementation based on simulating a random process",
			    randomized); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-pdd_rand_repeatable",
                            "Use a PDD implementation based on simulating a repeatable random process",
			    randomized_repeatable); CHKERRQ(ierr);
    ierr = PISMOptionsIsSet("-pdd_fausto",
                            "Set PDD parameters using formulas (6) and (7) in [Faustoetal2009]",
			    fausto_params); CHKERRQ(ierr);

    ierr = PISMOptionsReal("-pdd_factor_snow", "PDD snow factor",
                           base_ddf.snow, flag); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-pdd_factor_ice", "PDD ice factor",
                           base_ddf.ice, flag); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-pdd_refreeze", "PDD refreeze fraction",
                           base_ddf.refreezeFrac, flag); CHKERRQ(ierr);

    ierr = PISMOptionsReal("-pdd_std_dev", "PDD standard deviation",
                           base_pddStdDev, flag); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-pdd_positive_threshold_temp",
                           "PDD uses this temp in K to determine 'positive' temperatures",
                           base_pddThresholdTemp, flag); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (randomized_repeatable) {
    mbscheme = new PDDrandMassBalance(config, true);
  } else if (randomized) {
    mbscheme = new PDDrandMassBalance(config, false);
  } else {
    mbscheme = new PDDMassBalance(config);
  }

  if (fausto_params) {
    faustogreve = new FaustoGrevePDDObject(grid, config);
  }

  ierr = climatic_mass_balance.create(grid, "climatic_mass_balance", false); CHKERRQ(ierr);
  ierr = climatic_mass_balance.set_attrs("diagnostic",
					 "instantaneous ice-equivalent surface mass balance (accumulation/ablation) rate",
					 "m s-1",  // m *ice-equivalent* per second
					 "land_ice_surface_specific_mass_balance");  // CF standard_name
  CHKERRQ(ierr);
  ierr = climatic_mass_balance.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  climatic_mass_balance.write_in_glaciological_units = true;
  ierr = climatic_mass_balance.set_attr("comment", "positive values correspond to ice gain"); CHKERRQ(ierr);

  // diagnostic fields:

  ierr = accumulation_rate.create(grid, "saccum", false); CHKERRQ(ierr);
  ierr = accumulation_rate.set_attrs("diagnostic",
                                     "instantaneous ice-equivalent surface accumulation rate"
                                     " (precipitation minus rain)",
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

  ierr = snow_depth.create(grid, "snow_depth", false); CHKERRQ(ierr);
  ierr = snow_depth.set_attrs("diagnostic",
			      "snow cover depth (set to zero once a year)",
			      "m", ""); CHKERRQ(ierr);
  ierr = snow_depth.set(0.0); CHKERRQ(ierr);

  ice_surface_temp.init_2d("ice_surface_temp", grid);
  ice_surface_temp.set_string("pism_intent", "diagnostic");
  ice_surface_temp.set_string("long_name",
                  "ice temperature at the ice surface");
  ierr = ice_surface_temp.set_units("K"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSTemperatureIndex::init(PISMVars &vars) {
  PetscErrorCode ierr;

  t = dt = GSL_NAN;  // every re-init restarts the clock

  ierr = PISMSurfaceModel::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
    "* Initializing the default temperature-index, PDD-based surface processes scheme.\n"
    "  Precipitation and 2m air temperature provided by atmosphere are inputs.\n"
    "  Surface mass balance and ice upper surface temperature are outputs.\n"
    "  See PISM User's Manual for control of degree-day factors.\n");
    CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
    "  Computing number of positive degree-days by: "); CHKERRQ(ierr);
  if (randomized_repeatable) {
    ierr = verbPrintf(2, grid.com, "simulation of a random process.\n"); CHKERRQ(ierr);
  } else if (randomized) {
    ierr = verbPrintf(2, grid.com, "repeatable simulation of a random process.\n"); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, "an expectation integral.\n"); CHKERRQ(ierr);
  }

  if (config.get_flag("pdd_limit_timestep")) {
    ierr = verbPrintf(2, grid.com, "  NOTE: Limiting time-steps to 1 year.\n"); CHKERRQ(ierr);
  }

  if ((config.get("pdd_std_dev_lapse_lat_rate") != 0.0) || fausto_params) {
    lat = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
    if (!lat)
      SETERRQ(grid.com, 10, "ERROR: 'latitude' is not available or is wrong type in dictionary");
  } else {
    lat = NULL;
  }

  if (fausto_params) {
    ierr = verbPrintf(2, grid.com,
       "  Setting PDD parameters from [Faustoetal2009] ...\n");
       CHKERRQ(ierr);

    base_pddStdDev = 2.53;

    lon = dynamic_cast<IceModelVec2S*>(vars.get("longitude"));
    if (!lon)
      SETERRQ(grid.com, 11, "ERROR: 'longitude' is not available or is wrong type in dictionary");

    usurf = dynamic_cast<IceModelVec2S*>(vars.get("usurf"));
    if (!usurf)
      SETERRQ(grid.com, 12, "ERROR: 'usurf' is not available or is wrong type in dictionary");
  } else {
    // generally, this is the case in which degree day factors do not depend
    //   on location; we use base_ddf
    lon = NULL;
    usurf = NULL;
  }

  string input_file;
  bool regrid = false;
  int start = -1;
  
  // find PISM input file to read data from:
  ierr = find_pism_input(input_file, regrid, start); CHKERRQ(ierr);

  // read snow precipitation rate from file
  ierr = verbPrintf(2, grid.com,
		    "    reading snow depth (ice equivalent meters) from %s ... \n",
		    input_file.c_str()); CHKERRQ(ierr);
  ierr = snow_depth.regrid(input_file.c_str(), 0.0); CHKERRQ(ierr);
  
  return 0;
}

PetscErrorCode PSTemperatureIndex::max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict) {
  PetscErrorCode ierr;

  // compute the time corresponding to the beginning of the next balance year
  PetscReal one_year = convert(1.0, "years", "seconds"),
    one_day = convert(1.0, "days", "seconds"),
    year_start = my_t - grid.time->mod(my_t, one_year),
    balance_year_start = year_start + (config.get("pdd_balance_year_start_day") - 1.0) * one_day,
    next_balance_year_start = balance_year_start > my_t ? balance_year_start : balance_year_start + one_year;
  
  if (PetscAbs(my_t - next_balance_year_start) < 1e-12)
    my_dt = convert(1.0, "years", "seconds");
  else
    my_dt = next_balance_year_start - my_t;

  PetscReal dt_atmosphere;
  ierr = atmosphere->max_timestep(my_t, dt_atmosphere, restrict); CHKERRQ(ierr);

  if (restrict) {
    if (my_dt > 0)
      my_dt = PetscMin(my_dt, dt_atmosphere);
    else
      my_dt = dt_atmosphere;
  }

  if (my_dt > 0) {
    restrict = true;

    // try to avoid very small time steps:
    // (Necessary when driving this PDD model with monthly data, for example.)
    if (PetscAbs(my_t + my_dt - next_balance_year_start) < 1)
      my_dt = next_balance_year_start - my_t;
  } else {
    restrict = false;
  }

  return 0;
}

PetscErrorCode PSTemperatureIndex::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;

  if ((fabs(my_t - t) < 1e-12) &&
      (fabs(my_dt - dt) < 1e-12))
    return 0;

  t  = my_t;
  dt = my_dt;

  // This flag is set in pclimate to allow testing the model. It does not
  // affect the normal operation of PISM.
  if (config.get_flag("pdd_limit_timestep")) {
    PetscReal one_year = convert(1.0, "years", "seconds");
    my_dt = PetscMin(my_dt, one_year);
  }

  ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSTemperatureIndex::update_internal(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;
  PetscReal year_fraction = grid.time->year_fraction(my_t),
    balance_year_start_year_fraction = config.get("pdd_balance_year_start_day") / 365.0;

  if (fabs(year_fraction - balance_year_start_year_fraction) < 0.01) { // about 4 days
    ierr = verbPrintf(3, grid.com, "  PDD model: Re-setting snow depth to 0 meters.\n"); CHKERRQ(ierr);
    ierr = snow_depth.set(0.0); CHKERRQ(ierr);
  }

  // upate to ensure that temperature and precipitation time series
  // are correct:
  ierr = atmosphere->update(my_t, my_dt); CHKERRQ(ierr);

  // set up air temperature time series
  PetscInt Nseries;
  ierr = mbscheme->getNForTemperatureSeries(my_t, my_dt, Nseries); CHKERRQ(ierr);

  const PetscScalar dtseries = my_dt / ((PetscScalar) (Nseries - 1));
  vector<PetscScalar> ts(Nseries), T(Nseries), P(Nseries);
  for (PetscInt k = 0; k < Nseries; ++k)
    ts[k] = my_t + k * dtseries;

  if (lat != NULL) {
    ierr = lat->begin_access(); CHKERRQ(ierr);
  }

  if (faustogreve != NULL) {
    if (lat == NULL) { SETERRQ(grid.com, 1,"faustogreve object is allocated BUT lat==NULL"); }
    if (lon == NULL) { SETERRQ(grid.com, 2,"faustogreve object is allocated BUT lon==NULL"); }
    if (usurf == NULL) { SETERRQ(grid.com, 3,"faustogreve object is allocated BUT usurf==NULL"); }
    ierr = lon->begin_access(); CHKERRQ(ierr);
    ierr = usurf->begin_access(); CHKERRQ(ierr);
    ierr = faustogreve->update_temp_mj(usurf, lat, lon); CHKERRQ(ierr);
  }

  const PetscScalar sigmalapserate = config.get("pdd_std_dev_lapse_lat_rate"),
                    sigmabaselat   = config.get("pdd_std_dev_lapse_lat_base");
  if (sigmalapserate != 0.0) {
    if (lat == NULL) { SETERRQ(grid.com, 4,"pdd_std_dev_lapse_lat_rate is nonzero BUT lat==NULL"); }
  }

  DegreeDayFactors  ddf = base_ddf;

  ierr = atmosphere->begin_pointwise_access(); CHKERRQ(ierr);
  ierr = climatic_mass_balance.begin_access(); CHKERRQ(ierr);

  ierr = accumulation_rate.begin_access(); CHKERRQ(ierr);
  ierr = melt_rate.begin_access(); CHKERRQ(ierr);
  ierr = runoff_rate.begin_access(); CHKERRQ(ierr);
  ierr = snow_depth.begin_access(); CHKERRQ(ierr);

  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {

      // the temperature time series from the PISMAtmosphereModel and its modifiers
      ierr = atmosphere->temp_time_series(i, j, Nseries, &ts[0], &T[0]); CHKERRQ(ierr);

      // the precipitation time series from PISMAtmosphereModel and its modifiers
      ierr = atmosphere->precip_time_series(i, j, Nseries, &ts[0], &P[0]); CHKERRQ(ierr);

      if (faustogreve != NULL) {
	// we have been asked to set mass balance parameters according to
	//   formula (6) in [\ref Faustoetal2009]; they overwrite ddf set above
	ierr = faustogreve->setDegreeDayFactors(i, j, (*usurf)(i, j),
                                                (*lat)(i, j), (*lon)(i, j), ddf);
        CHKERRQ(ierr);
      }

      // use the temperature time series, the "positive" threshhold, and the
      //   standard deviation of the daily variability to get the number of
      //   positive degree days (PDDs)
      PetscScalar sigma = base_pddStdDev;
      if (sigmalapserate != 0.0) {
        sigma += sigmalapserate * ((*lat)(i,j) - sigmabaselat);
      }
      PetscScalar pddsum = mbscheme->getPDDSumFromTemperatureTimeSeries(
                                  sigma, base_pddThresholdTemp,
                                  my_t, dtseries, &T[0], Nseries);

      // use the temperature time series to remove the rainfall from the precipitation
      PetscScalar snow_accumulation = mbscheme->getSnowFromPrecipAndTemperatureTimeSeries(
                                  &P[0], // precipitation rate (input)
                                  my_t, dtseries, &T[0], Nseries);

      // use degree-day factors, and number of PDDs, and the snow precipitation, to
      //   get surface mass balance (and diagnostics: accumulation, melt, runoff)
      ierr = mbscheme->getMassFluxesFromPDDs(ddf, my_dt, pddsum,
					     snow_accumulation,
					     snow_depth(i,j), // input-output
                                             accumulation_rate(i,j), // output
                                             melt_rate(i,j), // output
                                             runoff_rate(i,j), // output
                                             climatic_mass_balance(i,j)); // climatic_mass_balance = smb (output)
                                             CHKERRQ(ierr);
    }
  }

  ierr = accumulation_rate.end_access(); CHKERRQ(ierr);
  ierr = melt_rate.end_access(); CHKERRQ(ierr);
  ierr = runoff_rate.end_access(); CHKERRQ(ierr);
  ierr = snow_depth.end_access(); CHKERRQ(ierr);

  ierr = climatic_mass_balance.end_access(); CHKERRQ(ierr);
  ierr = atmosphere->end_pointwise_access(); CHKERRQ(ierr);

  if (lat != NULL) {
    ierr = lat->end_access(); CHKERRQ(ierr);
  }

  if (faustogreve != NULL) {
    ierr = lon->end_access(); CHKERRQ(ierr);
    ierr = usurf->end_access(); CHKERRQ(ierr);
  }

  return 0;
}


PetscErrorCode PSTemperatureIndex::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = climatic_mass_balance.copy_to(result); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PSTemperatureIndex::ice_surface_temperature(IceModelVec2S &result) {

  PetscErrorCode ierr = atmosphere->mean_annual_temp(result); CHKERRQ(ierr);

  return 0;
}

void PSTemperatureIndex::add_vars_to_output(string keyword, map<string,NCSpatialVariable> &result) {

  PISMSurfaceModel::add_vars_to_output(keyword, result);

  result["snow_depth"] = snow_depth.get_metadata();
  
  if (keyword == "medium" || keyword == "big") {
    result["climatic_mass_balance"] = climatic_mass_balance.get_metadata();
    result["ice_surface_temp"] = ice_surface_temp;
  }

  if (keyword == "big") {
    result["saccum"]  = accumulation_rate.get_metadata();
    result["smelt"]   = melt_rate.get_metadata();
    result["srunoff"] = runoff_rate.get_metadata();
  }
}

PetscErrorCode PSTemperatureIndex::define_variables(set<string> vars, const PIO &nc, PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "ice_surface_temp")) {
    ierr = ice_surface_temp.define(nc, nctype, true); CHKERRQ(ierr);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    ierr = climatic_mass_balance.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "saccum")) {
    ierr = accumulation_rate.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "smelt")) {
    ierr = melt_rate.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "srunoff")) {
    ierr = runoff_rate.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "snow_depth")) {
    ierr = snow_depth.define(nc, nctype); CHKERRQ(ierr);
  }

  ierr = PISMSurfaceModel::define_variables(vars, nc, nctype); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode PSTemperatureIndex::write_variables(set<string> vars, const PIO &nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, "ice_surface_temp")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "ice_surface_temp", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(ice_surface_temp, 0); CHKERRQ(ierr);

    ierr = ice_surface_temperature(tmp); CHKERRQ(ierr);

    ierr = tmp.write(nc); CHKERRQ(ierr);
    vars.erase("ice_surface_temp");
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    ierr = climatic_mass_balance.write(nc); CHKERRQ(ierr);
    vars.erase("climatic_mass_balance");
  }

  if (set_contains(vars, "saccum")) {
    ierr = accumulation_rate.write(nc); CHKERRQ(ierr);
    vars.erase("saccum");
  }

  if (set_contains(vars, "smelt")) {
    ierr = melt_rate.write(nc); CHKERRQ(ierr);
    vars.erase("smelt");
  }

  if (set_contains(vars, "srunoff")) {
    ierr = runoff_rate.write(nc); CHKERRQ(ierr);
    vars.erase("srunoff");
  }

  if (set_contains(vars, "snow_depth")) {
    ierr = snow_depth.write(nc); CHKERRQ(ierr);
    vars.erase("snow_depth");
  }

  ierr = PISMSurfaceModel::write_variables(vars, nc); CHKERRQ(ierr);

  return 0;
}
