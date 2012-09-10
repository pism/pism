// Copyright (C) 2011, 2012 PISM Authors
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
  mbscheme = NULL;
  faustogreve = NULL;
  base_ddf.snow = config.get("pdd_factor_snow");
  base_ddf.ice  = config.get("pdd_factor_ice");
  base_ddf.refreezeFrac = config.get("pdd_refreeze");
  base_pddStdDev = config.get("pdd_std_dev");
  base_pddThresholdTemp = config.get("pdd_positive_threshold_temp");

  pdd_annualize = false;
}

PSTemperatureIndex::~PSTemperatureIndex() {
  delete mbscheme;
  delete faustogreve;
}

PetscErrorCode PSTemperatureIndex::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool           pdd_rand, pdd_rand_repeatable, fausto_params, pSet;

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
    ierr = PISMOptionsIsSet("-pdd_annualize",
                            "Compute annual mass balance, removing yearly variations",
                            pdd_annualize); CHKERRQ(ierr);

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
    "* Initializing the default temperature-index, PDD-based surface processes scheme.\n"
    "  Precipitation and 2m air temperature provided by atmosphere are inputs.\n"
    "  Surface mass balance and ice upper surface temperature are outputs.\n"
    "  See PISM User's Manual for control of degree-day factors.\n");
    CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
    "  Computing number of positive degree-days by: "); CHKERRQ(ierr);
  if (pdd_rand_repeatable) {
    ierr = verbPrintf(2, grid.com, "simulation of a random process.\n"); CHKERRQ(ierr);
    mbscheme = new PDDrandMassBalance(config, true);
  } else if (pdd_rand) {
    ierr = verbPrintf(2, grid.com, "repeatable simulation of a random process.\n");
      CHKERRQ(ierr);
    mbscheme = new PDDrandMassBalance(config, false);
  } else {
    ierr = verbPrintf(2, grid.com, "an expectation integral.\n"); CHKERRQ(ierr);
    mbscheme = new PDDMassBalance(config);
  }

  if (config.get_flag("pdd_limit_timestep")) {
    ierr = verbPrintf(2, grid.com, "  NOTE: Limiting time-steps to 1 year.\n"); CHKERRQ(ierr);
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

  if ((config.get("pdd_std_dev_lapse_lat_rate") != 0.0) || fausto_params) {
    lat = dynamic_cast<IceModelVec2S*>(vars.get("latitude"));
    if (!lat)
      SETERRQ(grid.com, 10, "ERROR: 'latitude' is not available or is wrong type in dictionary");
  } else
    lat = NULL;


  if (fausto_params) {
    ierr = verbPrintf(2, grid.com,
       "  Setting PDD parameters from [Faustoetal2009] ...\n");
       CHKERRQ(ierr);

    //FIXME: this seems not to work because config is "const"?:  config.set("pdd_std_dev",2.53);
    base_pddStdDev = 2.53;

    lon = dynamic_cast<IceModelVec2S*>(vars.get("longitude"));
    if (!lon)
      SETERRQ(grid.com, 11, "ERROR: 'longitude' is not available or is wrong type in dictionary");
    usurf = dynamic_cast<IceModelVec2S*>(vars.get("usurf"));
    if (!usurf)
      SETERRQ(grid.com, 12, "ERROR: 'usurf' is not available or is wrong type in dictionary");

    faustogreve = new FaustoGrevePDDObject(grid, config);
  } else {
    // generally, this is the case in which degree day factors do not depend
    //   on location; we use base_ddf
    lon = NULL;
    usurf = NULL;
  }

  // if -pdd_annualize is set, update mass balance immediately (at the
  // beginning of the run)
  next_pdd_update = grid.time->current();

  ice_surface_temp.init_2d("ice_surface_temp", grid);
  ice_surface_temp.set_string("pism_intent", "diagnostic");
  ice_surface_temp.set_string("long_name",
                  "ice temperature at the ice surface");
  ierr = ice_surface_temp.set_units("K"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSTemperatureIndex::max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict) {
  PetscErrorCode ierr;

  if (pdd_annualize) {
    if (PetscAbs(my_t - next_pdd_update) < 1e-12)
      my_dt = convert(1.0, "years", "seconds");
    else
      my_dt = next_pdd_update - my_t;
  } else {
    my_dt = -1;
  }

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
    if (pdd_annualize && PetscAbs(my_t + my_dt - next_pdd_update) < 1)
      my_dt = next_pdd_update - my_t;
  } else
    restrict = false;

  return 0;
}

PetscErrorCode PSTemperatureIndex::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;

  PetscReal one_year = convert(1.0, "years", "seconds");

  if ((fabs(my_t - t) < 1e-12) &&
      (fabs(my_dt - dt) < 1e-12))
    return 0;

  t  = my_t;
  dt = my_dt;

  // This flag is set in pclimate to allow testing the model. It does not
  // affect the normal operation of PISM.
  if (config.get_flag("pdd_limit_timestep")) {
    my_dt = PetscMin(my_dt, one_year);
  }

  if (pdd_annualize) {
    if (my_t + my_dt > next_pdd_update) {
      ierr = verbPrintf(3, grid.com,
                        "  Updating mass balance for one year starting at %s ...\n",
                        grid.time->date(my_t).c_str());
      ierr = update_internal(my_t, one_year); CHKERRQ(ierr);
      next_pdd_update = my_t + one_year;
    }
  } else {
    ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PSTemperatureIndex::update_internal(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;

  // to ensure that temperature time series are correct:
  ierr = atmosphere->update(my_t, my_dt); CHKERRQ(ierr);

  // This is a point-wise (local) computation, so we can use "climatic_mass_balance" to store
  // precipitation:
  ierr = atmosphere->mean_precipitation(climatic_mass_balance); CHKERRQ(ierr);

  // set up air temperature time series
  PetscInt Nseries;
  ierr = mbscheme->getNForTemperatureSeries(my_t, my_dt, Nseries); CHKERRQ(ierr);

  // time since the beginning of the year, in seconds
  const PetscScalar dtseries = my_dt / ((PetscScalar) (Nseries - 1));

  // times for the air temperature time-series, in years:
  vector<PetscScalar> ts(Nseries), T(Nseries);
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

  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {

      // the temperature time series from the PISMAtmosphereModel and its modifiers
      ierr = atmosphere->temp_time_series(i, j, Nseries, &ts[0], &T[0]); CHKERRQ(ierr);

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
      PetscScalar snow_amount = mbscheme->getSnowFromPrecipAndTemperatureTimeSeries(
                                  climatic_mass_balance(i,j), // precipitation rate (input)
                                  my_t, dtseries, &T[0], Nseries);

      // use degree-day factors, and number of PDDs, and the snow precipitation, to
      //   get surface mass balance (and diagnostics: accumulation, melt, runoff)
      ierr = mbscheme->getMassFluxesFromPDDs(ddf, my_dt, pddsum, snow_amount,
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

  ierr = PISMSurfaceModel::define_variables(vars, nc, nctype); CHKERRQ(ierr);

  return 0;

}

PetscErrorCode PSTemperatureIndex::write_variables(set<string> vars, string filename) {
  PetscErrorCode ierr;

  if (set_contains(vars, "ice_surface_temp")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "ice_surface_temp", false); CHKERRQ(ierr);
    ierr = tmp.set_metadata(ice_surface_temp, 0); CHKERRQ(ierr);

    ierr = ice_surface_temperature(tmp); CHKERRQ(ierr);

    ierr = tmp.write(filename.c_str()); CHKERRQ(ierr);
    vars.erase("ice_surface_temp");
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    ierr = climatic_mass_balance.write(filename.c_str()); CHKERRQ(ierr);
    vars.erase("climatic_mass_balance");
  }

  if (set_contains(vars, "saccum")) {
    ierr = accumulation_rate.write(filename.c_str()); CHKERRQ(ierr);
    vars.erase("saccum");
  }

  if (set_contains(vars, "smelt")) {
    ierr = melt_rate.write(filename.c_str()); CHKERRQ(ierr);
    vars.erase("smelt");
  }

  if (set_contains(vars, "srunoff")) {
    ierr = runoff_rate.write(filename.c_str()); CHKERRQ(ierr);
    vars.erase("srunoff");
  }

  ierr = PISMSurfaceModel::write_variables(vars, filename); CHKERRQ(ierr);

  return 0;
}
