// Copyright (C) 2011, 2012, 2014 PISM Authors
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

#include "PSTemperatureIndex_Old.hh"
#include "localMassBalance_old.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include "PISMVars.hh"
#include "PISMTime.hh"
#include "PISMAtmosphere.hh"
#include "PISMConfig.hh"

///// PISM surface model implementing a PDD scheme.

PSTemperatureIndex_Old::PSTemperatureIndex_Old(IceGrid &g, const PISMConfig &conf)
  : PISMSurfaceModel(g, conf), artm(g.get_unit_system()) {
  mbscheme = NULL;
  faustogreve = NULL;
  base_ddf.snow = config.get("pdd_factor_snow");
  base_ddf.ice  = config.get("pdd_factor_ice");
  base_ddf.refreezeFrac = config.get("pdd_refreeze");
  base_pddStdDev = config.get("pdd_std_dev");
  base_pddThresholdTemp = config.get("pdd_positive_threshold_temp");

  pdd_annualize = false;
}

PSTemperatureIndex_Old::~PSTemperatureIndex_Old() {
  delete mbscheme;
  delete faustogreve;
}

PetscErrorCode PSTemperatureIndex_Old::init(PISMVars &vars) {
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
    mbscheme = new PDDrandMassBalance_Old(config, true);
  } else if (pdd_rand) {
    ierr = verbPrintf(2, grid.com, "repeatable simulation of a random process.\n");
      CHKERRQ(ierr);
    mbscheme = new PDDrandMassBalance_Old(config, false);
  } else {
    ierr = verbPrintf(2, grid.com, "an expectation integral.\n"); CHKERRQ(ierr);
    mbscheme = new PDDMassBalance_Old(config);
  }

  ierr = acab.create(grid, "acab", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = acab.set_attrs("diagnostic",
			"instantaneous ice-equivalent surface mass balance (accumulation/ablation) rate",
			"m s-1",  // m *ice-equivalent* per second
			"land_ice_surface_specific_mass_balance");  // CF standard_name
                        CHKERRQ(ierr);
  ierr = acab.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  acab.write_in_glaciological_units = true;
  acab.metadata().set_string("comment", "positive values correspond to ice gain"); CHKERRQ(ierr);

  // diagnostic fields:

  ierr = accumulation_rate.create(grid, "saccum", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = accumulation_rate.set_attrs("diagnostic",
                                     "instantaneous ice-equivalent surface accumulation rate (precip minus rain)",
                                     "m s-1",
                                     ""); CHKERRQ(ierr);
  ierr = accumulation_rate.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  accumulation_rate.write_in_glaciological_units = true;

  ierr = melt_rate.create(grid, "smelt", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = melt_rate.set_attrs("diagnostic",
                             "instantaneous ice-equivalent surface melt rate",
                             "m s-1",
                             ""); CHKERRQ(ierr);
  ierr = melt_rate.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  melt_rate.write_in_glaciological_units = true;

  ierr = runoff_rate.create(grid, "srunoff", WITHOUT_GHOSTS); CHKERRQ(ierr);
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

    faustogreve = new FaustoGrevePDDObject_Old(grid, config);
  } else {
    // generally, this is the case in which degree day factors do not depend
    //   on location; we use base_ddf
    lon = NULL;
    usurf = NULL;
  }

  // if -pdd_annualize is set, update mass balance immediately (at the
  // beginning of the run)
  next_pdd_update = grid.time->current();

  artm.init_2d("artm", grid);
  artm.set_string("pism_intent", "diagnostic");
  artm.set_string("long_name",
                  "ice temperature at the ice surface");
  ierr = artm.set_units("K"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PSTemperatureIndex_Old::max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict) {
  PetscErrorCode ierr;

  if (pdd_annualize) {
    if (PetscAbs(my_t - next_pdd_update) < 1e-12)
      my_dt = grid.convert(1.0, "years", "seconds");
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

  if (my_dt > 0)
    restrict = true;
  else
    restrict = false;

  return 0;
}

PetscErrorCode PSTemperatureIndex_Old::update(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;

  PetscReal one_year = grid.convert(1.0, "years", "seconds");

  if ((fabs(my_t - m_t) < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12))
    return 0;

  m_t  = my_t;
  m_dt = my_dt;

  if (pdd_annualize) {
    if (my_t + my_dt > next_pdd_update) {
      ierr = verbPrintf(3, grid.com,
                        "  Updating mass balance for one year starting at %1.1f...\n",
                        grid.time->date(my_t).c_str());
      ierr = update_internal(my_t, one_year); CHKERRQ(ierr);
      next_pdd_update = my_t + one_year;
    }
  } else {
    ierr = update_internal(my_t, my_dt); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PSTemperatureIndex_Old::update_internal(PetscReal my_t, PetscReal my_dt) {
  PetscErrorCode ierr;

  const double ice_density = config.get("ice_density");

  // to ensure that temperature time series are correct:
  ierr = atmosphere->update(my_t, my_dt); CHKERRQ(ierr);

  // This is a point-wise (local) computation, so we can use "acab" to store
  // precipitation:
  ierr = atmosphere->mean_precipitation(acab); CHKERRQ(ierr);

  // set up air temperature time series
  PetscInt Nseries = 0;
  ierr = mbscheme->getNForTemperatureSeries(my_t, my_dt, Nseries); CHKERRQ(ierr);

  PetscReal one_year = grid.convert(1.0, "years", "seconds");

  // time since the beginning of the year, in seconds
  const PetscScalar tseries = grid.time->mod(my_t, one_year),
    dtseries = my_dt / ((PetscScalar) (Nseries - 1));

  // times for the air temperature time-series, in years:
  std::vector<PetscScalar> ts(Nseries), T(Nseries);
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

  DegreeDayFactors_Old  ddf = base_ddf;

  ierr = atmosphere->begin_pointwise_access(); CHKERRQ(ierr);
  ierr = acab.begin_access(); CHKERRQ(ierr);

  ierr = accumulation_rate.begin_access(); CHKERRQ(ierr);
  ierr = melt_rate.begin_access(); CHKERRQ(ierr);
  ierr = runoff_rate.begin_access(); CHKERRQ(ierr);

  ierr = atmosphere->init_timeseries(&ts[0], Nseries); CHKERRQ(ierr);

  for (PetscInt i = grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j<grid.ys+grid.ym; ++j) {

      // the temperature time series from the PISMAtmosphereModel and its modifiers
      ierr = atmosphere->temp_time_series(i, j, &T[0]); CHKERRQ(ierr);

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
      PetscScalar pddsum = mbscheme->getPDDSumFromTemperatureTimeSeries(sigma, base_pddThresholdTemp,
                                                                        tseries, dtseries, &T[0], Nseries);

      // use the temperature time series to remove the rainfall from the precipitation
      PetscScalar snow_amount = mbscheme->getSnowFromPrecipAndTemperatureTimeSeries(
                                  acab(i,j), // precipitation rate (input)
                                  tseries, dtseries, &T[0], Nseries);

      // use degree-day factors, and number of PDDs, and the snow precipitation, to
      //   get surface mass balance (and diagnostics: accumulation, melt, runoff)
      ierr = mbscheme->getMassFluxesFromPDDs(ddf, my_dt, pddsum, snow_amount,
                                             accumulation_rate(i,j), // output
                                             melt_rate(i,j), // output
                                             runoff_rate(i,j), // output
                                             acab(i,j)); // acab = smb (output)
      CHKERRQ(ierr);

      // convert from m/s to m/s * kg/m3 = (kg m-2)/s
      acab(i,j) *= ice_density;
    }
  }

  ierr = accumulation_rate.end_access(); CHKERRQ(ierr);
  ierr = melt_rate.end_access(); CHKERRQ(ierr);
  ierr = runoff_rate.end_access(); CHKERRQ(ierr);

  ierr = acab.end_access(); CHKERRQ(ierr);
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


PetscErrorCode PSTemperatureIndex_Old::ice_surface_mass_flux(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = acab.copy_to(result); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode PSTemperatureIndex_Old::ice_surface_temperature(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = atmosphere->mean_annual_temp(result); CHKERRQ(ierr);

  return 0;
}

void PSTemperatureIndex_Old::add_vars_to_output(std::string keyword, std::set<std::string> &result) {
  if (keyword == "medium" || keyword == "big") {
    result.insert("acab");
    result.insert("artm");
  }

  if (keyword == "big") {
    result.insert("saccum");
    result.insert("smelt");
    result.insert("srunoff");
  }

  atmosphere->add_vars_to_output(keyword, result);
}

PetscErrorCode PSTemperatureIndex_Old::define_variables(std::set<std::string> vars, const PIO &nc, PISM_IO_Type nctype) {
  PetscErrorCode ierr;

  ierr = PISMSurfaceModel::define_variables(vars, nc, nctype); CHKERRQ(ierr);

  if (set_contains(vars, "artm")) {
    ierr = artm.define(nc, nctype, true); CHKERRQ(ierr);
  }

  if (set_contains(vars, "acab")) {
    ierr = acab.define(nc, nctype); CHKERRQ(ierr);
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

  return 0;
}

PetscErrorCode PSTemperatureIndex_Old::write_variables(std::set<std::string> vars, const PIO &nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, "artm")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "artm", WITHOUT_GHOSTS); CHKERRQ(ierr);
    tmp.metadata(0) = artm;

    ierr = ice_surface_temperature(tmp); CHKERRQ(ierr);

    ierr = tmp.write(nc); CHKERRQ(ierr);
    vars.erase("artm");
  }

  if (set_contains(vars, "acab")) {
    ierr = acab.write(nc); CHKERRQ(ierr);
    vars.erase("acab");
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

  ierr = PISMSurfaceModel::write_variables(vars, nc); CHKERRQ(ierr);

  return 0;
}
