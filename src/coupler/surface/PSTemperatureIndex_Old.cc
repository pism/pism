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
#include "error_handling.hh"

namespace pism {

///// PISM surface model implementing a PDD scheme.

PSTemperatureIndex_Old::PSTemperatureIndex_Old(const IceGrid &g)
  : SurfaceModel(g), temperature_name("ice_surface_temp"),
    ice_surface_temp(g.config.get_unit_system(), temperature_name, g) {
  mbscheme = NULL;
  faustogreve = NULL;
  base_ddf.snow = m_config.get("pdd_factor_snow");
  base_ddf.ice  = m_config.get("pdd_factor_ice");
  base_ddf.refreezeFrac = m_config.get("pdd_refreeze");
  base_pddStdDev = m_config.get("pdd_std_dev");
  base_pddThresholdTemp = m_config.get("pdd_positive_threshold_temp");

  pdd_annualize = false;

  mass_balance_name = "climatic_mass_balance";

  climatic_mass_balance.create(m_grid, mass_balance_name, WITHOUT_GHOSTS);
  climatic_mass_balance.set_attrs("diagnostic",
                                  "instantaneous ice-equivalent surface mass balance (accumulation/ablation) rate",
                                  "kg m-2 s-1",
                                  "land_ice_surface_specific_mass_balance_flux");
  climatic_mass_balance.set_glaciological_units("kg m-2 year-1");
  climatic_mass_balance.write_in_glaciological_units = true;
  climatic_mass_balance.metadata().set_string("comment", "positive values correspond to ice gain");

  // diagnostic fields:

  accumulation_rate.create(m_grid, "saccum", WITHOUT_GHOSTS);
  accumulation_rate.set_attrs("diagnostic",
                              "instantaneous ice-equivalent surface accumulation rate (precip minus rain)",
                              "m s-1",
                              "");
  accumulation_rate.set_glaciological_units("m year-1");
  accumulation_rate.write_in_glaciological_units = true;

  melt_rate.create(m_grid, "smelt", WITHOUT_GHOSTS);
  melt_rate.set_attrs("diagnostic",
                      "instantaneous ice-equivalent surface melt rate",
                      "m s-1",
                      "");
  melt_rate.set_glaciological_units("m year-1");
  melt_rate.write_in_glaciological_units = true;

  runoff_rate.create(m_grid, "srunoff", WITHOUT_GHOSTS);
  runoff_rate.set_attrs("diagnostic",
                        "instantaneous ice-equivalent surface meltwater runoff rate",
                        "m s-1",
                        "");
  runoff_rate.set_glaciological_units("m year-1");
  runoff_rate.write_in_glaciological_units = true;
}

PSTemperatureIndex_Old::~PSTemperatureIndex_Old() {
  delete mbscheme;
  delete faustogreve;
}

void PSTemperatureIndex_Old::init() {
  bool           pdd_rand, pdd_rand_repeatable, fausto_params, pSet;

  SurfaceModel::init();

  {
    pdd_rand = options::Bool("-pdd_rand",
                             "Use a PDD implementation based on simulating a random process");
    pdd_rand_repeatable = options::Bool("-pdd_rand_repeatable",
                                        "Use a PDD implementation based on simulating a repeatable random process");
    fausto_params = options::Bool("-pdd_fausto",
                                  "Set PDD parameters using formulas (6) and (7) in [Faustoetal2009]");
    pdd_annualize = options::Bool("-pdd_annualize",
                                  "Compute annual mass balance, removing yearly variations");

    OptionsReal("-pdd_factor_snow", "PDD snow factor",
                base_ddf.snow, pSet);
    OptionsReal("-pdd_factor_ice", "PDD ice factor",
                base_ddf.ice, pSet);
    OptionsReal("-pdd_refreeze", "PDD refreeze fraction",
                base_ddf.refreezeFrac, pSet);

    OptionsReal("-pdd_std_dev", "PDD standard deviation",
                base_pddStdDev, pSet);
    OptionsReal("-pdd_positive_threshold_temp",
                "PDD uses this temp in K to determine 'positive' temperatures",
                base_pddThresholdTemp, pSet);
  }

  verbPrintf(2, m_grid.com,
             "* Initializing the default temperature-index, PDD-based surface processes scheme.\n"
             "  Precipitation and 2m air temperature provided by atmosphere are inputs.\n"
             "  Surface mass balance and ice upper surface temperature are outputs.\n"
             "  See PISM User's Manual for control of degree-day factors.\n");

  std::string method = "the method reported above.";
  if (mbscheme == NULL) {
    if (pdd_rand_repeatable) {
      method = "simulation of a random process.";
      mbscheme = new PDDrandMassBalance_Old(m_config, true);
    } else if (pdd_rand) {
      method = "repeatable simulation of a random process.";
      mbscheme = new PDDrandMassBalance_Old(m_config, false);
    } else {
      method = "an expectation integral.";
      mbscheme = new PDDMassBalance_Old(m_config);
    }
  }
  verbPrintf(2, m_grid.com,
             "  Computing number of positive degree-days by: %s\n",
             method.c_str());


  if ((m_config.get("pdd_std_dev_lapse_lat_rate") != 0.0) || fausto_params) {
    lat = m_grid.variables().get_2d_scalar("latitude");
  } else
    lat = NULL;


  if (fausto_params) {
    verbPrintf(2, m_grid.com,
               "  Setting PDD parameters from [Faustoetal2009] ...\n");

    //FIXME: this seems not to work because config is "const"?:  config.set("pdd_std_dev",2.53);
    base_pddStdDev = 2.53;

    lon   = m_grid.variables().get_2d_scalar("longitude");
    usurf = m_grid.variables().get_2d_scalar("usurf");

  if (faustogreve == NULL) {
    faustogreve = new FaustoGrevePDDObject_Old(m_grid);
  }

  } else {
    // generally, this is the case in which degree day factors do not depend
    //   on location; we use base_ddf
    lon = NULL;
    usurf = NULL;
  }

  // if -pdd_annualize is set, update mass balance immediately (at the
  // beginning of the run)
  next_pdd_update = m_grid.time->current();

  ice_surface_temp.set_string("pism_intent", "diagnostic");
  ice_surface_temp.set_string("long_name",
                  "ice temperature at the ice surface");
  ice_surface_temp.set_units("K");
}

void PSTemperatureIndex_Old::max_timestep(PetscReal my_t, PetscReal &my_dt, bool &restrict) {

  if (pdd_annualize) {
    if (fabs(my_t - next_pdd_update) < 1e-12) {
      my_dt = m_grid.convert(1.0, "years", "seconds");
    } else {
      my_dt = next_pdd_update - my_t;
    }
  } else {
    my_dt = -1;
  }

  PetscReal dt_atmosphere;
  atmosphere->max_timestep(my_t, dt_atmosphere, restrict);

  if (restrict) {
    if (my_dt > 0) {
      my_dt = std::min(my_dt, dt_atmosphere);
    } else {
      my_dt = dt_atmosphere;
    }
  }

  if (my_dt > 0) {
    restrict = true;
  } else {
    restrict = false;
  }
}

void PSTemperatureIndex_Old::update(PetscReal my_t, PetscReal my_dt) {

  PetscReal one_year = m_grid.convert(1.0, "years", "seconds");

  if ((fabs(my_t - m_t) < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = my_t;
  m_dt = my_dt;

  if (pdd_annualize) {
    if (my_t + my_dt > next_pdd_update) {
      verbPrintf(3, m_grid.com,
                 "  Updating mass balance for one year starting at %1.1f...\n",
                 m_grid.time->date(my_t).c_str());
      update_internal(my_t, one_year);
      next_pdd_update = my_t + one_year;
    }
  } else {
    update_internal(my_t, my_dt);
  }
}

void PSTemperatureIndex_Old::update_internal(PetscReal my_t, PetscReal my_dt) {

  const double ice_density = m_config.get("ice_density");

  // to ensure that temperature time series are correct:
  atmosphere->update(my_t, my_dt);

  // This is a point-wise (local) computation, so we can use "climatic_mass_balance" to store
  // precipitation:
  atmosphere->mean_precipitation(climatic_mass_balance);

  // set up air temperature time series
  PetscInt Nseries = 0;
  mbscheme->getNForTemperatureSeries(my_t, my_dt, Nseries);

  PetscReal one_year = m_grid.convert(1.0, "years", "seconds");

  // time since the beginning of the year, in seconds
  const PetscScalar tseries = m_grid.time->mod(my_t, one_year),
    dtseries = my_dt / ((PetscScalar) (Nseries - 1));

  // times for the air temperature time-series, in years:
  std::vector<PetscScalar> ts(Nseries), T(Nseries);
  for (PetscInt k = 0; k < Nseries; ++k) {
    ts[k] = my_t + k * dtseries;
  }

  IceModelVec::AccessList list;

  if (lat != NULL) {
    list.add(*lat);
  }

  if (faustogreve != NULL) {
    if (lat == NULL) {
      throw RuntimeError("faustogreve object is allocated BUT lat==NULL");
    }
    if (lon == NULL) {
      throw RuntimeError("faustogreve object is allocated BUT lon==NULL");
    }
    if (usurf == NULL) {
      throw RuntimeError("faustogreve object is allocated BUT usurf==NULL");
    }
    list.add(*lon);
    list.add(*usurf);
    faustogreve->update_temp_mj(usurf, lat, lon);
  }

  const PetscScalar sigmalapserate = m_config.get("pdd_std_dev_lapse_lat_rate"),
                    sigmabaselat   = m_config.get("pdd_std_dev_lapse_lat_base");
  if (sigmalapserate != 0.0) {
    if (lat == NULL) {
      throw RuntimeError("pdd_std_dev_lapse_lat_rate is nonzero BUT lat==NULL");
    }
  }

  DegreeDayFactors_Old  ddf = base_ddf;

  atmosphere->begin_pointwise_access();
  list.add(climatic_mass_balance);

  list.add(accumulation_rate);
  list.add(melt_rate);
  list.add(runoff_rate);

  atmosphere->init_timeseries(ts);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // the temperature time series from the AtmosphereModel and its modifiers
    atmosphere->temp_time_series(i, j, T);

    if (faustogreve != NULL) {
      // we have been asked to set mass balance parameters according to
      //   formula (6) in [\ref Faustoetal2009]; they overwrite ddf set above
      faustogreve->setDegreeDayFactors(i, j, (*usurf)(i, j),
                                       (*lat)(i, j), (*lon)(i, j), ddf);
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
    PetscScalar snow_amount = mbscheme->getSnowFromPrecipAndTemperatureTimeSeries(climatic_mass_balance(i,j), // precipitation rate (input)
                                                                                  tseries, dtseries, &T[0], Nseries);

    // use degree-day factors, and number of PDDs, and the snow precipitation, to
    //   get surface mass balance (and diagnostics: accumulation, melt, runoff)
    mbscheme->getMassFluxesFromPDDs(ddf, my_dt, pddsum, snow_amount,
                                    accumulation_rate(i,j), // output
                                    melt_rate(i,j), // output
                                    runoff_rate(i,j), // output
                                    climatic_mass_balance(i,j));

    // convert from m/s to m/s * kg/m3 = (kg m-2)/s
    climatic_mass_balance(i,j) *= ice_density;
  }

  atmosphere->end_pointwise_access();
}


void PSTemperatureIndex_Old::ice_surface_mass_flux(IceModelVec2S &result) {

  climatic_mass_balance.copy_to(result);
}


void PSTemperatureIndex_Old::ice_surface_temperature(IceModelVec2S &result) {
  atmosphere->mean_annual_temp(result);
}

void PSTemperatureIndex_Old::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  if (keyword == "medium" || keyword == "big") {
    result.insert(mass_balance_name);
    result.insert(temperature_name);
  }

  if (keyword == "big") {
    result.insert("saccum");
    result.insert("smelt");
    result.insert("srunoff");
  }

  atmosphere->add_vars_to_output(keyword, result);
}

void PSTemperatureIndex_Old::define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {

  SurfaceModel::define_variables(vars, nc, nctype);

  if (set_contains(vars, temperature_name)) {
    ice_surface_temp.define(nc, nctype, true);
  }

  if (set_contains(vars, mass_balance_name)) {
    climatic_mass_balance.define(nc, nctype);
  }

  if (set_contains(vars, "saccum")) {
    accumulation_rate.define(nc, nctype);
  }

  if (set_contains(vars, "smelt")) {
    melt_rate.define(nc, nctype);
  }

  if (set_contains(vars, "srunoff")) {
    runoff_rate.define(nc, nctype);
  }
}

void PSTemperatureIndex_Old::write_variables(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, temperature_name)) {
    IceModelVec2S tmp;
    tmp.create(m_grid, temperature_name, WITHOUT_GHOSTS);
    tmp.metadata(0) = ice_surface_temp;

    ice_surface_temperature(tmp);

    tmp.write(nc);
    vars.erase(temperature_name);
  }

  if (set_contains(vars, mass_balance_name)) {
    climatic_mass_balance.write(nc);
    vars.erase(mass_balance_name);
  }

  if (set_contains(vars, "saccum")) {
    accumulation_rate.write(nc);
    vars.erase("saccum");
  }

  if (set_contains(vars, "smelt")) {
    melt_rate.write(nc);
    vars.erase("smelt");
  }

  if (set_contains(vars, "srunoff")) {
    runoff_rate.write(nc);
    vars.erase("srunoff");
  }

  SurfaceModel::write_variables(vars, nc);
}

} // end of namespace pism
