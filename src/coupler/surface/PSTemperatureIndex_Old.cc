// Copyright (C) 2011, 2012, 2014, 2015, 2016 PISM Authors
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
#include "base/util/IceGrid.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMTime.hh"
#include "base/util/PISMVars.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_options.hh"
#include "coupler/PISMAtmosphere.hh"
#include "localMassBalance_old.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/pism_utilities.hh"

namespace pism {
namespace surface {

///// PISM surface model implementing a PDD scheme.

TemperatureIndex_Old::TemperatureIndex_Old(IceGrid::ConstPtr g)
  : SurfaceModel(g), temperature_name("ice_surface_temp"),
    ice_surface_temp(m_sys, temperature_name) {
  mbscheme = NULL;
  faustogreve = NULL;
  base_ddf.snow = m_config->get_double("pdd_factor_snow");
  base_ddf.ice  = m_config->get_double("pdd_factor_ice");
  base_ddf.refreezeFrac = m_config->get_double("pdd_refreeze");
  base_pddStdDev = m_config->get_double("pdd_std_dev");
  base_pddThresholdTemp = m_config->get_double("pdd_positive_threshold_temp");

  pdd_annualize = false;

  mass_balance_name = "climatic_mass_balance";

  climatic_mass_balance.create(m_grid, mass_balance_name, WITHOUT_GHOSTS);
  climatic_mass_balance.set_attrs("diagnostic",
                                  "instantaneous ice-equivalent surface mass balance"
                                  " (accumulation/ablation) rate",
                                  "kg m-2 s-1",
                                  "land_ice_surface_specific_mass_balance_flux");
  climatic_mass_balance.metadata().set_string("glaciological_units", "kg m-2 year-1");
  climatic_mass_balance.write_in_glaciological_units = true;
  climatic_mass_balance.metadata().set_string("comment", "positive values correspond to ice gain");

  // diagnostic fields:

  accumulation_rate.create(m_grid, "saccum", WITHOUT_GHOSTS);
  accumulation_rate.set_attrs("diagnostic",
                              "instantaneous ice-equivalent surface accumulation rate (precip minus rain)",
                              "m s-1",
                              "");
  accumulation_rate.metadata().set_string("glaciological_units", "m year-1");
  accumulation_rate.write_in_glaciological_units = true;

  melt_rate.create(m_grid, "smelt", WITHOUT_GHOSTS);
  melt_rate.set_attrs("diagnostic",
                      "instantaneous ice-equivalent surface melt rate",
                      "m s-1",
                      "");
  melt_rate.metadata().set_string("glaciological_units", "m year-1");
  melt_rate.write_in_glaciological_units = true;

  runoff_rate.create(m_grid, "srunoff", WITHOUT_GHOSTS);
  runoff_rate.set_attrs("diagnostic",
                        "instantaneous ice-equivalent surface meltwater runoff rate",
                        "m s-1",
                        "");
  runoff_rate.metadata().set_string("glaciological_units", "m year-1");
  runoff_rate.write_in_glaciological_units = true;
}

TemperatureIndex_Old::~TemperatureIndex_Old() {
  delete mbscheme;
  delete faustogreve;
}

void TemperatureIndex_Old::init_impl() {
  bool pdd_rand, pdd_rand_repeatable, fausto_params;

  // call the default implementation (not the interface method init())
  SurfaceModel::init_impl();

  {
    pdd_rand = options::Bool("-pdd_rand",
                             "Use a PDD implementation based on simulating a random process");
    pdd_rand_repeatable = options::Bool("-pdd_rand_repeatable",
                                        "Use a PDD implementation based on simulating a repeatable random process");
    fausto_params = options::Bool("-pdd_fausto",
                                  "Set PDD parameters using formulas (6) and (7) in [Faustoetal2009]");
    pdd_annualize = options::Bool("-pdd_annualize",
                                  "Compute annual mass balance, removing yearly variations");

    base_ddf.snow = options::Real("-pdd_factor_snow", "PDD snow factor",
                                  base_ddf.snow);
    base_ddf.ice = options::Real("-pdd_factor_ice", "PDD ice factor",
                                 base_ddf.ice);
    base_ddf.refreezeFrac = options::Real("-pdd_refreeze", "PDD refreeze fraction",
                                          base_ddf.refreezeFrac);

    base_pddStdDev = options::Real("-pdd_std_dev", "PDD standard deviation",
                                   base_pddStdDev);
    base_pddThresholdTemp = options::Real("-pdd_positive_threshold_temp",
                                          "PDD uses this temp in K to determine 'positive' temperatures",
                                          base_pddThresholdTemp);
  }

  m_log->message(2,
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
  m_log->message(2,
             "  Computing number of positive degree-days by: %s\n",
             method.c_str());


  if (fausto_params) {
    m_log->message(2,
               "  Setting PDD parameters from [Faustoetal2009] ...\n");

    base_pddStdDev = 2.53;

    if (faustogreve == NULL) {
      faustogreve = new FaustoGrevePDDObject_Old(m_grid);
    }

  }

  // if -pdd_annualize is set, update mass balance immediately (at the
  // beginning of the run)
  next_pdd_update = m_grid->ctx()->time()->current();

  ice_surface_temp.set_string("pism_intent", "diagnostic");
  ice_surface_temp.set_string("long_name",
                  "ice temperature at the ice surface");
  ice_surface_temp.set_string("units", "K");
}

MaxTimestep TemperatureIndex_Old::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
}

MaxTimestep TemperatureIndex_Old::max_timestep(PetscReal my_t) {

  MaxTimestep max_dt;
  if (pdd_annualize) {
    if (fabs(my_t - next_pdd_update) < 1e-12) {
      max_dt = MaxTimestep(units::convert(m_sys, 1.0, "years", "seconds"));
    } else {
      max_dt = MaxTimestep(next_pdd_update - my_t);
    }
  }

  MaxTimestep dt_atmosphere = m_atmosphere->max_timestep(my_t);

  return std::min(max_dt, dt_atmosphere);
}

void TemperatureIndex_Old::update_impl(PetscReal my_t, PetscReal my_dt) {

  PetscReal one_year = units::convert(m_sys, 1.0, "years", "seconds");

  if ((fabs(my_t - m_t) < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = my_t;
  m_dt = my_dt;

  if (pdd_annualize) {
    if (my_t + my_dt > next_pdd_update) {
      m_log->message(3,
                     "  Updating mass balance for one year starting at %s...\n",
                     m_grid->ctx()->time()->date(my_t).c_str());
      update_internal(my_t, one_year);
      next_pdd_update = my_t + one_year;
    }
  } else {
    update_internal(my_t, my_dt);
  }
}

void TemperatureIndex_Old::update_internal(PetscReal my_t, PetscReal my_dt) {

  const double ice_density = m_config->get_double("ice_density");

  // to ensure that temperature time series are correct:
  m_atmosphere->update(my_t, my_dt);

  // This is a point-wise (local) computation, so we can use "climatic_mass_balance" to store
  // precipitation:
  m_atmosphere->mean_precipitation(climatic_mass_balance);

  // set up air temperature time series
  int Nseries = 0;
  mbscheme->getNForTemperatureSeries(my_t, my_dt, Nseries);

  PetscReal one_year = units::convert(m_sys, 1.0, "years", "seconds");

  // time since the beginning of the year, in seconds
  const PetscScalar tseries = m_grid->ctx()->time()->mod(my_t, one_year),
    dtseries = my_dt / ((PetscScalar) (Nseries - 1));

  // times for the air temperature time-series, in years:
  std::vector<PetscScalar> ts(Nseries), T(Nseries);
  for (int k = 0; k < Nseries; ++k) {
    ts[k] = my_t + k * dtseries;
  }

  const IceModelVec2S *lat = NULL, *lon = NULL, *usurf = NULL;

  IceModelVec::AccessList list;

  if (faustogreve != NULL) {
    lon   = m_grid->variables().get_2d_scalar("longitude");
    usurf = m_grid->variables().get_2d_scalar("usurf");

    list.add(*lon);
    list.add(*usurf);
    faustogreve->update_temp_mj(*usurf, *lat, *lon);
  }

  const PetscScalar sigmalapserate = m_config->get_double("pdd_std_dev_lapse_lat_rate"),
                    sigmabaselat   = m_config->get_double("pdd_std_dev_lapse_lat_base");
  if (sigmalapserate != 0.0) {
    lat = m_grid->variables().get_2d_scalar("latitude");
    list.add(*lat);
  }

  LocalMassBalance_Old::DegreeDayFactors  ddf = base_ddf;

  m_atmosphere->begin_pointwise_access();
  list.add(climatic_mass_balance);

  list.add(accumulation_rate);
  list.add(melt_rate);
  list.add(runoff_rate);

  m_atmosphere->init_timeseries(ts);

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // the temperature time series from the AtmosphereModel and its modifiers
      m_atmosphere->temp_time_series(i, j, T);

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
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_atmosphere->end_pointwise_access();
}


void TemperatureIndex_Old::ice_surface_mass_flux_impl(IceModelVec2S &result) {

  result.copy_from(climatic_mass_balance);
}


void TemperatureIndex_Old::ice_surface_temperature_impl(IceModelVec2S &result) {
  m_atmosphere->mean_annual_temp(result);
}

void TemperatureIndex_Old::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
  if (keyword == "medium" || keyword == "big" || keyword == "2dbig") {
    result.insert(mass_balance_name);
    result.insert(temperature_name);
  }

  if (keyword == "big" || keyword == "2dbig") {
    result.insert("saccum");
    result.insert("smelt");
    result.insert("srunoff");
  }

  m_atmosphere->add_vars_to_output(keyword, result);
}

void TemperatureIndex_Old::define_variables_impl(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {

  SurfaceModel::define_variables_impl(vars, nc, nctype);

  if (set_contains(vars, temperature_name)) {
    std::string order = m_config->get_string("output_variable_order");
    io::define_spatial_variable(ice_surface_temp, *m_grid, nc, nctype, order, true);
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

void TemperatureIndex_Old::write_variables_impl(const std::set<std::string> &vars_input, const PIO &nc) {
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

  SurfaceModel::write_variables_impl(vars, nc);
}

} // end of namespace surface
} // end of namespace pism
