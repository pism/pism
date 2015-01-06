// Copyright (C) 2011, 2012, 2013, 2014, 2015 PISM Authors
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

#include "PSTemperatureIndex.hh"
#include "localMassBalance.hh"
#include "IceGrid.hh"
#include "pism_options.hh"
#include "PISMVars.hh"
#include "PISMTime.hh"
#include "PISMAtmosphere.hh"
#include "Mask.hh"

#include <algorithm>            // std::min
#include <cassert>

#include "error_handling.hh"

namespace pism {

///// PISM surface model implementing a PDD scheme.

PSTemperatureIndex::PSTemperatureIndex(const IceGrid &g)
  : SurfaceModel(g),
    ice_surface_temp(g.config.get_unit_system(), "ice_surface_temp", m_grid) {

  mbscheme              = NULL;
  faustogreve           = NULL;
  sd_period             = 0;
  sd_ref_time           = 0.0;
  base_ddf.snow         = m_config.get("pdd_factor_snow");
  base_ddf.ice          = m_config.get("pdd_factor_ice");
  base_ddf.refreezeFrac = m_config.get("pdd_refreeze");
  base_pddThresholdTemp = m_config.get("pdd_positive_threshold_temp");
  base_pddStdDev        = m_config.get("pdd_std_dev");
  sd_use_param          = m_config.get_flag("pdd_std_dev_use_param");
  sd_param_a            = m_config.get("pdd_std_dev_param_a");
  sd_param_b            = m_config.get("pdd_std_dev_param_b");


  randomized = options::Bool("-pdd_rand",
                             "Use a PDD implementation based on simulating a random process");
  randomized_repeatable = options::Bool("-pdd_rand_repeatable",
                                        "Use a PDD implementation based on simulating a"
                                        " repeatable random process");
  fausto_params = options::Bool("-pdd_fausto",
                                "Set PDD parameters using formulas (6) and (7)"
                                " in [Faustoetal2009]");

  options::String file("-pdd_sd_file", "Read standard deviation from file");
  sd_file_set = file.is_set();

  options::Integer period("-pdd_sd_period",
                          "Length of the standard deviation data period in years", 0);
  sd_period = period;

  options::Integer sd_ref_year("-pdd_sd_reference_year",
                               "Standard deviation data reference year", 0);

  if (randomized_repeatable) {
    mbscheme = new PDDrandMassBalance(m_config, true);
  } else if (randomized) {
    mbscheme = new PDDrandMassBalance(m_config, false);
  } else {
    mbscheme = new PDDMassBalance(m_config);
  }

  if (fausto_params) {
    faustogreve = new FaustoGrevePDDObject(m_grid);
  }

  if (sd_ref_year.is_set()) {
    sd_ref_time = m_grid.convert(sd_ref_year, "years", "seconds");
  }

  if (sd_file_set) {
    // find out how many records there are in the file and set the
    // air_temp_sd buffer size

    unsigned int n_records = 0;
    std::string short_name = "air_temp_sd";
    unsigned int buffer_size = (unsigned int) m_config.get("climate_forcing_buffer_size");

    PIO nc(m_grid.com, "netcdf3", m_grid.config.get_unit_system());
    nc.open(file, PISM_READONLY);
    n_records = nc.inq_nrecords(short_name, "");
    nc.close();

    // If -..._period is not set, make ..._n_records the minimum of the
    // buffer size and the number of available records. Otherwise try
    // to keep all available records in memory.
    if (sd_period == 0) {
      n_records = std::min(n_records, buffer_size);
    }

    if (n_records < 1) {
      throw RuntimeError::formatted("Can't find '%s' in %s.",
                                    short_name.c_str(), file->c_str());
    }

    air_temp_sd.set_n_records(n_records);

  } else {
    // using constant standard deviation, so set buffer size to 1
    air_temp_sd.set_n_records(1);
  }

  air_temp_sd.create(m_grid, "air_temp_sd", false);
  air_temp_sd.set_attrs("climate_forcing",
                        "standard deviation of near-surface air temperature",
                        "Kelvin", "");

  climatic_mass_balance.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  climatic_mass_balance.set_attrs("diagnostic",
                                  "instantaneous surface mass balance (accumulation/ablation) rate",
                                  "kg m-2 s-1",
                                  "land_ice_surface_specific_mass_balance_flux");
  climatic_mass_balance.set_glaciological_units("kg m-2 year-1");
  climatic_mass_balance.write_in_glaciological_units = true;
  climatic_mass_balance.metadata().set_string("comment", "positive values correspond to ice gain");

  // diagnostic fields:

  accumulation_rate.create(m_grid, "saccum", WITHOUT_GHOSTS);
  accumulation_rate.set_attrs("diagnostic",
                              "instantaneous surface accumulation rate"
                              " (precipitation minus rain)",
                              "kg m-2 s-1",
                              "");
  accumulation_rate.set_glaciological_units("kg m-2 year-1");
  accumulation_rate.write_in_glaciological_units = true;

  melt_rate.create(m_grid, "smelt", WITHOUT_GHOSTS);
  melt_rate.set_attrs("diagnostic",
                      "instantaneous surface melt rate",
                      "kg m-2 s-1",
                      "");
  melt_rate.set_glaciological_units("kg m-2 year-1");
  melt_rate.write_in_glaciological_units = true;

  runoff_rate.create(m_grid, "srunoff", WITHOUT_GHOSTS);
  runoff_rate.set_attrs("diagnostic",
                        "instantaneous surface meltwater runoff rate",
                        "kg m-2 s-1",
                        "");
  runoff_rate.set_glaciological_units("kg m-2 year-1");
  runoff_rate.write_in_glaciological_units = true;

  snow_depth.create(m_grid, "snow_depth", WITHOUT_GHOSTS);
  snow_depth.set_attrs("diagnostic",
                       "snow cover depth (set to zero once a year)",
                       "m", "");
  snow_depth.set(0.0);

  ice_surface_temp.set_string("pism_intent", "diagnostic");
  ice_surface_temp.set_string("long_name",
                              "ice temperature at the ice surface");
  ice_surface_temp.set_units("K");
}

PSTemperatureIndex::~PSTemperatureIndex() {
  delete mbscheme;
  delete faustogreve;
}

void PSTemperatureIndex::init() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  SurfaceModel::init();

  verbPrintf(2, m_grid.com,
             "* Initializing the default temperature-index, PDD-based surface processes scheme.\n"
             "  Precipitation and 2m air temperature provided by atmosphere are inputs.\n"
             "  Surface mass balance and ice upper surface temperature are outputs.\n"
             "  See PISM User's Manual for control of degree-day factors.\n");

  verbPrintf(2, m_grid.com,
             "  Computing number of positive degree-days by: ");
  if (randomized) {
    verbPrintf(2, m_grid.com, "simulation of a random process.\n");
  } else if (randomized_repeatable) {
    verbPrintf(2, m_grid.com, "repeatable simulation of a random process.\n");
  } else {
    verbPrintf(2, m_grid.com, "an expectation integral.\n");
  }

  mask = m_grid.variables().get_2d_mask("mask");

  if ((m_config.get("pdd_std_dev_lapse_lat_rate") != 0.0) || fausto_params) {
    lat = m_grid.variables().get_2d_scalar("latitude");
  } else {
    lat = NULL;
  }

  if (fausto_params) {
    verbPrintf(2, m_grid.com,
               "  Setting PDD parameters from [Faustoetal2009] ...\n");

    base_pddStdDev = 2.53;

    lon   = m_grid.variables().get_2d_scalar("longitude");
    usurf = m_grid.variables().get_2d_scalar("usurf");
  } else {
    // generally, this is the case in which degree day factors do not depend
    //   on location; we use base_ddf
    lon   = NULL;
    usurf = NULL;
  }

  options::String file("-pdd_sd_file", "Read standard deviation from file");
  if (file.is_set()) {
    verbPrintf(2, m_grid.com,
               "  Reading standard deviation of near-surface temperature from '%s'...\n",
               file->c_str());
    air_temp_sd.init(file, sd_period, sd_ref_time);
  } else {
    verbPrintf(2, m_grid.com,
               "  Option -pdd_sd_file is not set. Using a constant value.\n");
    air_temp_sd.init_constant(base_pddStdDev);
  }

  std::string input_file;
  bool do_regrid = false;
  int start = -1;

  // find PISM input file to read data from:
  find_pism_input(input_file, do_regrid, start);

  // read snow precipitation rate from file
  verbPrintf(2, m_grid.com,
             "    reading snow depth (ice equivalent meters) from %s ... \n",
             input_file.c_str());
  snow_depth.regrid(input_file, OPTIONAL, 0.0);

  m_next_balance_year_start = compute_next_balance_year_start(m_grid.time->current());
}

void PSTemperatureIndex::max_timestep(double my_t, double &my_dt, bool &restrict) {
  atmosphere->max_timestep(my_t, my_dt, restrict);
}

double PSTemperatureIndex::compute_next_balance_year_start(double time) {
    // compute the time corresponding to the beginning of the next balance year
    double
      balance_year_start_day = m_config.get("pdd_balance_year_start_day"),
      one_day                = m_grid.convert(1.0, "days", "seconds"),
      year_start             = m_grid.time->calendar_year_start(time),
      balance_year_start     = year_start + (balance_year_start_day - 1.0) * one_day;

    if (balance_year_start > time) {
      return balance_year_start;
    }
    return m_grid.time->increment_date(balance_year_start, 1);
}


void PSTemperatureIndex::update(double my_t, double my_dt) {

  if ((fabs(my_t - m_t) < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = my_t;
  m_dt = my_dt;

  // update to ensure that temperature and precipitation time series
  // are correct:
  atmosphere->update(my_t, my_dt);

  // set up air temperature and precipitation time series
  int Nseries = mbscheme->get_timeseries_length(my_dt);

  const double dtseries = my_dt / Nseries;
  std::vector<double> ts(Nseries), T(Nseries), S(Nseries), P(Nseries), PDDs(Nseries);
  for (int k = 0; k < Nseries; ++k) {
    ts[k] = my_t + k * dtseries;
  }

  // update standard deviation time series
  if (sd_file_set == true) {
    air_temp_sd.update(my_t, my_dt);
    air_temp_sd.init_interpolation(ts);
  }

  MaskQuery m(*mask);

  IceModelVec::AccessList list(*mask);

  if (lat != NULL) {
    list.add(*lat);
  }

  if (faustogreve != NULL) {
    assert(lat != NULL && lon != NULL && usurf != NULL);
    list.add(*lon);
    list.add(*usurf);
    faustogreve->update_temp_mj(*usurf, *lat, *lon);
  }

  const double sigmalapserate = m_config.get("pdd_std_dev_lapse_lat_rate"),
    sigmabaselat   = m_config.get("pdd_std_dev_lapse_lat_base");
  if (sigmalapserate != 0.0) {
    assert(lat != NULL);
  }

  DegreeDayFactors  ddf = base_ddf;

  atmosphere->begin_pointwise_access();
  list.add(air_temp_sd);
  list.add(climatic_mass_balance);

  list.add(accumulation_rate);
  list.add(melt_rate);
  list.add(runoff_rate);
  list.add(snow_depth);

  atmosphere->init_timeseries(ts);

  const double ice_density = m_config.get("ice_density");

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    // the temperature time series from the AtmosphereModel and its modifiers
    atmosphere->temp_time_series(i, j, T);

    // the precipitation time series from AtmosphereModel and its modifiers
    atmosphere->precip_time_series(i, j, P);

    // interpolate temperature standard deviation time series
    if (sd_file_set == true) {
      air_temp_sd.interp(i, j, S);
    } else {
      for (int k = 0; k < Nseries; ++k) {
        S[k] = air_temp_sd(i, j);
      }
    }

    if (faustogreve != NULL) {
      // we have been asked to set mass balance parameters according to
      //   formula (6) in [\ref Faustoetal2009]; they overwrite ddf set above
      faustogreve->setDegreeDayFactors(i, j, (*usurf)(i, j),
                                       (*lat)(i, j), (*lon)(i, j), ddf);
    }

    // apply standard deviation lapse rate on top of prescribed values
    if (sigmalapserate != 0.0) {
      for (int k = 0; k < Nseries; ++k) {
        S[k] += sigmalapserate * ((*lat)(i,j) - sigmabaselat);
      }
      air_temp_sd(i, j) = S[0]; // ensure correct SD reporting
    }

    // apply standard deviation param over ice if in use
    if (sd_use_param && m.icy(i,j)) {
      for (int k = 0; k < Nseries; ++k) {
        S[k] = sd_param_a * (T[k] - 273.15) + sd_param_b;
        if (S[k] < 0.0) {
          S[k] = 0.0 ;
        }
      }
      air_temp_sd(i, j) = S[0]; // ensure correct SD reporting
    }

    // Use temperature time series, the "positive" threshhold, and
    // the standard deviation of the daily variability to get the
    // number of positive degree days (PDDs)
    mbscheme->get_PDDs(&S[0], dtseries, &T[0], Nseries, &PDDs[0]);

    // Use temperature time series to remove rainfall from precipitation
    mbscheme->get_snow_accumulation(&P[0], // precipitation rate (input-output)
                                    &T[0], // air temperature (input)
                                    Nseries);

    // Use degree-day factors, and number of PDDs, and the snow
    // precipitation, to get surface mass balance (and diagnostics:
    // accumulation, melt, runoff)
    {
      double next_snow_depth_reset = m_next_balance_year_start;
      accumulation_rate(i,j)     = 0.0;
      melt_rate(i,j)             = 0.0;
      runoff_rate(i,j)           = 0.0;
      climatic_mass_balance(i,j) = 0.0;
      for (int k = 0; k < Nseries; ++k) {
        if (ts[k] >= next_snow_depth_reset) {
          snow_depth(i,j)       = 0.0;
          while (next_snow_depth_reset <= ts[k]) {
            next_snow_depth_reset = m_grid.time->increment_date(next_snow_depth_reset, 1);
          }
        }

        double accumulation     = P[k] * dtseries;
        accumulation_rate(i,j) += accumulation;

        mbscheme->step(ddf, PDDs[k], accumulation,
                       snow_depth(i,j), melt_rate(i,j), runoff_rate(i,j),
                       climatic_mass_balance(i,j));
      }
      // convert from [m during the current time-step] to kg m-2 s-1
      accumulation_rate(i,j)     *= (ice_density/m_dt);
      melt_rate(i,j)             *= (ice_density/m_dt);
      runoff_rate(i,j)           *= (ice_density/m_dt);
      climatic_mass_balance(i,j) *= (ice_density/m_dt);
    }

    if (m.ocean(i,j)) {
      snow_depth(i,j) = 0.0;  // snow over the ocean does not stick
    }
  }

  atmosphere->end_pointwise_access();

  m_next_balance_year_start = compute_next_balance_year_start(m_grid.time->current());
}


void PSTemperatureIndex::ice_surface_mass_flux(IceModelVec2S &result) {
  climatic_mass_balance.copy_to(result);
}


void PSTemperatureIndex::ice_surface_temperature(IceModelVec2S &result) {

  atmosphere->mean_annual_temp(result);
}

void PSTemperatureIndex::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {

  SurfaceModel::add_vars_to_output(keyword, result);

  result.insert("snow_depth");

  if (keyword == "medium" || keyword == "big") {
    result.insert("climatic_mass_balance");
    result.insert("ice_surface_temp");
  }

  if (keyword == "big") {
    result.insert("air_temp_sd");
    result.insert("saccum");
    result.insert("smelt");
    result.insert("srunoff");
  }
}

void PSTemperatureIndex::define_variables(const std::set<std::string> &vars, const PIO &nc, IO_Type nctype) {

  if (set_contains(vars, "ice_surface_temp")) {
    ice_surface_temp.define(nc, nctype, true);
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    climatic_mass_balance.define(nc, nctype);
  }

  if (set_contains(vars, "air_temp_sd")) {
    air_temp_sd.define(nc, nctype);
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

  if (set_contains(vars, "snow_depth")) {
    snow_depth.define(nc, nctype);
  }

  SurfaceModel::define_variables(vars, nc, nctype);
}

void PSTemperatureIndex::write_variables(const std::set<std::string> &vars_input, const PIO &nc) {
  std::set<std::string> vars = vars_input;

  if (set_contains(vars, "ice_surface_temp")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "ice_surface_temp", WITHOUT_GHOSTS);
    tmp.metadata() = ice_surface_temp;

    ice_surface_temperature(tmp);

    tmp.write(nc);
    vars.erase("ice_surface_temp");
  }

  if (set_contains(vars, "climatic_mass_balance")) {
    climatic_mass_balance.write(nc);
    vars.erase("climatic_mass_balance");
  }

  if (set_contains(vars, "air_temp_sd")) {
    air_temp_sd.average(m_t, m_dt);
    air_temp_sd.write(nc);
    vars.erase("air_temp_sd");
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

  if (set_contains(vars, "snow_depth")) {
    snow_depth.write(nc);
    vars.erase("snow_depth");
  }

  SurfaceModel::write_variables(vars, nc);
}

} // end of namespace pism
