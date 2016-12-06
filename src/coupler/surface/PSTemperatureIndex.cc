// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#include <algorithm>            // std::min
#include <cassert>
#include <gsl/gsl_math.h>

#include "PSTemperatureIndex.hh"
#include "localMassBalance.hh"
#include "base/util/IceGrid.hh"
#include "base/util/pism_options.hh"
#include "base/util/PISMVars.hh"
#include "base/util/PISMTime.hh"
#include "coupler/PISMAtmosphere.hh"
#include "base/util/Mask.hh"
#include "base/util/io/PIO.hh"

#include "base/util/error_handling.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/pism_utilities.hh"
#include "base/util/IceModelVec2CellType.hh"

namespace pism {
namespace surface {

// Diagnostic classes.

/*! @brief Surface accumulation averaged over reporting intervals. */
class PDD_saccum_average : public Diag_average<TemperatureIndex>
{
public:
  PDD_saccum_average(const TemperatureIndex *m);
protected:
  const IceModelVec2S &cumulative_value() const;
};

/*! @brief Surface melt averaged over reporting intervals. */
class PDD_smelt_average : public Diag_average<TemperatureIndex>
{
public:
  PDD_smelt_average(const TemperatureIndex *m);
protected:
  const IceModelVec2S &cumulative_value() const;
};

/*! @brief Surface runoff averaged over reporting intervals. */
class PDD_srunoff_average : public Diag_average<TemperatureIndex>
{
public:
  PDD_srunoff_average(const TemperatureIndex *m);
protected:
  const IceModelVec2S &cumulative_value() const;
};

///// PISM surface model implementing a PDD scheme.

TemperatureIndex::TemperatureIndex(IceGrid::ConstPtr g)
  : SurfaceModel(g) {

  m_mbscheme              = NULL;
  m_faustogreve           = NULL;
  m_sd_period             = 0;
  m_sd_ref_time           = 0.0;
  m_base_ddf.snow         = m_config->get_double("surface.pdd.factor_snow");
  m_base_ddf.ice          = m_config->get_double("surface.pdd.factor_ice");
  m_base_ddf.refreezeFrac = m_config->get_double("surface.pdd.refreeze");
  m_base_pddThresholdTemp = m_config->get_double("surface.pdd.positive_threshold_temp");
  m_base_pddStdDev        = m_config->get_double("surface.pdd.std_dev");
  m_sd_use_param          = m_config->get_boolean("surface.pdd.std_dev_use_param");
  m_sd_param_a            = m_config->get_double("surface.pdd.std_dev_param_a");
  m_sd_param_b            = m_config->get_double("surface.pdd.std_dev_param_b");


  m_randomized = options::Bool("-pdd_rand",
                               "Use a PDD implementation based on simulating a random process");
  m_randomized_repeatable = options::Bool("-pdd_rand_repeatable",
                                          "Use a PDD implementation based on simulating a"
                                          " repeatable random process");
  m_use_fausto_params = options::Bool("-pdd_fausto",
                                      "Set PDD parameters using formulas (6) and (7)"
                                      " in [Faustoetal2009]");

  options::String file("-pdd_sd_file", "Read standard deviation from file");
  m_sd_file_set = file.is_set();

  options::Integer period("-pdd_sd_period",
                          "Length of the standard deviation data period in years", 0);
  m_sd_period = period;

  options::Integer sd_ref_year("-pdd_sd_reference_year",
                               "Standard deviation data reference year", 0);

  if (m_randomized_repeatable) {
    m_mbscheme = new PDDrandMassBalance(m_config, m_sys, true);
  } else if (m_randomized) {
    m_mbscheme = new PDDrandMassBalance(m_config, m_sys, false);
  } else {
    m_mbscheme = new PDDMassBalance(m_config, m_sys);
  }

  if (m_use_fausto_params) {
    m_faustogreve = new FaustoGrevePDDObject(m_grid);
  }

  if (sd_ref_year.is_set()) {
    m_sd_ref_time = units::convert(m_sys, sd_ref_year, "years", "seconds");
  }

  if (m_sd_file_set) {
    // find out how many records there are in the file and set the
    // air_temp_sd buffer size

    unsigned int n_records = 0;
    std::string short_name = "air_temp_sd";
    unsigned int buffer_size = (unsigned int) m_config->get_double("climate_forcing.buffer_size");

    {
      PIO nc(m_grid->com, "netcdf3", file, PISM_READONLY);
      n_records = nc.inq_nrecords(short_name, "", m_grid->ctx()->unit_system());
    }

    // If -..._period is not set, make ..._n_records the minimum of the
    // buffer size and the number of available records. Otherwise try
    // to keep all available records in memory.
    if (m_sd_period == 0) {
      n_records = std::min(n_records, buffer_size);
    }

    if (n_records < 1) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Can't find '%s' in %s.",
                                    short_name.c_str(), file->c_str());
    }

    m_air_temp_sd.set_n_records(n_records);

  } else {
    // using constant standard deviation, so set buffer size to 1
    m_air_temp_sd.set_n_records(1);
  }

  m_air_temp_sd.create(m_grid, "air_temp_sd");
  m_air_temp_sd.set_attrs("climate_forcing",
                          "standard deviation of near-surface air temperature",
                          "Kelvin", "");

  m_climatic_mass_balance.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  m_climatic_mass_balance.set_attrs("diagnostic",
                                    "instantaneous surface mass balance (accumulation/ablation) rate",
                                    "kg m-2 s-1",
                                    "land_ice_surface_specific_mass_balance_flux");
  m_climatic_mass_balance.metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_climatic_mass_balance.write_in_glaciological_units = true;
  m_climatic_mass_balance.metadata().set_string("comment", "positive values correspond to ice gain");

  // diagnostic fields:

  m_accumulation_rate.create(m_grid, "saccum", WITHOUT_GHOSTS);
  m_accumulation_rate.set_attrs("diagnostic",
                                "instantaneous surface accumulation rate"
                                " (precipitation minus rain)",
                                "kg m-2 s-1",
                                "");
  m_accumulation_rate.metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_accumulation_rate.write_in_glaciological_units = true;

  m_melt_rate.create(m_grid, "smelt", WITHOUT_GHOSTS);
  m_melt_rate.set_attrs("diagnostic",
                        "instantaneous surface melt rate",
                        "kg m-2 s-1",
                        "");
  m_melt_rate.metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_melt_rate.write_in_glaciological_units = true;

  m_runoff_rate.create(m_grid, "srunoff", WITHOUT_GHOSTS);
  m_runoff_rate.set_attrs("diagnostic",
                          "instantaneous surface meltwater runoff rate",
                          "kg m-2 s-1",
                          "");
  m_runoff_rate.metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_runoff_rate.write_in_glaciological_units = true;

  {
    m_cumulative_accumulation.create(m_grid, "saccum_cumulative", WITHOUT_GHOSTS);
    m_cumulative_accumulation.set_attrs("diagnostic",
                                       "cumulative surface accumulation"
                                       " (precipitation minus rain)",
                                       "kg m-2", "");

    m_cumulative_melt.create(m_grid, "smelt_cumulative", WITHOUT_GHOSTS);
    m_cumulative_melt.set_attrs("diagnostic", "cumulative surface melt",
                               "kg m-2", "");

    m_cumulative_runoff.create(m_grid, "srunoff_cumulative", WITHOUT_GHOSTS);
    m_cumulative_runoff.set_attrs("diagnostic",
                                 "cumulative surface meltwater runoff",
                                 "kg m-2", "");
  }

  m_snow_depth.create(m_grid, "snow_depth", WITHOUT_GHOSTS);
  m_snow_depth.set_attrs("diagnostic",
                         "snow cover depth (set to zero once a year)",
                         "m", "");
  m_snow_depth.set(0.0);
}

TemperatureIndex::~TemperatureIndex() {
  delete m_mbscheme;
  delete m_faustogreve;
}

void TemperatureIndex::init_impl() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  // call the default implementation (not the interface method init())
  SurfaceModel::init_impl();

  m_log->message(2,
             "* Initializing the default temperature-index, PDD-based surface processes scheme.\n"
             "  Precipitation and 2m air temperature provided by atmosphere are inputs.\n"
             "  Surface mass balance and ice upper surface temperature are outputs.\n"
             "  See PISM User's Manual for control of degree-day factors.\n");

  m_log->message(2,
             "  Computing number of positive degree-days by: ");
  if (m_randomized) {
    m_log->message(2, "simulation of a random process.\n");
  } else if (m_randomized_repeatable) {
    m_log->message(2, "repeatable simulation of a random process.\n");
  } else {
    m_log->message(2, "an expectation integral.\n");
  }

  if (m_use_fausto_params) {
    m_log->message(2,
               "  Setting PDD parameters from [Faustoetal2009] ...\n");

    m_base_pddStdDev = 2.53;

  }

  options::String file("-pdd_sd_file", "Read standard deviation from file");
  if (file.is_set()) {
    m_log->message(2,
               "  Reading standard deviation of near-surface temperature from '%s'...\n",
               file->c_str());
    m_air_temp_sd.init(file, m_sd_period, m_sd_ref_time);
  } else {
    m_log->message(2,
               "  Option -pdd_sd_file is not set. Using a constant value.\n");
    m_air_temp_sd.init_constant(m_base_pddStdDev);
  }

  std::string input_file = process_input_options(m_grid->com).filename;

  // read snow precipitation rate from file
  m_log->message(2,
                 "    reading snow depth (ice equivalent meters) from %s ... \n",
                 input_file.c_str());
  m_snow_depth.regrid(input_file, OPTIONAL, 0.0);

  m_next_balance_year_start = compute_next_balance_year_start(m_grid->ctx()->time()->current());

  m_cumulative_accumulation.set(0.0);
  m_cumulative_melt.set(0.0);
  m_cumulative_runoff.set(0.0);
}

MaxTimestep TemperatureIndex::max_timestep_impl(double my_t) const {
  return m_atmosphere->max_timestep(my_t);
}

double TemperatureIndex::compute_next_balance_year_start(double time) {
  // compute the time corresponding to the beginning of the next balance year
  double
    balance_year_start_day = m_config->get_double("surface.pdd.balance_year_start_day"),
    one_day                = units::convert(m_sys, 1.0, "days", "seconds"),
    year_start             = m_grid->ctx()->time()->calendar_year_start(time),
    balance_year_start     = year_start + (balance_year_start_day - 1.0) * one_day;

  if (balance_year_start > time) {
    return balance_year_start;
  }
  return m_grid->ctx()->time()->increment_date(balance_year_start, 1);
}

void TemperatureIndex::update_impl(double my_t, double my_dt) {

  if ((fabs(my_t - m_t) < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = my_t;
  m_dt = my_dt;

  // update to ensure that temperature and precipitation time series
  // are correct:
  m_atmosphere->update(my_t, my_dt);

  // set up air temperature and precipitation time series
  int Nseries = m_mbscheme->get_timeseries_length(my_dt);

  const double dtseries = my_dt / Nseries;
  std::vector<double> ts(Nseries), T(Nseries), S(Nseries), P(Nseries), PDDs(Nseries);
  for (int k = 0; k < Nseries; ++k) {
    ts[k] = my_t + k * dtseries;
  }

  // update standard deviation time series
  if (m_sd_file_set == true) {
    m_air_temp_sd.update(my_t, my_dt);
    m_air_temp_sd.init_interpolation(ts);
  }

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");
  const IceModelVec2S
    *surface_altitude = NULL,
    *latitude         = NULL,
    *longitude        = NULL;

  IceModelVec::AccessList list(mask);

  if (m_faustogreve != NULL) {
    surface_altitude = m_grid->variables().get_2d_scalar("surface_altitude");
    latitude         = m_grid->variables().get_2d_scalar("latitude");
    longitude        = m_grid->variables().get_2d_scalar("longitude");

    list.add(*latitude);
    list.add(*longitude);
    list.add(*surface_altitude);
    m_faustogreve->update_temp_mj(*surface_altitude, *latitude, *longitude);
  }

  const double
    sigmalapserate = m_config->get_double("surface.pdd.std_dev_lapse_lat_rate"),
    sigmabaselat   = m_config->get_double("surface.pdd.std_dev_lapse_lat_base");

  if (sigmalapserate != 0.0) {
    latitude = m_grid->variables().get_2d_scalar("latitude");
    list.add(*latitude);
  }

  LocalMassBalance::DegreeDayFactors  ddf = m_base_ddf;

  m_atmosphere->init_timeseries(ts);

  m_atmosphere->begin_pointwise_access();
  list.add(m_air_temp_sd);
  list.add(m_climatic_mass_balance);

  list.add(m_accumulation_rate);
  list.add(m_melt_rate);
  list.add(m_runoff_rate);
  list.add(m_snow_depth);
  list.add(m_cumulative_accumulation);
  list.add(m_cumulative_melt);
  list.add(m_cumulative_runoff);

  const double ice_density = m_config->get_double("constants.ice.density");

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // the temperature time series from the AtmosphereModel and its modifiers
      m_atmosphere->temp_time_series(i, j, T);

      // the precipitation time series from AtmosphereModel and its modifiers
      m_atmosphere->precip_time_series(i, j, P);

      // convert precipitation from "kg m-2 second-1" to "m second-1" (PDDMassBalance expects
      // accumulation in m/second ice equivalent)
      for (int k = 0; k < Nseries; ++k) {
        P[k] = P[k] / ice_density;
        // kg / (m^2 * second) / (kg / m^3) = m / second
      }

      // interpolate temperature standard deviation time series
      if (m_sd_file_set == true) {
        m_air_temp_sd.interp(i, j, S);
      } else {
        for (int k = 0; k < Nseries; ++k) {
          S[k] = m_air_temp_sd(i, j);
        }
      }

      if (m_faustogreve != NULL) {
        // we have been asked to set mass balance parameters according to
        //   formula (6) in [\ref Faustoetal2009]; they overwrite ddf set above
        m_faustogreve->setDegreeDayFactors(i, j, (*surface_altitude)(i, j),
                                           (*latitude)(i, j), (*longitude)(i, j), ddf);
      }

      // apply standard deviation lapse rate on top of prescribed values
      if (sigmalapserate != 0.0) {
        for (int k = 0; k < Nseries; ++k) {
          S[k] += sigmalapserate * ((*latitude)(i,j) - sigmabaselat);
        }
        m_air_temp_sd(i, j) = S[0]; // ensure correct SD reporting
      }

      // apply standard deviation param over ice if in use
      if (m_sd_use_param && mask.icy(i,j)) {
        for (int k = 0; k < Nseries; ++k) {
          S[k] = m_sd_param_a * (T[k] - 273.15) + m_sd_param_b;
          if (S[k] < 0.0) {
            S[k] = 0.0 ;
          }
        }
        m_air_temp_sd(i, j) = S[0]; // ensure correct SD reporting
      }

      // Use temperature time series, the "positive" threshhold, and
      // the standard deviation of the daily variability to get the
      // number of positive degree days (PDDs)
      m_mbscheme->get_PDDs(&S[0], dtseries, &T[0], Nseries, &PDDs[0]);

      // Use temperature time series to remove rainfall from precipitation
      m_mbscheme->get_snow_accumulation(&P[0], // precipitation rate (input-output)
                                        &T[0], // air temperature (input)
                                        Nseries);

      // Use degree-day factors, and number of PDDs, and the snow
      // precipitation, to get surface mass balance (and diagnostics:
      // accumulation, melt, runoff)
      {
        double next_snow_depth_reset = m_next_balance_year_start;
        m_accumulation_rate(i,j)     = 0.0;
        m_melt_rate(i,j)             = 0.0;
        m_runoff_rate(i,j)           = 0.0;
        m_climatic_mass_balance(i,j) = 0.0;
        for (int k = 0; k < Nseries; ++k) {
          if (ts[k] >= next_snow_depth_reset) {
            m_snow_depth(i,j)       = 0.0;
            while (next_snow_depth_reset <= ts[k]) {
              next_snow_depth_reset = m_grid->ctx()->time()->increment_date(next_snow_depth_reset, 1);
            }
          }

          double accumulation     = P[k] * dtseries;
          m_accumulation_rate(i,j) += accumulation;

          m_mbscheme->step(ddf, PDDs[k], accumulation,
                           m_snow_depth(i,j), m_melt_rate(i,j), m_runoff_rate(i,j),
                           m_climatic_mass_balance(i,j));
        }

        // convert from [m during the current time-step] to kg m-2 s-1
        m_accumulation_rate(i,j)     *= (ice_density/m_dt);
        m_melt_rate(i,j)             *= (ice_density/m_dt);
        m_runoff_rate(i,j)           *= (ice_density/m_dt);
        m_climatic_mass_balance(i,j) *= (ice_density/m_dt);
      }

      // update cumulative quantities
      {
        m_cumulative_accumulation(i, j) += m_accumulation_rate(i, j) * m_dt;
        m_cumulative_melt(i, j)         += m_melt_rate(i, j) * m_dt;
        m_cumulative_runoff(i, j)       += m_runoff_rate(i, j) * m_dt;
      }

      if (mask.ocean(i,j)) {
        m_snow_depth(i,j) = 0.0;  // snow over the ocean does not stick
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_atmosphere->end_pointwise_access();

  m_next_balance_year_start = compute_next_balance_year_start(m_grid->ctx()->time()->current());
}

void TemperatureIndex::ice_surface_mass_flux_impl(IceModelVec2S &result) const {
  result.copy_from(m_climatic_mass_balance);
}


void TemperatureIndex::ice_surface_temperature_impl(IceModelVec2S &result) const {
  m_atmosphere->mean_annual_temp(result);
}

const IceModelVec2S& TemperatureIndex::cumulative_surface_accumulation() const {
  return m_cumulative_accumulation;
}

const IceModelVec2S& TemperatureIndex::cumulative_surface_melt() const {
  return m_cumulative_melt;
}

const IceModelVec2S& TemperatureIndex::cumulative_surface_runoff() const {
  return m_cumulative_runoff;
}

const IceModelVec2S& TemperatureIndex::surface_accumulation() const {
  return m_accumulation_rate;
}

const IceModelVec2S& TemperatureIndex::surface_melt() const {
  return m_melt_rate;
}

const IceModelVec2S& TemperatureIndex::surface_runoff() const {
  return m_runoff_rate;
}

const IceModelVec2S& TemperatureIndex::snow_depth() const {
  return m_snow_depth;
}

const IceModelVec2S& TemperatureIndex::air_temp_sd() const {
  return m_air_temp_sd;
}

void TemperatureIndex::define_model_state_impl(const PIO &output) const {
  SurfaceModel::define_model_state_impl(output);
  m_snow_depth.define(output, PISM_DOUBLE);
}

void TemperatureIndex::write_model_state_impl(const PIO &output) const {
  SurfaceModel::write_model_state_impl(output);
  m_snow_depth.write(output);
}

std::map<std::string, Diagnostic::Ptr> TemperatureIndex::diagnostics_impl() const {
  std::map<std::string, Diagnostic::Ptr> result = {
    {"saccum",          Diagnostic::Ptr(new PDD_saccum(this))},
    {"smelt",           Diagnostic::Ptr(new PDD_smelt(this))},
    {"srunoff",         Diagnostic::Ptr(new PDD_srunoff(this))},
    {"saccum_average",  Diagnostic::Ptr(new PDD_saccum_average(this))},
    {"smelt_average",   Diagnostic::Ptr(new PDD_smelt_average(this))},
    {"srunoff_average", Diagnostic::Ptr(new PDD_srunoff_average(this))},
    {"air_temp_sd",     Diagnostic::Ptr(new PDD_air_temp_sd(this))}
  };

  result = pism::combine(result, SurfaceModel::diagnostics_impl());

  return result;
}

PDD_saccum::PDD_saccum(const TemperatureIndex *m)
  : Diag<TemperatureIndex>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "saccum")};

  set_attrs("instantaneous surface accumulation rate (precipitation minus rain)", "",
            "kg m-2 s-1", "kg m-2 year-1", 0);
}

IceModelVec::Ptr PDD_saccum::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "saccum", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  result->copy_from(model->surface_accumulation());

  return result;
}

PDD_smelt::PDD_smelt(const TemperatureIndex *m)
  : Diag<TemperatureIndex>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "smelt")};

  set_attrs("instantaneous surface melt rate", "",
            "kg m-2 s-1", "kg m-2 year-1", 0);
}

IceModelVec::Ptr PDD_smelt::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "smelt", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  result->copy_from(model->surface_melt());

  return result;
}

PDD_srunoff::PDD_srunoff(const TemperatureIndex *m)
  : Diag<TemperatureIndex>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "srunoff")};

  set_attrs("instantaneous surface meltwater runoff rate", "",
            "kg m-2 s-1", "kg m-2 year-1", 0);
}

IceModelVec::Ptr PDD_srunoff::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "srunoff", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  result->copy_from(model->surface_runoff());

  return result;
}

PDD_snow_depth::PDD_snow_depth(const TemperatureIndex *m)
  : Diag<TemperatureIndex>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "snow_depth")};

  set_attrs("snow cover depth (set to zero once a year)", "",
            "m", "m", 0);
}

IceModelVec::Ptr PDD_snow_depth::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "snow_depth", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  result->copy_from(model->snow_depth());

  return result;
}

PDD_air_temp_sd::PDD_air_temp_sd(const TemperatureIndex *m)
  : Diag<TemperatureIndex>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "air_temp_sd")};

  set_attrs("standard deviation of near-surface air temperature", "",
            "Kelvin", "Kelvin", 0);
}

IceModelVec::Ptr PDD_air_temp_sd::compute_impl() {

  IceModelVec2S::Ptr result(new IceModelVec2S);
  result->create(m_grid, "air_temp_sd", WITHOUT_GHOSTS);
  result->metadata(0) = m_vars[0];

  result->copy_from(model->air_temp_sd());

  return result;
}

PDD_saccum_average::PDD_saccum_average(const TemperatureIndex *m)
  : Diag_average<TemperatureIndex>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "saccum_average")};

  set_attrs("surface accumulation averaged over reporting intervals", "",
            "kg m-2 s-1", "kg m-2 year-1", 0);
  double fill_value = units::convert(m_sys, m_fill_value,
                                     "kg year-1", "kg second-1");
  m_vars[0].set_double("_FillValue", fill_value);
}

const IceModelVec2S & PDD_saccum_average::cumulative_value() const {
  return model->cumulative_surface_accumulation();
}

PDD_smelt_average::PDD_smelt_average(const TemperatureIndex *m)
  : Diag_average<TemperatureIndex>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "smelt_average")};

  set_attrs("surface melt averaged over reporting intervals", "",
            "kg m-2 s-1", "kg m-2 year-1", 0);
  double fill_value = units::convert(m_sys, m_fill_value,
                                     "kg year-1", "kg second-1");
  m_vars[0].set_double("_FillValue", fill_value);
}

const IceModelVec2S & PDD_smelt_average::cumulative_value() const {
  return model->cumulative_surface_melt();
}

PDD_srunoff_average::PDD_srunoff_average(const TemperatureIndex *m)
  : Diag_average<TemperatureIndex>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "srunoff_average")};

  set_attrs("surface runoff averaged over reporting intervals", "",
            "kg m-2 s-1", "kg m-2 year-1", 0);
  double fill_value = units::convert(m_sys, m_fill_value,
                                     "kg year-1", "kg second-1");
  m_vars[0].set_double("_FillValue", fill_value);
}

const IceModelVec2S & PDD_srunoff_average::cumulative_value() const {
  return model->cumulative_surface_runoff();
}

} // end of namespace surface
} // end of namespace pism
