// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017 PISM Authors
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
#include <gsl/gsl_math.h>

#include "TemperatureIndex.hh"
#include "localMassBalance.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/Vars.hh"
#include "pism/util/Time.hh"
#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/Mask.hh"
#include "pism/util/io/PIO.hh"

#include "pism/util/error_handling.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/IceModelVec2CellType.hh"

namespace pism {
namespace surface {

namespace diagnostics {

/*! @brief Snow cover depth. */
class PDD_snow_depth : public Diag<TemperatureIndex>
{
public:
  PDD_snow_depth(const TemperatureIndex *m);
protected:
  IceModelVec::Ptr compute_impl() const;
};

/*! @brief Standard deviation of near-surface air temperature. */
class PDD_air_temp_sd : public Diag<TemperatureIndex>
{
public:
  PDD_air_temp_sd(const TemperatureIndex *m);
protected:
  IceModelVec::Ptr compute_impl() const;
};

PDD_snow_depth::PDD_snow_depth(const TemperatureIndex *m)
  : Diag<TemperatureIndex>(m) {

  /* set metadata: */
  m_vars = {SpatialVariableMetadata(m_sys, "snow_depth")};

  set_attrs("snow cover depth (set to zero once a year)", "",
            "m", "m", 0);
}

IceModelVec::Ptr PDD_snow_depth::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "snow_depth", WITHOUT_GHOSTS));
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

IceModelVec::Ptr PDD_air_temp_sd::compute_impl() const {

  IceModelVec2S::Ptr result(new IceModelVec2S(m_grid, "air_temp_sd", WITHOUT_GHOSTS));
  result->metadata(0) = m_vars[0];

  result->copy_from(model->air_temp_sd());

  return result;
}

} // end of namespace diagnostics

///// PISM surface model implementing a PDD scheme.

TemperatureIndex::TemperatureIndex(IceGrid::ConstPtr g)
  : SurfaceModel(g) {

  m_mbscheme                   = NULL;
  m_faustogreve                = NULL;
  m_sd_period                  = 0;
  m_base_ddf.snow              = m_config->get_double("surface.pdd.factor_snow");
  m_base_ddf.ice               = m_config->get_double("surface.pdd.factor_ice");
  m_base_ddf.refreeze_fraction = m_config->get_double("surface.pdd.refreeze");
  m_base_pddStdDev             = m_config->get_double("surface.pdd.std_dev");
  m_sd_use_param               = m_config->get_boolean("surface.pdd.std_dev_use_param");
  m_sd_param_a                 = m_config->get_double("surface.pdd.std_dev_param_a");
  m_sd_param_b                 = m_config->get_double("surface.pdd.std_dev_param_b");

  bool randomized = options::Bool("-pdd_rand",
                                  "Use a PDD implementation based on simulating a random process");
  bool randomized_repeatable = options::Bool("-pdd_rand_repeatable",
                                             "Use a PDD implementation based on simulating a"
                                             " repeatable random process");
  bool use_fausto_params = options::Bool("-pdd_fausto",
                                         "Set PDD parameters using formulas (6) and (7)"
                                         " in [Faustoetal2009]");

  std::string sd_file = m_config->get_string("surface.pdd.temperature_standard_deviation_file");

  m_sd_file_set = not sd_file.empty();

  options::Integer period("-pdd_sd_period",
                          "Length of the standard deviation data period in years", 0);
  m_sd_period = period;

  if (randomized_repeatable) {
    m_mbscheme = new PDDrandMassBalance(m_config, m_sys, PDDrandMassBalance::REPEATABLE);
  } else if (randomized) {
    m_mbscheme = new PDDrandMassBalance(m_config, m_sys, PDDrandMassBalance::NOT_REPEATABLE);
  } else {
    m_mbscheme = new PDDMassBalance(m_config, m_sys);
  }

  if (use_fausto_params) {
    m_faustogreve = new FaustoGrevePDDObject(m_grid);
    m_base_pddStdDev = 2.53;
  }

  if (m_sd_file_set) {
    // find out how many records there are in the file and set the
    // air_temp_sd buffer size

    unsigned int n_records = 0;
    std::string short_name = "air_temp_sd";
    unsigned int buffer_size = (unsigned int) m_config->get_double("climate_forcing.buffer_size");

    {
      PIO nc(m_grid->com, "netcdf3", sd_file, PISM_READONLY);
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
                                    short_name.c_str(), sd_file.c_str());
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
  m_climatic_mass_balance.metadata().set_string("comment", "positive values correspond to ice gain");

  // diagnostic fields:

  {
    m_accumulation.create(m_grid, "saccum", WITHOUT_GHOSTS);
    m_accumulation.set_attrs("diagnostic", "surface accumulation (precipitation minus rain)",
                             "kg m-2", "");

    m_melt.create(m_grid, "smelt", WITHOUT_GHOSTS);
    m_melt.set_attrs("diagnostic", "surface melt", "kg m-2", "");

    m_runoff.create(m_grid, "srunoff", WITHOUT_GHOSTS);
    m_runoff.set_attrs("diagnostic", "surface meltwater runoff",
                       "kg m-2", "");
  }

  m_snow_depth.create(m_grid, "snow_depth", WITHOUT_GHOSTS);
  m_snow_depth.set_attrs("diagnostic",
                         "snow cover depth (set to zero once a year)",
                         "m", "");
  m_snow_depth.set(0.0);

  m_firn_depth.create(m_grid, "firn_depth", WITHOUT_GHOSTS);
  m_firn_depth.set_attrs("diagnostic",
                         "firn cover depth",
                         "m", "");
  m_firn_depth.set(0.0);
}

TemperatureIndex::~TemperatureIndex() {
  delete m_mbscheme;
  delete m_faustogreve;
}

void TemperatureIndex::init_impl() {
  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  // call the default implementation (not the interface method init())
  SurfaceModel::init_impl();

  // report user's modeling choices
  {
    m_log->message(2,
                   "* Initializing the default temperature-index, PDD-based surface processes scheme.\n"
                   "  Precipitation and 2m air temperature provided by atmosphere are inputs.\n"
                   "  Surface mass balance and ice upper surface temperature are outputs.\n"
                   "  See PISM User's Manual for control of degree-day factors.\n");

    m_log->message(2,
                   "  Computing number of positive degree-days by: %s.\n",
                   m_mbscheme->method().c_str());

    if (m_faustogreve) {
      m_log->message(2,
                     "  Setting PDD parameters from [Faustoetal2009].\n");
    } else {
      m_log->message(2,
                     "  Using default PDD parameters.\n");
    }
  }

  // initialize the spatially-variable air temperature standard deviation
  {
    std::string sd_file = m_config->get_string("surface.pdd.temperature_standard_deviation_file");
    if (sd_file.empty()) {
      m_log->message(2,
                     "  Using constant standard deviation of near-surface temperature.\n");
      m_air_temp_sd.init_constant(m_base_pddStdDev);
    } else {
      m_log->message(2,
                     "  Reading standard deviation of near-surface temperature from '%s'...\n",
                     sd_file.c_str());

      options::Integer sd_ref_year("-pdd_sd_reference_year",
                                   "Standard deviation data reference year", 0);

      double sd_ref_time = units::convert(m_sys, sd_ref_year, "years", "seconds");

      m_air_temp_sd.init(sd_file, m_sd_period, sd_ref_time);
    }
  }

  // initializing the model state
  InputOptions input = process_input_options(m_grid->com);

  std::string firn_file = m_config->get_string("surface.pdd.firn_depth_file");

  if (input.type == INIT_RESTART) {
    if (not firn_file.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "surface.pdd.firn_depth_file is not allowed when"
                                    " re-starting from a PISM output file.");
    }

    m_firn_depth.read(input.filename, input.record);
    m_snow_depth.read(input.filename, input.record);
  } else if (input.type == INIT_BOOTSTRAP) {

    m_snow_depth.regrid(input.filename, OPTIONAL, 0.0);

    if (firn_file.empty()) {
      m_firn_depth.regrid(input.filename, OPTIONAL, 0.0);
    } else {
      m_firn_depth.regrid(firn_file, CRITICAL);
    }
  } else {

    m_snow_depth.set(0.0);

    if (firn_file.empty()) {
      m_firn_depth.set(0.0);
    } else {
      m_firn_depth.regrid(firn_file, CRITICAL);
    }
  }

  {
    regrid("PDD surface model", m_snow_depth);
    regrid("PDD surface model", m_firn_depth);
  }

  // finish up
  {
    m_next_balance_year_start = compute_next_balance_year_start(m_grid->ctx()->time()->current());

    m_accumulation.set(0.0);
    m_melt.set(0.0);
    m_runoff.set(0.0);
  }
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

void TemperatureIndex::update_impl(double t, double dt) {

  if ((fabs(t - m_t) < 1e-12) &&
      (fabs(dt - m_dt) < 1e-12)) {
    return;
  }

  // make a copy of the pointer to convince clang static analyzer that its value does not
  // change during the call
  FaustoGrevePDDObject *fausto_greve = m_faustogreve;

  m_t  = t;
  m_dt = dt;

  // update to ensure that temperature and precipitation time series are correct:
  m_atmosphere->update(t, dt);

  // set up air temperature and precipitation time series
  int N = m_mbscheme->get_timeseries_length(dt);

  const double dtseries = dt / N;
  std::vector<double> ts(N), T(N), S(N), P(N), PDDs(N);
  for (int k = 0; k < N; ++k) {
    ts[k] = t + k * dtseries;
  }

  // update standard deviation time series
  if (m_sd_file_set) {
    m_air_temp_sd.update(t, dt);
    m_air_temp_sd.init_interpolation(ts);
  }

  const IceModelVec2CellType &mask = *m_grid->variables().get_2d_cell_type("mask");

  IceModelVec::AccessList list{&mask, &m_air_temp_sd, &m_climatic_mass_balance,
      &m_firn_depth, &m_snow_depth, &m_accumulation, &m_melt, &m_runoff};

  const double
    sigmalapserate = m_config->get_double("surface.pdd.std_dev_lapse_lat_rate"),
    sigmabaselat   = m_config->get_double("surface.pdd.std_dev_lapse_lat_base");

  const IceModelVec2S *latitude = NULL;
  if (fausto_greve or sigmalapserate != 0.0) {
    latitude = m_grid->variables().get_2d_scalar("latitude");

    list.add({latitude});
  }

  if (fausto_greve) {
    const IceModelVec2S
      *longitude        = m_grid->variables().get_2d_scalar("longitude"),
      *surface_altitude = m_grid->variables().get_2d_scalar("surface_altitude");

    fausto_greve->update_temp_mj(*surface_altitude, *latitude, *longitude);
  }

  LocalMassBalance::DegreeDayFactors ddf = m_base_ddf;

  m_atmosphere->init_timeseries(ts);

  m_atmosphere->begin_pointwise_access();

  const double ice_density = m_config->get_double("constants.ice.density");

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      // reset total accumulation, melt, and runoff, and SMB
      {
        m_accumulation(i, j)          = 0.0;
        m_melt(i, j)                  = 0.0;
        m_runoff(i, j)                = 0.0;
        m_climatic_mass_balance(i, j) = 0.0;
      }

      // the temperature time series from the AtmosphereModel and its modifiers
      m_atmosphere->temp_time_series(i, j, T);

      // the precipitation time series from AtmosphereModel and its modifiers
      m_atmosphere->precip_time_series(i, j, P);


      // convert precipitation from "kg m-2 second-1" to "m second-1" (PDDMassBalance expects
      // accumulation in m/second ice equivalent)
      for (int k = 0; k < N; ++k) {
        P[k] = P[k] / ice_density;
        // kg / (m^2 * second) / (kg / m^3) = m / second
      }

      // interpolate temperature standard deviation time series
      if (m_sd_file_set) {
        m_air_temp_sd.interp(i, j, S);
      } else {
        for (int k = 0; k < N; ++k) {
          S[k] = m_air_temp_sd(i, j);
        }
      }

      if (fausto_greve) {
        // we have been asked to set mass balance parameters according to
        //   formula (6) in [\ref Faustoetal2009]; they overwrite ddf set above
        ddf = fausto_greve->degree_day_factors(i, j, (*latitude)(i, j));
      }

      // apply standard deviation lapse rate on top of prescribed values
      if (sigmalapserate != 0.0) {
        for (int k = 0; k < N; ++k) {
          S[k] += sigmalapserate * ((*latitude)(i,j) - sigmabaselat);
        }
        m_air_temp_sd(i, j) = S[0]; // ensure correct SD reporting
      }

      // apply standard deviation param over ice if in use
      if (m_sd_use_param and mask.icy(i, j)) {
        for (int k = 0; k < N; ++k) {
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
      m_mbscheme->get_PDDs(dtseries, S, T, // inputs
                           PDDs);          // output

      // Use temperature time series to remove rainfall from precipitation
      m_mbscheme->get_snow_accumulation(T,  // air temperature (input)
                                        P); // precipitation rate (input-output)

      // Use degree-day factors, the number of PDDs, and the snow precipitation to get surface mass
      // balance (and diagnostics: accumulation, melt, runoff)
      {
        double next_snow_depth_reset = m_next_balance_year_start;

        for (int k = 0; k < N; ++k) {
          if (ts[k] >= next_snow_depth_reset) {
            m_snow_depth(i,j)       = 0.0;
            while (next_snow_depth_reset <= ts[k]) {
              next_snow_depth_reset = m_grid->ctx()->time()->increment_date(next_snow_depth_reset, 1);
            }
          }

          const double accumulation = P[k] * dtseries;

          LocalMassBalance::Changes changes;
          changes = m_mbscheme->step(ddf, PDDs[k], m_firn_depth(i, j), m_snow_depth(i, j), accumulation);

          // update firn depth
          m_firn_depth(i, j) += changes.firn_depth;
          // update snow depth
          m_snow_depth(i, j) += changes.snow_depth;

          // update total accumulation, melt, and runoff, converting from "meters, ice equivalent"
          // to "kg / meter^2"
          {
            m_accumulation(i, j) += accumulation * ice_density;
            m_melt(i, j)         += changes.melt * ice_density;
            m_runoff(i, j)       += changes.runoff * ice_density;
          }

          // m_climatic_mass_balance (unlike m_accumulation, m_melt, and m_runoff), is a rate.
          // m * (kg / m^3) / second = kg / m^2 / second
          m_climatic_mass_balance(i, j) += changes.smb * ice_density / m_dt;
        } // end of the time-stepping loop
      }

      if (mask.ocean(i,j)) {
        m_firn_depth(i,j) = 0.0;  // no firn over the ocean
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

void TemperatureIndex::mass_flux_impl(IceModelVec2S &result) const {
  result.copy_from(m_climatic_mass_balance);
}

void TemperatureIndex::temperature_impl(IceModelVec2S &result) const {
  m_atmosphere->mean_annual_temp(result);
}

const IceModelVec2S& TemperatureIndex::accumulation() const {
  return m_accumulation;
}

const IceModelVec2S& TemperatureIndex::melt() const {
  return m_melt;
}

const IceModelVec2S& TemperatureIndex::runoff() const {
  return m_runoff;
}

const IceModelVec2S& TemperatureIndex::firn_depth() const {
  return m_firn_depth;
}

const IceModelVec2S& TemperatureIndex::snow_depth() const {
  return m_snow_depth;
}

const IceModelVec2S& TemperatureIndex::air_temp_sd() const {
  return m_air_temp_sd;
}

void TemperatureIndex::define_model_state_impl(const PIO &output) const {
  SurfaceModel::define_model_state_impl(output);
  m_firn_depth.define(output, PISM_DOUBLE);
  m_snow_depth.define(output, PISM_DOUBLE);
}

void TemperatureIndex::write_model_state_impl(const PIO &output) const {
  SurfaceModel::write_model_state_impl(output);
  m_firn_depth.write(output);
  m_snow_depth.write(output);
}

namespace diagnostics {

/*! @brief Report surface melt, averaged over the reporting interval */
class SurfaceMelt : public DiagAverageRate<TemperatureIndex>
{
public:
  SurfaceMelt(const TemperatureIndex *m)
    : DiagAverageRate<TemperatureIndex>(m, "smelt", TOTAL_CHANGE) {

    m_vars = {SpatialVariableMetadata(m_sys, "smelt")};
    m_accumulator.metadata().set_string("units", "kg m-2");

    set_attrs("surface melt, averaged over the reporting interval", "",
              "kg m-2 s-1", "kg m-2 year-1", 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value,
                                       m_vars[0].get_string("glaciological_units"),
                                       m_vars[0].get_string("units"));
    m_vars[0].set_double("_FillValue", fill_value);
  }

protected:
  const IceModelVec2S& model_input() {
    return model->melt();
  }
};

/*! @brief Report surface runoff, averaged over the reporting interval */
class SurfaceRunoff : public DiagAverageRate<TemperatureIndex>
{
public:
  SurfaceRunoff(const TemperatureIndex *m)
    : DiagAverageRate<TemperatureIndex>(m, "srunoff", TOTAL_CHANGE) {

    m_vars = {SpatialVariableMetadata(m_sys, "srunoff")};
    m_accumulator.metadata().set_string("units", "kg m-2");

    set_attrs("surface runoff, averaged over the reporting interval", "",
              "kg m-2 s-1", "kg m-2 year-1", 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value,
                                       m_vars[0].get_string("glaciological_units"),
                                       m_vars[0].get_string("units"));
    m_vars[0].set_double("_FillValue", fill_value);
  }

protected:
  const IceModelVec2S& model_input() {
    return model->runoff();
  }
};

/*! @brief Report accumulation (precipitation minus rain), averaged over the reporting interval */
class Accumulation : public DiagAverageRate<TemperatureIndex>
{
public:
  Accumulation(const TemperatureIndex *m)
    : DiagAverageRate<TemperatureIndex>(m, "saccum", TOTAL_CHANGE) {

    m_vars = {SpatialVariableMetadata(m_sys, "saccum")};
    m_accumulator.metadata().set_string("units", "kg m-2");

    set_attrs("accumulation (precipitation minus rain), averaged over the reporting interval", "",
              "kg m-2 s-1", "kg m-2 year-1", 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value,
                                       m_vars[0].get_string("glaciological_units"),
                                       m_vars[0].get_string("units"));
    m_vars[0].set_double("_FillValue", fill_value);
  }

protected:
  const IceModelVec2S& model_input() {
    return model->accumulation();
  }
};

} // end of namespace diagnostics

std::map<std::string, Diagnostic::Ptr> TemperatureIndex::diagnostics_impl() const {
  using namespace diagnostics;

  std::map<std::string, Diagnostic::Ptr> result = {
    {"saccum",      Diagnostic::Ptr(new Accumulation(this))},
    {"smelt",       Diagnostic::Ptr(new SurfaceMelt(this))},
    {"srunoff",     Diagnostic::Ptr(new SurfaceRunoff(this))},
    {"air_temp_sd", Diagnostic::Ptr(new PDD_air_temp_sd(this))},
    {"snow_depth",  Diagnostic::Ptr(new PDD_snow_depth(this))},
    {"firn_depth",  Diagnostic::wrap(m_firn_depth)},
  };

  result = pism::combine(result, SurfaceModel::diagnostics_impl());

  return result;
}

} // end of namespace surface
} // end of namespace pism
