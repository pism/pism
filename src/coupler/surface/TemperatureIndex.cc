// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace surface {

///// PISM surface model implementing a PDD scheme.

TemperatureIndex::TemperatureIndex(IceGrid::ConstPtr g,
                                   std::shared_ptr<atmosphere::AtmosphereModel> input)
  : SurfaceModel(g, input) {

  m_sd_period                  = 0;
  m_base_ddf.snow              = m_config->get_double("surface.pdd.factor_snow");
  m_base_ddf.ice               = m_config->get_double("surface.pdd.factor_ice");
  m_base_ddf.refreeze_fraction = m_config->get_double("surface.pdd.refreeze");
  m_base_pddStdDev             = m_config->get_double("surface.pdd.std_dev");
  m_sd_use_param               = m_config->get_boolean("surface.pdd.std_dev_use_param");
  m_sd_param_a                 = m_config->get_double("surface.pdd.std_dev_param_a");
  m_sd_param_b                 = m_config->get_double("surface.pdd.std_dev_param_b");

  bool use_fausto_params = options::Bool("-pdd_fausto",
                                         "Set PDD parameters using formulas (6) and (7)"
                                         " in [Faustoetal2009]");

  options::Integer period("-pdd_sd_period",
                          "Length of the standard deviation data period in years", 0);
  m_sd_period = period;

  std::string method = m_config->get_string("surface.pdd.method");

  if (method == "repeatable_random_process") {
    m_mbscheme.reset(new PDDrandMassBalance(m_config, m_sys, PDDrandMassBalance::REPEATABLE));
  } else if (method == "random_process") {
    m_mbscheme.reset(new PDDrandMassBalance(m_config, m_sys, PDDrandMassBalance::NOT_REPEATABLE));
  } else {
    m_mbscheme.reset(new PDDMassBalance(m_config, m_sys));
  }

  if (use_fausto_params) {
    m_faustogreve.reset(new FaustoGrevePDDObject(m_grid));
    m_base_pddStdDev = 2.53;
  }

  std::string sd_file = m_config->get_string("surface.pdd.temperature_standard_deviation_file");

  if (not sd_file.empty()) {
    int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");
    int max_buffer_size = (unsigned int) m_config->get_double("climate_forcing.buffer_size");

    PIO file(m_grid->com, "netcdf3", sd_file, PISM_READONLY);
    m_air_temp_sd = IceModelVec2T::ForcingField(m_grid, file,
                                                "air_temp_sd", "",
                                                max_buffer_size,
                                                evaluations_per_year,
                                                m_sd_period > 0);
    m_sd_file_set = true;
  } else {
    m_air_temp_sd.reset(new IceModelVec2T(m_grid, "air_temp_sd", 1, 1));
    m_sd_file_set = false;
  }

  m_air_temp_sd->set_attrs("climate_forcing",
                           "standard deviation of near-surface air temperature",
                           "Kelvin", "");

  m_mass_flux.create(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS);
  m_mass_flux.set_attrs("diagnostic",
                                    "instantaneous surface mass balance (accumulation/ablation) rate",
                                    "kg m-2 s-1",
                                    "land_ice_surface_specific_mass_balance_flux");
  m_mass_flux.metadata().set_string("glaciological_units", "kg m-2 year-1");
  m_mass_flux.metadata().set_string("comment", "positive values correspond to ice gain");

  // diagnostic fields:

  {
    m_accumulation.create(m_grid, "surface_accumulation_flux", WITHOUT_GHOSTS);
    m_accumulation.set_attrs("diagnostic", "surface accumulation (precipitation minus rain)",
                             "kg m-2", "");

    m_melt.create(m_grid, "surface_melt_flux", WITHOUT_GHOSTS);
    m_melt.set_attrs("diagnostic", "surface melt", "kg m-2", "");

    m_runoff.create(m_grid, "surface_runoff_flux", WITHOUT_GHOSTS);
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
  m_firn_depth.metadata().set_double("valid_min", 0.0);
  m_firn_depth.set(0.0);

  m_temperature = allocate_temperature(g);
}

TemperatureIndex::~TemperatureIndex() {
  // empty
}

void TemperatureIndex::init_impl(const Geometry &geometry) {

  // call the default implementation (not the interface method init())
  SurfaceModel::init_impl(geometry);

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
      m_air_temp_sd->init_constant(m_base_pddStdDev);
    } else {
      m_log->message(2,
                     "  Reading standard deviation of near-surface temperature from '%s'...\n",
                     sd_file.c_str());

      options::Integer sd_ref_year("-pdd_sd_reference_year",
                                   "Standard deviation data reference year", 0);

      double sd_ref_time = units::convert(m_sys, sd_ref_year, "years", "seconds");

      m_air_temp_sd->init(sd_file, m_sd_period, sd_ref_time);
    }
  }

  // initializing the model state
  InputOptions input = process_input_options(m_grid->com, m_config);

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

void TemperatureIndex::update_impl(const Geometry &geometry, double t, double dt) {

  // make a copy of the pointer to convince clang static analyzer that its value does not
  // change during the call
  FaustoGrevePDDObject *fausto_greve = m_faustogreve.get();

  // update to ensure that temperature and precipitation time series are correct:
  m_atmosphere->update(geometry, t, dt);

  m_temperature->copy_from(m_atmosphere->mean_annual_temp());

  // set up air temperature and precipitation time series
  int N = m_mbscheme->get_timeseries_length(dt);

  const double dtseries = dt / N;
  std::vector<double> ts(N), T(N), S(N), P(N), PDDs(N);
  for (int k = 0; k < N; ++k) {
    ts[k] = t + k * dtseries;
  }

  // update standard deviation time series
  if (m_sd_file_set) {
    m_air_temp_sd->update(t, dt);
    m_air_temp_sd->init_interpolation(ts);
  }

  const IceModelVec2CellType &mask = geometry.cell_type;
  const IceModelVec2S        &H    = geometry.ice_thickness;

  IceModelVec::AccessList list{&mask, &H, m_air_temp_sd.get(), &m_mass_flux,
      &m_firn_depth, &m_snow_depth, &m_accumulation, &m_melt, &m_runoff};

  const double
    sigmalapserate = m_config->get_double("surface.pdd.std_dev_lapse_lat_rate"),
    sigmabaselat   = m_config->get_double("surface.pdd.std_dev_lapse_lat_base");

  const IceModelVec2S *latitude = nullptr;
  if (fausto_greve or sigmalapserate != 0.0) {
    latitude = &geometry.latitude;

    list.add(*latitude);
  }

  if (fausto_greve) {
    const IceModelVec2S
      *longitude        = &geometry.latitude,
      *surface_altitude = &geometry.ice_surface_elevation;

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

      // the temperature time series from the AtmosphereModel and its modifiers
      m_atmosphere->temp_time_series(i, j, T);

      if (mask.ice_free_ocean(i, j)) {
        // ignore precipitation over ice-free ocean
        for (int k = 0; k < N; ++k) {
          P[k] = 0.0;
        }
      } else {
        // elsewhere, get precipitation from the atmosphere model
        m_atmosphere->precip_time_series(i, j, P);
      }

      // convert precipitation from "kg m-2 second-1" to "m second-1" (PDDMassBalance expects
      // accumulation in m/second ice equivalent)
      for (int k = 0; k < N; ++k) {
        P[k] = P[k] / ice_density;
        // kg / (m^2 * second) / (kg / m^3) = m / second
      }

      // interpolate temperature standard deviation time series
      if (m_sd_file_set) {
        m_air_temp_sd->interp(i, j, S);
      } else {
        double tmp = (*m_air_temp_sd)(i, j);
        for (int k = 0; k < N; ++k) {
          S[k] = tmp;
        }
      }

      if (fausto_greve) {
        // we have been asked to set mass balance parameters according to
        //   formula (6) in [\ref Faustoetal2009]; they overwrite ddf set above
        ddf = fausto_greve->degree_day_factors(i, j, (*latitude)(i, j));
      }

      // apply standard deviation lapse rate on top of prescribed values
      if (sigmalapserate != 0.0) {
        double lat = (*latitude)(i, j);
        for (int k = 0; k < N; ++k) {
          S[k] += sigmalapserate * (lat - sigmabaselat);
        }
        (*m_air_temp_sd)(i, j) = S[0]; // ensure correct SD reporting
      }

      // apply standard deviation param over ice if in use
      if (m_sd_use_param and mask.icy(i, j)) {
        for (int k = 0; k < N; ++k) {
          S[k] = m_sd_param_a * (T[k] - 273.15) + m_sd_param_b;
          if (S[k] < 0.0) {
            S[k] = 0.0 ;
          }
        }
        (*m_air_temp_sd)(i, j) = S[0]; // ensure correct SD reporting
      }

      // Use temperature time series, the "positive" threshhold, and
      // the standard deviation of the daily variability to get the
      // number of positive degree days (PDDs)
      if (mask.ice_free_ocean(i, j)) {
        for (int k = 0; k < N; ++k) {
          PDDs[k] = 0.0;
        }
      } else {
        m_mbscheme->get_PDDs(dtseries, S, T, // inputs
                             PDDs);          // output
      }

      // Use temperature time series to remove rainfall from precipitation
      m_mbscheme->get_snow_accumulation(T,  // air temperature (input)
                                        P); // precipitation rate (input-output)

      // Use degree-day factors, the number of PDDs, and the snow precipitation to get surface mass
      // balance (and diagnostics: accumulation, melt, runoff)
      {
        double next_snow_depth_reset = m_next_balance_year_start;

        // make copies of firn and snow depth values at this point to avoid accessing 2D
        // fields in the inner loop
        double
          ice  = H(i, j),
          firn = m_firn_depth(i, j),
          snow = m_snow_depth(i, j);

        // accumulation, melt, runoff over this time-step
        double
          A   = 0.0,
          M   = 0.0,
          R   = 0.0,
          SMB = 0.0;

        for (int k = 0; k < N; ++k) {
          if (ts[k] >= next_snow_depth_reset) {
            snow = 0.0;
            while (next_snow_depth_reset <= ts[k]) {
              next_snow_depth_reset = m_grid->ctx()->time()->increment_date(next_snow_depth_reset, 1);
            }
          }

          const double accumulation = P[k] * dtseries;

          LocalMassBalance::Changes changes;
          changes = m_mbscheme->step(ddf, PDDs[k],
                                     ice, firn, snow, accumulation);

          // update ice thickness
          ice += changes.smb;
          assert(ice >= 0);

          // update firn depth
          firn += changes.firn_depth;
          assert(firn >= 0);

          // update snow depth
          snow += changes.snow_depth;
          assert(snow >= 0);

          // update total accumulation, melt, and runoff
          {
            A   += accumulation;
            M   += changes.melt;
            R   += changes.runoff;
            SMB += changes.smb;
          }
        } // end of the time-stepping loop

        // set firn and snow depths
        m_firn_depth(i, j) = firn;
        m_snow_depth(i, j) = snow;

        // set total accumulation, melt, and runoff, and SMB at this point, converting
        // from "meters, ice equivalent" to "kg / m^2"
        {
          m_accumulation(i, j)          = A * ice_density;
          m_melt(i, j)                  = M * ice_density;
          m_runoff(i, j)                = R * ice_density;
          // m_mass_flux (unlike m_accumulation, m_melt, and m_runoff), is a
          // rate. m * (kg / m^3) / second = kg / m^2 / second
          m_mass_flux(i, j) = SMB * ice_density / dt;
        }
      }

      if (mask.ice_free_ocean(i, j)) {
        m_firn_depth(i, j) = 0.0;  // no firn in the ocean
        m_snow_depth(i, j) = 0.0;  // snow over the ocean does not stick
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_atmosphere->end_pointwise_access();

  m_next_balance_year_start = compute_next_balance_year_start(m_grid->ctx()->time()->current());
}

const IceModelVec2S &TemperatureIndex::mass_flux_impl() const {
  return m_mass_flux;
}

const IceModelVec2S &TemperatureIndex::temperature_impl() const {
  return *m_temperature;
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
  return *m_air_temp_sd;
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

enum AmountKind {AMOUNT, MASS};

/*! @brief Report surface melt, averaged over the reporting interval */
class SurfaceMelt : public DiagAverageRate<TemperatureIndex>
{
public:
  SurfaceMelt(const TemperatureIndex *m, AmountKind kind)
    : DiagAverageRate<TemperatureIndex>(m,
                                        kind == AMOUNT
                                        ? "surface_melt_flux"
                                        : "surface_melt_rate",
                                        TOTAL_CHANGE),
    m_kind(kind) {

    std::string
      name              = "surface_melt_flux",
      long_name         = "surface melt, averaged over the reporting interval",
      standard_name     = "surface_snow_and_ice_melt_flux",
      accumulator_units = "kg m-2",
      internal_units    = "kg m-2 second-1",
      external_units    = "kg m-2 year-1";
    if (kind == MASS) {
      name              = "surface_melt_rate";
      standard_name     = "";
      accumulator_units = "kg",
      internal_units    = "kg second-1";
      external_units    = "Gt year-1" ;

      m_melt_mass.create(m_grid, "melt_mass", WITHOUT_GHOSTS);
    }

    m_vars = {SpatialVariableMetadata(m_sys, name)};
    m_accumulator.metadata().set_string("units", accumulator_units);

    set_attrs(long_name, standard_name, internal_units, external_units, 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value, external_units, internal_units);
    m_vars[0].set_double("_FillValue", fill_value);
  }

protected:
  const IceModelVec2S& model_input() {
    const IceModelVec2S &melt_amount = model->melt();

    if (m_kind == MASS) {
      double cell_area = m_grid->cell_area();

      IceModelVec::AccessList list{&m_melt_mass, &melt_amount};

      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();
        m_melt_mass(i, j) = melt_amount(i, j) * cell_area;
      }
      return m_melt_mass;
    } else {
      return melt_amount;
    }
  }
private:
  IceModelVec2S m_melt_mass;
  AmountKind m_kind;
};

/*! @brief Report surface runoff, averaged over the reporting interval */
class SurfaceRunoff : public DiagAverageRate<TemperatureIndex>
{
public:
  SurfaceRunoff(const TemperatureIndex *m, AmountKind kind)
    : DiagAverageRate<TemperatureIndex>(m,
                                        kind == AMOUNT
                                        ? "surface_runoff_flux"
                                        : "surface_runoff_rate",
                                        TOTAL_CHANGE),
    m_kind(kind) {

    std::string
      name              = "surface_runoff_flux",
      long_name         = "surface runoff, averaged over the reporting interval",
      standard_name     = "surface_runoff_flux",
      accumulator_units = "kg m-2",
      internal_units    = "kg m-2 second-1",
      external_units    = "kg m-2 year-1";
    if (kind == MASS) {
      name              = "surface_runoff_rate";
      standard_name     = "",
      accumulator_units = "kg",
      internal_units    = "kg second-1";
      external_units    = "Gt year-1" ;

      m_runoff_mass.create(m_grid, "runoff_mass", WITHOUT_GHOSTS);
    }

    m_vars = {SpatialVariableMetadata(m_sys, name)};
    m_accumulator.metadata().set_string("units", accumulator_units);

    set_attrs(long_name, standard_name, internal_units, external_units, 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value, external_units, internal_units);
    m_vars[0].set_double("_FillValue", fill_value);
  }

protected:
  const IceModelVec2S& model_input() {
    const IceModelVec2S &runoff_amount = model->runoff();

    if (m_kind == MASS) {
      double cell_area = m_grid->cell_area();

      IceModelVec::AccessList list{&m_runoff_mass, &runoff_amount};

      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();
        m_runoff_mass(i, j) = runoff_amount(i, j) * cell_area;
      }
      return m_runoff_mass;
    } else {
      return runoff_amount;
    }
  }
private:
  AmountKind m_kind;
  IceModelVec2S m_runoff_mass;
};

/*! @brief Report accumulation (precipitation minus rain), averaged over the reporting interval */
class Accumulation : public DiagAverageRate<TemperatureIndex>
{
public:
  Accumulation(const TemperatureIndex *m, AmountKind kind)
    : DiagAverageRate<TemperatureIndex>(m,
                                        kind == AMOUNT
                                        ? "surface_accumulation_flux"
                                        : "surface_accumulation_rate",
                                        TOTAL_CHANGE),
    m_kind(kind) {

    // possible standard name: surface_accumulation_flux
    std::string
      name              = "surface_accumulation_flux",
      long_name         = "accumulation (precipitation minus rain), averaged over the reporting interval",
      accumulator_units = "kg m-2",
      internal_units    = "kg m-2 second-1",
      external_units    = "kg m-2 year-1";
    if (kind == MASS) {
      name              = "surface_accumulation_rate";
      accumulator_units = "kg",
      internal_units    = "kg second-1";
      external_units    = "Gt year-1" ;

      m_accumulation_mass.create(m_grid, "accumulation_mass", WITHOUT_GHOSTS);
    }


    m_vars = {SpatialVariableMetadata(m_sys, name)};
    m_accumulator.metadata().set_string("units", accumulator_units);

    set_attrs(long_name, "", internal_units, external_units, 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value, external_units, internal_units);
    m_vars[0].set_double("_FillValue", fill_value);
  }

protected:
  const IceModelVec2S& model_input() {
    const IceModelVec2S &accumulation_amount = model->accumulation();

    if (m_kind == MASS) {
      double cell_area = m_grid->cell_area();

      IceModelVec::AccessList list{&m_accumulation_mass, &accumulation_amount};

      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();
        m_accumulation_mass(i, j) = accumulation_amount(i, j) * cell_area;
      }
      return m_accumulation_mass;
    } else {
      return accumulation_amount;
    }
  }
private:
  AmountKind m_kind;
  IceModelVec2S m_accumulation_mass;
};

/*!
 * Integrate a field over the computational domain.
 *
 * If the input has units kg/m^2, the output will be in kg.
 */
static double integrate(const IceModelVec2S &input) {
  IceGrid::ConstPtr grid = input.grid();

  double cell_area = grid->cell_area();

  IceModelVec::AccessList list{&input};

  double result = 0.0;

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    result += input(i, j) * cell_area;
  }

  return GlobalSum(grid->com, result);
}


//! \brief Reports the total accumulation rate.
class TotalSurfaceAccumulation : public TSDiag<TSFluxDiagnostic, TemperatureIndex>
{
public:
  TotalSurfaceAccumulation(const TemperatureIndex *m)
    : TSDiag<TSFluxDiagnostic, TemperatureIndex>(m, "surface_accumulation_rate") {

    m_ts.variable().set_string("units", "kg s-1");
    m_ts.variable().set_string("glaciological_units", "kg year-1");
    m_ts.variable().set_string("long_name", "surface accumulation rate (PDD model)");
  }

  double compute() {
    return integrate(model->accumulation());
  }
};


//! \brief Reports the total melt rate.
class TotalSurfaceMelt : public TSDiag<TSFluxDiagnostic, TemperatureIndex>
{
public:
  TotalSurfaceMelt(const TemperatureIndex *m)
    : TSDiag<TSFluxDiagnostic, TemperatureIndex>(m, "surface_melt_rate") {

    m_ts.variable().set_string("units", "kg s-1");
    m_ts.variable().set_string("glaciological_units", "kg year-1");
    m_ts.variable().set_string("long_name", "surface melt rate (PDD model)");
  }

  double compute() {
    return integrate(model->melt());
  }
};


//! \brief Reports the total top surface ice flux.
class TotalSurfaceRunoff : public TSDiag<TSFluxDiagnostic, TemperatureIndex>
{
public:
  TotalSurfaceRunoff(const TemperatureIndex *m)
    : TSDiag<TSFluxDiagnostic, TemperatureIndex>(m, "surface_runoff_rate") {

    m_ts.variable().set_string("units", "kg s-1");
    m_ts.variable().set_string("glaciological_units", "kg year-1");
    m_ts.variable().set_string("long_name", "surface runoff rate (PDD model)");
  }

  double compute() {
    return integrate(model->runoff());
  }
};

} // end of namespace diagnostics

DiagnosticList TemperatureIndex::diagnostics_impl() const {
  using namespace diagnostics;

  DiagnosticList result = {
    {"surface_accumulation_flux", Diagnostic::Ptr(new Accumulation(this, AMOUNT))},
    {"surface_accumulation_rate", Diagnostic::Ptr(new Accumulation(this, MASS))},
    {"surface_melt_flux",         Diagnostic::Ptr(new SurfaceMelt(this, AMOUNT))},
    {"surface_melt_rate",         Diagnostic::Ptr(new SurfaceMelt(this, MASS))},
    {"surface_runoff_flux",       Diagnostic::Ptr(new SurfaceRunoff(this, AMOUNT))},
    {"surface_runoff_rate",       Diagnostic::Ptr(new SurfaceRunoff(this, MASS))},
    {"air_temp_sd",               Diagnostic::wrap(*m_air_temp_sd)},
    {"snow_depth",                Diagnostic::wrap(m_snow_depth)},
    {"firn_depth",                Diagnostic::wrap(m_firn_depth)},
  };

  result = pism::combine(result, SurfaceModel::diagnostics_impl());

  return result;
}

TSDiagnosticList TemperatureIndex::ts_diagnostics_impl() const {
  using namespace diagnostics;

  TSDiagnosticList result = {
    {"surface_accumulation_rate", TSDiagnostic::Ptr(new TotalSurfaceAccumulation(this))},
    {"surface_melt_rate",         TSDiagnostic::Ptr(new TotalSurfaceMelt(this))},
    {"surface_runoff_rate",       TSDiagnostic::Ptr(new TotalSurfaceRunoff(this))},
  };

  result = pism::combine(result, SurfaceModel::ts_diagnostics_impl());

  return result;
}

} // end of namespace surface
} // end of namespace pism
