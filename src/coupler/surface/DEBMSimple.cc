// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2022, 2023 PISM Authors
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

#include <algorithm> // std::min

#include "pism/coupler/surface/DEBMSimple.hh"

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Time.hh"
#include "pism/util/Vars.hh"
#include "pism/util/array/CellType.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/interpolation.hh"
#include "pism/util/array/Forcing.hh"

namespace pism {
namespace surface {

///// PISM surface model implementing a dEBM-Simple scheme.

DEBMSimple::DEBMSimple(std::shared_ptr<const Grid> g, std::shared_ptr<atmosphere::AtmosphereModel> input)
  : SurfaceModel(g, std::move(input)),
      m_model(*g->ctx()),
      m_mass_flux(m_grid, "climatic_mass_balance"),
      m_snow_depth(m_grid, "snow_depth"),
      m_temperature_driven_melt(m_grid, "debm_temperature_driven_melt_flux"),
      m_insolation_driven_melt(m_grid, "debm_insolation_driven_melt_flux"),
      m_offset_melt(m_grid, "debm_offset_melt_flux"),
      m_surface_albedo(m_grid, "surface_albedo"),
      m_transmissivity(m_grid, "atmosphere_transmissivity") {

  m_sd_use_param = m_config->get_flag("surface.debm_simple.std_dev.param.enabled");
  m_sd_param_a   = m_config->get_number("surface.debm_simple.std_dev.param.a");
  m_sd_param_b   = m_config->get_number("surface.debm_simple.std_dev.param.b");

  m_precip_as_snow = m_config->get_flag("surface.debm_simple.interpret_precip_as_snow");
  m_Tmax           = m_config->get_number("surface.debm_simple.air_temp_all_precip_as_rain");
  m_Tmin           = m_config->get_number("surface.debm_simple.air_temp_all_precip_as_snow");

  if (m_Tmax <= m_Tmin) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "surface.debm_simple.air_temp_all_precip_as_rain has to exceed "
                                  "surface.debm_simple.air_temp_all_precip_as_snow");
  }

  // note: does not need to be calendar-aware
  m_year_length = units::convert(g->ctx()->unit_system(), 1.0, "years", "seconds");

  m_n_per_year = static_cast<unsigned int>(m_config->get_number("surface.debm_simple.max_evals_per_year"));

  auto albedo_input = m_config->get_string("surface.debm_simple.albedo_input.file");
  if (not albedo_input.empty()) {

    m_log->message(2, " Surface albedo is read in from %s...", albedo_input.c_str());

    File file(m_grid->com, albedo_input, io::PISM_GUESS, io::PISM_READONLY);

    int buffer_size = static_cast<int>(m_config->get_number("input.forcing.buffer_size"));
    bool periodic = m_config->get_flag("surface.debm_simple.albedo_input.periodic");
    m_input_albedo  = std::make_shared<array::Forcing>(m_grid,
                                                       file,
                                                       "surface_albedo",
                                                       "surface_albedo",
                                                       buffer_size,
                                                       periodic,
                                                       LINEAR);
  } else {
    m_input_albedo  = nullptr;
  }

  // initialize the spatially-variable air temperature standard deviation

  auto air_temp_sd = m_config->get_string("surface.debm_simple.std_dev.file");
  if (not air_temp_sd.empty()) {
    m_log->message(2, "  Reading standard deviation of near-surface air temperature from '%s'...\n",
                   air_temp_sd.c_str());

    int buffer_size = static_cast<int>(m_config->get_number("input.forcing.buffer_size"));

    File file(m_grid->com, air_temp_sd, io::PISM_GUESS, io::PISM_READONLY);

    bool periodic = m_config->get_flag("surface.debm_simple.std_dev.periodic");
    m_air_temp_sd = std::make_shared<array::Forcing>(m_grid, file, "air_temp_sd",
                                                     "", // no standard name
                                                     buffer_size, periodic, LINEAR);
  } else {
    double temp_std_dev = m_config->get_number("surface.debm_simple.std_dev");

    m_air_temp_sd = array::Forcing::Constant(m_grid, "air_temp_sd", temp_std_dev);
    m_log->message(2, "  Using constant standard deviation of near-surface air temperature.\n");
  }

  m_air_temp_sd->metadata(0)
      .long_name("standard deviation of near-surface air temperature")
      .units("Kelvin");

  m_mass_flux.metadata(0)
      .long_name("instantaneous surface mass balance (accumulation/ablation) rate")
      .units("kg m-2 s-1")
      .standard_name("land_ice_surface_specific_mass_balance_flux");
  m_mass_flux.metadata().set_string("comment", "positive values correspond to ice gain");

  // diagnostic fields:

  {
    m_accumulation = allocate_accumulation(g);
    m_melt         = allocate_melt(g);
    m_runoff       = allocate_runoff(g);
    m_temperature  = allocate_temperature(g);

    m_temperature_driven_melt.metadata(0)
        .long_name("temperature-driven melt in dEBM-simple")
        .units("kg m-2");
    m_temperature_driven_melt.set(0.0);

    m_insolation_driven_melt.metadata(0)
        .long_name("insolation-driven melt in dEBM-simple")
        .units("kg m-2");
    m_insolation_driven_melt.set(0.0);

    m_offset_melt.metadata(0)
        .long_name("offset melt in dEBM-simple")
        .units("kg m-2");
    m_offset_melt.set(0.0);
  }

  m_snow_depth.metadata(0)
      .long_name("snow cover depth (set to zero once a year)")
      .units("m");
  m_snow_depth.set(0.0);

  m_surface_albedo.metadata(0)
      .long_name("surface_albedo")
      .units("1")
      .standard_name("surface_albedo");
  m_surface_albedo.set(0.0);

  m_transmissivity.metadata(0)
      .long_name("atmosphere_transmissivity")
      .units("1");

  m_transmissivity.set(0.0);
}

void DEBMSimple::init_impl(const Geometry &geometry) {

  // call the default implementation (not the interface method init())
  SurfaceModel::init_impl(geometry);

  {
    m_log->message(2,
                   "* Initializing dEBM-simple, the diurnal Energy Balance Model (simple version).\n"
                   "  Inputs:  precipitation and 2m air temperature from an atmosphere model.\n"
                   "  Outputs: SMB and ice upper surface temperature.\n");
  }

  // initializing the model state
  InputOptions input = process_input_options(m_grid->com, m_config);

  auto default_albedo =  m_config->get_number("surface.debm_simple.albedo_max");
  if (input.type == INIT_RESTART) {
    m_snow_depth.read(input.filename, input.record);
    m_surface_albedo.read(input.filename, input.record);
  } else if (input.type == INIT_BOOTSTRAP) {
    m_snow_depth.regrid(input.filename, io::Default(0.0));
    m_surface_albedo.regrid(input.filename, io::Default(default_albedo));
  } else {
    m_snow_depth.set(0.0);
    m_surface_albedo.set(default_albedo);
  }

  {
    regrid("dEBM-simple surface model", m_snow_depth);
    regrid("dEBM-simple surface model", m_surface_albedo);
  }

  if ((bool)m_input_albedo) {
    auto filename = m_config->get_string("surface.debm_simple.albedo_input.file");
    bool periodic = m_config->get_flag("surface.debm_simple.albedo_input.periodic");
    m_input_albedo->init(filename, periodic);
  }

  // finish up
  {
    m_next_balance_year_start = compute_next_balance_year_start(time().current());

    m_accumulation->set(0.0);
    m_melt->set(0.0);
    m_runoff->set(0.0);
  }
}

MaxTimestep DEBMSimple::max_timestep_impl(double my_t) const {
  return m_atmosphere->max_timestep(my_t);
}

double DEBMSimple::compute_next_balance_year_start(double time) {
  // compute the time corresponding to the beginning of the next balance year
  double balance_year_start_day = m_config->get_number("surface.mass_balance_year_start_day"),
         one_day                = units::convert(m_sys, 1.0, "days", "seconds"),
         year_start             = this->time().calendar_year_start(time),
         balance_year_start     = year_start + (balance_year_start_day - 1.0) * one_day;

  if (balance_year_start > time) {
    return balance_year_start;
  }
  return this->time().increment_date(balance_year_start, 1);
}

/** @brief Extracts snow accumulation from mixed (snow and rain) precipitation using a
  *  temperature threshold with a linear transition.
  *
  * Rain is removed entirely from the surface mass balance, and will not be included in
  * the computed runoff, which is meltwater runoff.
  *
  * There is an linear transition for Tmin below which all precipitation is interpreted as
  * snow, and Tmax above which all precipitation is rain (see, e.g. [\ref Hock2005b]).
  *
  * Returns the *solid* (snow) accumulation *rate*.
  *
  * @param[in] T air temperature
  * @param[in] P precipitation rate
  */
double DEBMSimple::snow_accumulation(double T, double P) const {

  // do not allow negative precipitation
  if (P < 0.0) {
    return 0.0;
  }

  if (m_precip_as_snow or T <= m_Tmin) {
    // T <= Tmin, all precip is snow
    return P;
  }

  if (T < m_Tmax) { // linear transition from Tmin to Tmax
    return P * (m_Tmax - T) / (m_Tmax - m_Tmin);
  }

  // T >= Tmax, all precip is rain -- ignore it
  return 0.0;
}


void DEBMSimple::update_impl(const Geometry &geometry, double t, double dt) {

  const double melting_point = 273.15;

  // update to ensure that temperature and precipitation time series are correct:
  m_atmosphere->update(geometry, t, dt);

  // Use near-surface air temperature as the top-of-the-ice temperature:
  m_temperature->copy_from(m_atmosphere->air_temperature());

  // Set up air temperature and precipitation time series
  int N = static_cast<int>(timeseries_length(dt));

  const double dtseries = dt / N;
  std::vector<double> ts(N), T(N), S(N), P(N), Alb(N);
  std::vector<DEBMSimplePointwise::OrbitalParameters> orbital(N);

  for (int k = 0; k < N; ++k) {
    ts[k] = t + k * dtseries;

    // pre-compute orbital parameters which depend on time and *not* on the map-plane
    // location
    orbital[k] = m_model.orbital_parameters(ts[k]);
  }

  // update standard deviation time series
  m_air_temp_sd->update(t, dt);
  m_air_temp_sd->init_interpolation(ts);

  const auto &mask             = geometry.cell_type;
  const auto &H                = geometry.ice_thickness;
  const auto &surface_altitude = geometry.ice_surface_elevation;

  array::AccessScope list
    { &mask,
      &H,
      m_air_temp_sd.get(),
      &m_mass_flux,
      &m_snow_depth,
      m_accumulation.get(),
      m_melt.get(),
      m_runoff.get(),
      &m_temperature_driven_melt,
      &m_insolation_driven_melt,
      &m_offset_melt,
      &m_surface_albedo,
      &m_transmissivity,
      &geometry.latitude,
      &surface_altitude
    };

  if ((bool)m_input_albedo) {
    m_input_albedo->update(t, dt);
    m_input_albedo->init_interpolation(ts);
    list.add(*m_input_albedo);
  }

  double
    ice_density    = m_config->get_number("constants.ice.density"),
    sigmalapserate = m_config->get_number("surface.pdd.std_dev.lapse_lat_rate"),
    sigmabaselat   = m_config->get_number("surface.pdd.std_dev.lapse_lat_base");

  m_atmosphere->init_timeseries(ts);
  m_atmosphere->begin_pointwise_access();

  ParallelSection loop(m_grid->com);
  try {
    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      double latitude = geometry.latitude(i, j);

      // Get the temperature time series from an atmosphere model and its modifiers
      m_atmosphere->temp_time_series(i, j, T);

      if (mask.ice_free_water(i, j)) {
        // ignore precipitation over ice-free ocean
        for (int k = 0; k < N; ++k) {
          P[k] = 0.0;
        }
      } else {
        // elsewhere, get precipitation from the atmosphere model
        m_atmosphere->precip_time_series(i, j, P);

        // Use temperature time series to remove rainfall from precipitation and convert to
        // m/s ice equivalent.
        for (int k = 0; k < N; ++k) {
          P[k] = snow_accumulation(T[k],  // air temperature (input)
                                   P[k] / ice_density); // precipitation rate (input, gets overwritten)
        }
      }

      if ((bool)m_input_albedo) {
        m_input_albedo->interp(i, j, Alb);
      }

      // standard deviation of daily variability of air temperature
      {
        // interpolate temperature standard deviation time series
        //
        // Note: this works when m_air_temp_sd is constant in time.
        m_air_temp_sd->interp(i, j, S);

        if (sigmalapserate != 0.0) {
          // apply standard deviation lapse rate on top of prescribed values
          for (int k = 0; k < N; ++k) {
            S[k] += sigmalapserate * (latitude - sigmabaselat);
          }
          (*m_air_temp_sd)(i, j) = S[0]; // ensure correct SD reporting
        } else if (m_sd_use_param and mask.icy(i, j)) {
          // apply standard deviation parameterization over ice if in use
          for (int k = 0; k < N; ++k) {
            S[k] = std::max(m_sd_param_a * (T[k] - melting_point) + m_sd_param_b, 0.0);
          }
          (*m_air_temp_sd)(i, j) = S[0]; // ensure correct SD reporting
        }
      }

      {
        double next_snow_depth_reset = m_next_balance_year_start;

        // make copies of firn and snow depth values at this point to avoid accessing 2D
        // fields in the inner loop
        double
          ice_thickness = H(i, j),
          snow          = m_snow_depth(i, j),
          surfelev      = surface_altitude(i, j),
          albedo        = m_surface_albedo(i, j);

        auto cell_type = static_cast<cell_type::Value>(mask.as_int(i, j));

        double
          A   = 0.0,            // accumulation
          M   = 0.0,            // melt
          R   = 0.0,            // runoff
          SMB = 0.0,            // resulting mass balance
          Mi  = 0.0,            // insolation melt contribution
          Mt  = 0.0,            // temperature melt contribution
          Mc  = 0.0,            // offset melt contribution
          Al  = 0.0;            // albedo

        // beginning of the loop over small time steps:
        for (int k = 0; k < N; ++k) {

          if (ts[k] >= next_snow_depth_reset) {
            snow = 0.0;
            while (next_snow_depth_reset <= ts[k]) {
              next_snow_depth_reset = time().increment_date(next_snow_depth_reset, 1);
            }
          }

          auto accumulation = P[k] * dtseries;

          DEBMSimpleMelt melt_info{};
          if (not cell_type::ice_free_water(cell_type)) {

            melt_info = m_model.melt(orbital[k].declination,
                                     orbital[k].distance_factor,
                                     dtseries,
                                     S[k],
                                     T[k],
                                     surfelev,
                                     latitude,
                                     (bool)m_input_albedo ? Alb[k] : albedo);
          }

          auto changes = m_model.step(ice_thickness,
                                      melt_info.total_melt,
                                      snow,
                                      accumulation);

          if ((bool) m_input_albedo) {
            albedo = Alb[k];
          } else {
            albedo = m_model.albedo(changes.melt / dtseries, cell_type);
          }

          // update ice thickness
          ice_thickness += changes.smb;
          assert(ice_thickness >= 0);
          // update snow depth
          snow += changes.snow_depth;
          assert(snow >= 0);
          // update total accumulation, melt, and runoff
          {
            A   += accumulation;
            M   += changes.melt;
            Mt  += melt_info.temperature_melt;
            Mi  += melt_info.insolation_melt;
            Mc  += melt_info.offset_melt;
            R   += changes.runoff;
            SMB += changes.smb;
            Al  += albedo;
          }
        } // end of the time-stepping loop

        // set firn and snow depths
        m_snow_depth(i, j)     = snow;
        m_surface_albedo(i, j) = Al / N;
        m_transmissivity(i, j) = m_model.atmosphere_transmissivity(surfelev);

        // set melt terms at this point, converting
        // from "meters, ice equivalent" to "kg / m^2"
        m_temperature_driven_melt(i, j) = Mt * ice_density;
        m_insolation_driven_melt(i, j)  = Mi * ice_density;
        m_offset_melt(i, j)         = Mc * ice_density;

        // set total accumulation, melt, and runoff, and SMB at this point, converting
        // from "meters, ice equivalent" to "kg / m^2"
        {
          (*m_accumulation)(i, j) = A * ice_density;
          (*m_melt)(i, j)         = M * ice_density;
          (*m_runoff)(i, j)       = R * ice_density;
          // m_mass_flux (unlike m_accumulation, m_melt, and m_runoff), is a
          // rate. m * (kg / m^3) / second = kg / m^2 / second
          m_mass_flux(i, j) = SMB * ice_density / dt;
        }
      }

      if (mask.ice_free_water(i, j)) {
        m_snow_depth(i, j) = 0.0; // snow over the ocean does not stick
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_atmosphere->end_pointwise_access();

  m_next_balance_year_start =
    compute_next_balance_year_start(time().current());
}


const array::Scalar &DEBMSimple::mass_flux_impl() const {
  return m_mass_flux;
}

const array::Scalar &DEBMSimple::temperature_impl() const {
  return *m_temperature;
}

const array::Scalar &DEBMSimple::accumulation_impl() const {
  return *m_accumulation;
}

const array::Scalar &DEBMSimple::melt_impl() const {
  return *m_melt;
}

const array::Scalar &DEBMSimple::runoff_impl() const {
  return *m_runoff;
}

const array::Scalar &DEBMSimple::snow_depth() const {
  return m_snow_depth;
}

const array::Scalar &DEBMSimple::air_temp_sd() const {
  return *m_air_temp_sd;
}

const array::Scalar &DEBMSimple::insolation_driven_melt() const {
  return m_insolation_driven_melt;
}

const array::Scalar &DEBMSimple::temperature_driven_melt() const {
  return m_temperature_driven_melt;
}

const array::Scalar &DEBMSimple::offset_melt() const {
  return m_offset_melt;
}

const array::Scalar &DEBMSimple::surface_albedo() const {
  return m_surface_albedo;
}

const array::Scalar &DEBMSimple::atmosphere_transmissivity() const {
  return m_transmissivity;
}

void DEBMSimple::define_model_state_impl(const File &output) const {
  SurfaceModel::define_model_state_impl(output);
  m_snow_depth.define(output, io::PISM_DOUBLE);
  m_surface_albedo.define(output, io::PISM_DOUBLE);
}

void DEBMSimple::write_model_state_impl(const File &output) const {
  SurfaceModel::write_model_state_impl(output);
  m_snow_depth.write(output);
  m_surface_albedo.write(output);
}

const DEBMSimplePointwise& DEBMSimple::pointwise_model() const {
  return m_model;
}

namespace diagnostics {

/*! @brief Report mean top of atmosphere insolation */
class DEBMSInsolation : public Diag<DEBMSimple>
{
public:
  DEBMSInsolation(const DEBMSimple *m) : Diag<DEBMSimple>(m) {
    m_vars = { { m_sys, "insolation" } };
    m_vars[0]
        .long_name(
            "mean top of atmosphere insolation during the period when the sun is above the critical angle Phi")
        .units("W m-2");
  }

protected:
  std::shared_ptr<array::Array> compute_impl() const {

    auto result = allocate<array::Scalar>("insolation");

    const auto *latitude = m_grid->variables().get_2d_scalar("latitude");
    auto ctx = m_grid->ctx();

    {
      const auto& M = model->pointwise_model();

      auto orbital = M.orbital_parameters(ctx->time()->current());

      array::AccessScope list{latitude, result.get()};

      for (auto p = m_grid->points(); p; p.next()) {
        const int i = p.i(), j = p.j();

        (*result)(i, j) = M.insolation(orbital.declination,
                                       orbital.distance_factor,
                                       (*latitude)(i, j));
      }
    }

    return result;
  }
};

enum AmountKind { AMOUNT, MASS };

/*! @brief Report surface insolation melt, averaged over the reporting interval */
class DEBMSInsolationMelt : public DiagAverageRate<DEBMSimple> {
public:
  DEBMSInsolationMelt(const DEBMSimple *m, AmountKind kind)
    : DiagAverageRate<DEBMSimple>(m,
                                  kind == AMOUNT
                                  ? "debm_insolation_driven_melt_flux"
                                  : "debm_insolation_driven_melt_rate",
                                  TOTAL_CHANGE),
        m_kind(kind),
        m_melt_mass(m_grid, "debm_insolation_driven_melt_mass") {

    std::string
      name          = "debm_insolation_driven_melt_flux",
      long_name     = "surface insolation melt, averaged over the reporting interval",
      accumulator_units = "kg m-2",
      internal_units = "kg m-2 second-1",
      external_units = "kg m-2 year-1";
    if (kind == MASS) {
      name          = "debm_insolation_driven_melt_rate";
      accumulator_units = "kg";
      internal_units = "kg second-1";
      external_units = "Gt year-1";
    }

    m_accumulator.metadata().units(accumulator_units);

    m_vars = { { m_sys, name } };
    m_vars[0]
        .long_name(long_name)
        .units(internal_units)
        .output_units(external_units);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value, external_units, internal_units);
    m_vars[0].set_number("_FillValue", fill_value);
  }

protected:
  const array::Scalar &model_input() {
    const auto &melt_amount = model->insolation_driven_melt();

    if (m_kind == MASS) {
      m_melt_mass.copy_from(melt_amount);
      m_melt_mass.scale(m_grid->cell_area());
      return m_melt_mass;
    }

    return melt_amount;
  }

private:
  AmountKind m_kind;
  array::Scalar m_melt_mass;
};

/*! @brief Report surface temperature melt, averaged over the reporting interval */
class DEBMSTemperatureMelt : public DiagAverageRate<DEBMSimple> {
public:
  DEBMSTemperatureMelt(const DEBMSimple *m, AmountKind kind)
      : DiagAverageRate<DEBMSimple>(m,
                                    kind == AMOUNT
                                    ? "debm_temperature_driven_melt_flux"
                                    : "debm_temperature_driven_melt_rate",
                                    TOTAL_CHANGE),
        m_kind(kind),
        m_melt_mass(m_grid, "temperature_melt_mass") {

    std::string
      name          = "debm_temperature_driven_melt_flux",
      long_name     = "temperature-driven melt, averaged over the reporting interval",
      accumulator_units = "kg m-2",
      internal_units = "kg m-2 second-1",
      external_units = "kg m-2 year-1";
    if (kind == MASS) {
      name          = "debm_temperature_driven_melt_rate";
      accumulator_units = "kg";
      internal_units = "kg second-1";
      external_units = "Gt year-1";
    }

    m_accumulator.metadata().units(accumulator_units);

    m_vars = { { m_sys, name } };
    m_vars[0]
        .long_name(long_name)
        .units(internal_units)
        .output_units(external_units);
    m_vars[0].set_string("cell_methods", "time: mean");
    m_vars[0].set_number("_FillValue", to_internal(m_fill_value));
  }

protected:
  const array::Scalar &model_input() {
    const auto &melt_amount = model->temperature_driven_melt();

    if (m_kind == MASS) {
      m_melt_mass.copy_from(melt_amount);
      m_melt_mass.scale(m_grid->cell_area());
      return m_melt_mass;
    }

    return melt_amount;
  }

private:
  AmountKind m_kind;
  array::Scalar m_melt_mass;
};

/*! @brief Report surface backround melt, averaged over the reporting interval */
class DEBMSBackroundMelt : public DiagAverageRate<DEBMSimple> {
public:
  DEBMSBackroundMelt(const DEBMSimple *m, AmountKind kind)
      : DiagAverageRate<DEBMSimple>(m,
                                    kind == AMOUNT
                                    ? "debm_offset_melt_flux"
                                    : "debm_offset_melt_rate",
                                    TOTAL_CHANGE),
        m_kind(kind),
        m_melt_mass(m_grid, "backround_melt_mass") {

    std::string name              = "debm_offset_melt_flux",
                long_name         = "offset melt, averaged over the reporting interval",
                accumulator_units = "kg m-2",
                internal_units    = "kg m-2 second-1",
                external_units    = "kg m-2 year-1";

    if (kind == MASS) {
      name              = "debm_offset_melt_rate";
      accumulator_units = "kg";
      internal_units    = "kg second-1";
      external_units    = "Gt year-1";
    }
    m_accumulator.metadata().units(accumulator_units);

    m_vars = { { m_sys, name } };
    m_vars[0]
        .long_name(long_name)
        .units(internal_units)
        .output_units(external_units);
    m_vars[0].set_string("cell_methods", "time: mean");
    m_vars[0].set_number("_FillValue", to_internal(m_fill_value));
  }

protected:
  const array::Scalar &model_input() {
    const auto &melt_amount = model->offset_melt();

    if (m_kind == MASS) {
      m_melt_mass.copy_from(melt_amount);
      m_melt_mass.scale(m_grid->cell_area());
      return m_melt_mass;
    }

    return melt_amount;
  }

private:
  AmountKind m_kind;
  array::Scalar m_melt_mass;
};

} // end of namespace diagnostics

/*! @brief The number of points for temperature and precipitation time-series.
 */
unsigned int DEBMSimple::timeseries_length(double dt) const {
  double dt_years = dt / m_year_length;

  return std::max(1U, static_cast<unsigned int>(ceil(m_n_per_year * dt_years)));
}

DiagnosticList DEBMSimple::diagnostics_impl() const {
  using namespace diagnostics;

  DiagnosticList result = {
    { "debm_insolation_driven_melt_flux", Diagnostic::Ptr(new DEBMSInsolationMelt(this, AMOUNT)) },
    { "debm_insolation_driven_melt_rate", Diagnostic::Ptr(new DEBMSInsolationMelt(this, MASS)) },
    { "debm_temperature_driven_melt_flux", Diagnostic::Ptr(new DEBMSTemperatureMelt(this, AMOUNT)) },
    { "debm_temperature_driven_melt_rate", Diagnostic::Ptr(new DEBMSTemperatureMelt(this, MASS)) },
    { "debm_offset_melt_flux", Diagnostic::Ptr(new DEBMSBackroundMelt(this, AMOUNT)) },
    { "debm_offset_melt_rate", Diagnostic::Ptr(new DEBMSBackroundMelt(this, MASS)) },
    { "air_temp_sd", Diagnostic::wrap(*m_air_temp_sd) },
    { "snow_depth", Diagnostic::wrap(m_snow_depth) },
    { "surface_albedo", Diagnostic::wrap(m_surface_albedo) },
    { "atmosphere_transmissivity", Diagnostic::wrap(m_transmissivity) },
    { "insolation", Diagnostic::Ptr(new DEBMSInsolation(this)) }
  };

  result = pism::combine(result, SurfaceModel::diagnostics_impl());

  return result;
}

} // end of namespace surface
} // end of namespace pism
