// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2022 PISM Authors
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

#include "DEBMSimple.hh"

#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Time.hh"
#include "pism/util/io/File.hh"

#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/iceModelVec2T.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace surface {

///// PISM surface model implementing a dEBM-Simple scheme.

DEBMSimple::DEBMSimple(IceGrid::ConstPtr g, std::shared_ptr<atmosphere::AtmosphereModel> input)
  : SurfaceModel(g, std::move(input)),
      m_model(*g->ctx()),
      m_mass_flux(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS),
      m_firn_depth(m_grid, "firn_depth", WITHOUT_GHOSTS),
      m_snow_depth(m_grid, "snow_depth", WITHOUT_GHOSTS),
      m_temperature_driven_melt(m_grid, "debms_temperature_driven_melt_flux", WITHOUT_GHOSTS),
      m_insolation_driven_melt(m_grid, "debms_insolation_driven_melt_flux", WITHOUT_GHOSTS),
      m_background_melt(m_grid, "debms_background_melt_flux", WITHOUT_GHOSTS),
      m_surface_albedo(m_grid, "surface_albedo", WITHOUT_GHOSTS),
      m_transmissivity(m_grid, "atmosphere_transmissivity", WITHOUT_GHOSTS),
      m_insolation(m_grid, "insolation", WITHOUT_GHOSTS) {

  m_sd_use_param = m_config->get_flag("surface.debm_simple.std_dev_param.enabled");
  m_sd_param_a   = m_config->get_number("surface.debm_simple.std_dev_param.a");
  m_sd_param_b   = m_config->get_number("surface.debm_simple.std_dev_param.b");

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

    File file(m_grid->com, albedo_input, PISM_GUESS, PISM_READONLY);

    int buffer_size = static_cast<int>(m_config->get_number("input.forcing.buffer_size"));
    bool periodic = m_config->get_flag("surface.debm_simple.albedo_input.periodic");
    m_input_albedo  = IceModelVec2T::ForcingField(m_grid,
                                                  file,
                                                  "albedo",
                                                  "surface_albedo",
                                                  buffer_size,
                                                  periodic,
                                                  LINEAR);
  } else {
    m_input_albedo  = nullptr;
  }

  // initialize the spatially-variable air temperature standard deviation

  ForcingOptions air_temp_sd(*m_grid->ctx(), "surface.debm_simple.std_dev");
  if (not air_temp_sd.filename.empty()) {
    m_log->message(2, "  Reading standard deviation of near-surface air temperature from '%s'...\n",
                   air_temp_sd.filename.c_str());

    int buffer_size = static_cast<int>(m_config->get_number("input.forcing.buffer_size"));

    File file(m_grid->com, air_temp_sd.filename, PISM_GUESS, PISM_READONLY);

    m_air_temp_sd          = IceModelVec2T::ForcingField(m_grid, file, "air_temp_sd",
                                                         "", // no standard name
                                                         buffer_size, air_temp_sd.periodic, LINEAR);
    m_use_air_temp_sd_file = true;
  } else {
    double temp_std_dev = m_config->get_number("surface.debm_simple.std_dev");

    m_air_temp_sd = IceModelVec2T::Constant(m_grid, "air_temp_sd", temp_std_dev);
    m_log->message(2, "  Using constant standard deviation of near-surface air temperature.\n");
    m_use_air_temp_sd_file = false;
  }

  m_air_temp_sd->set_attrs("climate_forcing", "standard deviation of near-surface air temperature", "Kelvin", "Kelvin",
                           "", 0);

  m_mass_flux.set_attrs("diagnostic", "instantaneous surface mass balance (accumulation/ablation) rate", "kg m-2 s-1",
                        "kg m-2 s-1", "land_ice_surface_specific_mass_balance_flux", 0);
  m_mass_flux.metadata().set_string("comment", "positive values correspond to ice gain");

  // diagnostic fields:

  {
    m_accumulation = allocate_accumulation(g);
    m_melt         = allocate_melt(g);
    m_runoff       = allocate_runoff(g);
    m_temperature  = allocate_temperature(g);

    m_temperature_driven_melt.set_attrs("diagnostic",
                                        "temperature-driven melt in dEBM-simple",
                                        "kg m-2", "kg m-2", "", 0);
    m_temperature_driven_melt.set(0.0);

    m_insolation_driven_melt.set_attrs("diagnostic",
                                       "insolation-driven melt in dEBM-simple",
                                       "kg m-2", "kg m-2", "", 0);
    m_insolation_driven_melt.set(0.0);

    m_background_melt.set_attrs("diagnostic",
                                "background melt in dEBM-simple",
                                "kg m-2", "kg m-2", "", 0);
    m_background_melt.set(0.0);
  }

  m_snow_depth.set_attrs("diagnostic", "snow cover depth (set to zero once a year)", "m", "m", "", 0);
  m_snow_depth.set(0.0);

  m_firn_depth.set_attrs("diagnostic", "firn cover depth", "m", "m", "", 0);
  m_firn_depth.metadata().set_number("valid_min", 0.0);
  m_firn_depth.set(0.0);

  m_surface_albedo.set_attrs("diagnostic", "surface_albedo", "1", "1", "surface_albedo", 0);
  m_surface_albedo.set(0.0);

  m_transmissivity.set_attrs("diagnostic", "atmosphere_transmissivity", "", "", "", 0);
  m_transmissivity.set(0.0);

  m_insolation.set_attrs("diagnostic",
                         "average topf of the atmosphere insolation "
                         "during the period when the sun is above the critical angle Phi",
                         "W m-2",
                         "W m-2",
                         "",
                         0);
  m_insolation.set(0.0);
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

  std::string firn_file = m_config->get_string("surface.debm_simple.firn_depth_file");

  if (input.type == INIT_RESTART) {
    if (not firn_file.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "surface.debm_simple.firn_depth_file is not allowed when"
                                    " re-starting from a PISM output file.");
    }

    m_firn_depth.read(input.filename, input.record);
    m_snow_depth.read(input.filename, input.record);
    m_surface_albedo.read(input.filename, input.record);
  } else if (input.type == INIT_BOOTSTRAP) {

    m_snow_depth.regrid(input.filename, OPTIONAL, 0.0);
    m_surface_albedo.regrid(input.filename, OPTIONAL, m_config->get_number("surface.debm_simple.albedo_snow"));

    if (firn_file.empty()) {
      m_firn_depth.regrid(input.filename, OPTIONAL, 0.0);
    } else {
      m_firn_depth.regrid(firn_file, CRITICAL);
    }
  } else {

    m_snow_depth.set(0.0);
    m_surface_albedo.set(m_config->get_number("surface.debm_simple.albedo_snow"));

    if (firn_file.empty()) {
      m_firn_depth.set(0.0);
    } else {
      m_firn_depth.regrid(firn_file, CRITICAL);
    }
  }

  {
    regrid("dEBM-simple surface model", m_snow_depth);
    regrid("dEBM-simple surface model", m_firn_depth);
    regrid("dEBM-simple surface model", m_surface_albedo);
  }

  // finish up
  {
    m_next_balance_year_start = compute_next_balance_year_start(m_grid->ctx()->time()->current());

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
  double balance_year_start_day = m_config->get_number("surface.pdd.balance_year_start_day"),
         one_day                = units::convert(m_sys, 1.0, "days", "seconds"),
         year_start             = m_grid->ctx()->time()->calendar_year_start(time),
         balance_year_start     = year_start + (balance_year_start_day - 1.0) * one_day;

  if (balance_year_start > time) {
    return balance_year_start;
  }
  return m_grid->ctx()->time()->increment_date(balance_year_start, 1);
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

  for (int k = 0; k < N; ++k) {
    ts[k] = t + k * dtseries;
  }

  // update standard deviation time series
  if (m_use_air_temp_sd_file) {
    m_air_temp_sd->update(t, dt);
    m_air_temp_sd->init_interpolation(ts);
  }

  const auto &mask             = geometry.cell_type;
  const auto &H                = geometry.ice_thickness;
  const auto &surface_altitude = geometry.ice_surface_elevation;

  IceModelVec::AccessList list
    { &mask,
      &H,
      m_air_temp_sd.get(),
      &m_mass_flux,
      &m_firn_depth,
      &m_snow_depth,
      m_accumulation.get(),
      m_melt.get(),
      m_runoff.get(),
      &m_temperature_driven_melt,
      &m_insolation_driven_melt,
      &m_background_melt,
      &m_surface_albedo,
      &m_transmissivity,
      &m_insolation,
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
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double latitude = geometry.latitude(i, j);

      // Get the temperature time series from an atmosphere model and its modifiers
      m_atmosphere->temp_time_series(i, j, T);

      if (mask.ice_free_ocean(i, j)) {
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
          P[i] = snow_accumulation(T[i],  // air temperature (input)
                                   P[i] / ice_density); // precipitation rate (input, gets overwritten)
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
            S[k] = m_sd_param_a * (T[k] - melting_point) + m_sd_param_b;
            if (S[k] < 0.0) {
              S[k] = 0.0;
            }
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
          firn          = m_firn_depth(i, j),
          snow          = m_snow_depth(i, j),
          surfelev      = surface_altitude(i, j),
          albedo        = m_surface_albedo(i, j);

        auto cell_type = static_cast<MaskValue>(mask.as_int(i, j));

        double
          A   = 0.0,            // accumulation
          M   = 0.0,            // melt
          R   = 0.0,            // runoff
          SMB = 0.0,            // resulting mass balance
          Mi  = 0.0,            // insolation melt contribution
          Mt  = 0.0,            // temperature melt contribution
          Mc  = 0.0,            // background melt contribution
          Qi  = 0.0,            // insolation averaged over \Delta t_Phi
          Al  = 0.0;            // albedo

        // beginning of the loop over small time steps:
        for (int k = 0; k < N; ++k) {

          if (ts[k] >= next_snow_depth_reset) {
            snow = 0.0;
            while (next_snow_depth_reset <= ts[k]) {
              next_snow_depth_reset = m_grid->ctx()->time()->increment_date(next_snow_depth_reset, 1);
            }
          }

          auto accumulation = P[k] * dtseries;

          DEBMSimpleMelt melt_info{};
          if (not mask::ice_free_ocean(cell_type)) {
            melt_info = m_model.melt(ts[k],
                                     dtseries,
                                     S[k],
                                     T[k],
                                     surfelev,
                                     latitude,
                                     (bool)m_input_albedo ? Alb[k] : albedo);
          }

          auto changes = m_model.step(ice_thickness,
                                      melt_info.total_melt,
                                      firn,
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
            Mt  += melt_info.temperature_melt;
            Mi  += melt_info.insolation_melt;
            Mc  += melt_info.background_melt;
            R   += changes.runoff;
            SMB += changes.smb;
            Qi  += melt_info.insolation;
            Al  += albedo;
          }
        } // end of the time-stepping loop

        // set firn and snow depths
        m_firn_depth(i, j)     = firn;
        m_snow_depth(i, j)     = snow;
        m_surface_albedo(i, j) = Al / N;
        m_transmissivity(i, j) = m_model.atmosphere_transmissivity(surfelev);
        m_insolation(i, j)     = Qi / N;

        // set melt terms at this point, converting
        // from "meters, ice equivalent" to "kg / m^2"
        m_temperature_driven_melt(i, j)  = Mt * ice_density;
        m_insolation_driven_melt(i, j) = Mi * ice_density;
        m_background_melt(i, j)     = Mc * ice_density;

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

      if (mask.ice_free_ocean(i, j)) {
        m_firn_depth(i, j) = 0.0; // no firn in the ocean
        m_snow_depth(i, j) = 0.0; // snow over the ocean does not stick
      }
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  m_atmosphere->end_pointwise_access();

  m_next_balance_year_start =
    compute_next_balance_year_start(m_grid->ctx()->time()->current());
}


const IceModelVec2S &DEBMSimple::mass_flux_impl() const {
  return m_mass_flux;
}

const IceModelVec2S &DEBMSimple::temperature_impl() const {
  return *m_temperature;
}

const IceModelVec2S &DEBMSimple::accumulation_impl() const {
  return *m_accumulation;
}

const IceModelVec2S &DEBMSimple::melt_impl() const {
  return *m_melt;
}

const IceModelVec2S &DEBMSimple::runoff_impl() const {
  return *m_runoff;
}

const IceModelVec2S &DEBMSimple::firn_depth() const {
  return m_firn_depth;
}

const IceModelVec2S &DEBMSimple::snow_depth() const {
  return m_snow_depth;
}

const IceModelVec2S &DEBMSimple::air_temp_sd() const {
  return *m_air_temp_sd;
}

const IceModelVec2S &DEBMSimple::insolation_driven_melt() const {
  return m_insolation_driven_melt;
}

const IceModelVec2S &DEBMSimple::temperature_driven_melt() const {
  return m_temperature_driven_melt;
}

const IceModelVec2S &DEBMSimple::background_melt() const {
  return m_background_melt;
}

const IceModelVec2S &DEBMSimple::surface_albedo() const {
  return m_surface_albedo;
}

const IceModelVec2S &DEBMSimple::atmosphere_transmissivity() const {
  return m_transmissivity;
}

const IceModelVec2S &DEBMSimple::insolation() const {
  return m_insolation;
}

void DEBMSimple::define_model_state_impl(const File &output) const {
  SurfaceModel::define_model_state_impl(output);
  m_firn_depth.define(output, PISM_DOUBLE);
  m_snow_depth.define(output, PISM_DOUBLE);
  m_surface_albedo.define(output, PISM_DOUBLE);
}

void DEBMSimple::write_model_state_impl(const File &output) const {
  SurfaceModel::write_model_state_impl(output);
  m_firn_depth.write(output);
  m_snow_depth.write(output);
  m_surface_albedo.write(output);
}

namespace diagnostics {

enum AmountKind { AMOUNT, MASS };

/*! @brief Report surface insolation melt, averaged over the reporting interval */
class DEBMSInsolationMelt : public DiagAverageRate<DEBMSimple> {
public:
  DEBMSInsolationMelt(const DEBMSimple *m, AmountKind kind)
    : DiagAverageRate<DEBMSimple>(m,
                                  kind == AMOUNT
                                  ? "debms_insolation_driven_melt_flux"
                                  : "debms_insolation_driven_melt_flux",
                                  TOTAL_CHANGE),
        m_kind(kind),
        m_melt_mass(m_grid, "debm_insolation_driven_melt_mass", WITHOUT_GHOSTS) {

    std::string
      name          = "debms_insolation_driven_melt_flux",
      long_name     = "surface insolation melt, averaged over the reporting interval",
      standard_name = "",
      accumulator_units = "kg m-2",
      internal_units = "kg m-2 second-1",
      external_units = "kg m-2 year-1";
    if (kind == MASS) {
      name          = "debms_insolation_driven_melt_rate";
      standard_name = "";
      accumulator_units = "kg";
      internal_units = "kg second-1";
      external_units = "Gt year-1";
    }

    m_vars = { SpatialVariableMetadata(m_sys, name) };
    m_accumulator.metadata().set_string("units", accumulator_units);

    set_attrs(long_name, standard_name, internal_units, external_units, 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value, external_units, internal_units);
    m_vars[0].set_number("_FillValue", fill_value);
  }

protected:
  const IceModelVec2S &model_input() {
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
  IceModelVec2S m_melt_mass;
};

/*! @brief Report surface temperature melt, averaged over the reporting interval */
class DEBMSTemperatureMelt : public DiagAverageRate<DEBMSimple> {
public:
  DEBMSTemperatureMelt(const DEBMSimple *m, AmountKind kind)
      : DiagAverageRate<DEBMSimple>(m,
                                    kind == AMOUNT
                                    ? "debms_temperature_driven_melt_flux"
                                    : "debms_temperature_driven_melt_rate",
                                    TOTAL_CHANGE),
        m_kind(kind),
        m_melt_mass(m_grid, "temperature_melt_mass", WITHOUT_GHOSTS) {

    std::string
      name          = "debms_temperature_driven_melt_flux",
      long_name     = "temperature-driven melt, averaged over the reporting interval",
      standard_name = "",
      accumulator_units = "kg m-2",
      internal_units = "kg m-2 second-1",
      external_units = "kg m-2 year-1";
    if (kind == MASS) {
      name          = "debms_temperature_driven_melt_rate";
      standard_name = "";
      accumulator_units = "kg";
      internal_units = "kg second-1";
      external_units = "Gt year-1";
    }

    m_vars = { SpatialVariableMetadata(m_sys, name) };
    m_accumulator.metadata().set_string("units", accumulator_units);

    set_attrs(long_name, standard_name, internal_units, external_units, 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value, external_units, internal_units);
    m_vars[0].set_number("_FillValue", fill_value);
  }

protected:
  const IceModelVec2S &model_input() {
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
  IceModelVec2S m_melt_mass;
};

/*! @brief Report surface backround melt, averaged over the reporting interval */
class DEBMSBackroundMelt : public DiagAverageRate<DEBMSimple> {
public:
  DEBMSBackroundMelt(const DEBMSimple *m, AmountKind kind)
      : DiagAverageRate<DEBMSimple>(m,
                                    kind == AMOUNT
                                    ? "debms_background_melt_flux"
                                    : "debms_background_melt_rate",
                                    TOTAL_CHANGE),
        m_kind(kind),
        m_melt_mass(m_grid, "backround_melt_mass", WITHOUT_GHOSTS) {

    std::string name          = "debms_background_melt_flux",
                long_name     = "background melt, averaged over the reporting interval",
                standard_name = "debms_background_melt_flux", accumulator_units = "kg m-2",
                internal_units = "kg m-2 second-1", external_units = "kg m-2 year-1";
    if (kind == MASS) {
      name          = "debms_background_melt_rate";
      standard_name = "", accumulator_units = "kg", internal_units = "kg second-1";
      external_units = "Gt year-1";
    }

    m_vars = { SpatialVariableMetadata(m_sys, name) };
    m_accumulator.metadata().set_string("units", accumulator_units);

    set_attrs(long_name, standard_name, internal_units, external_units, 0);
    m_vars[0].set_string("cell_methods", "time: mean");

    double fill_value = units::convert(m_sys, m_fill_value, external_units, internal_units);
    m_vars[0].set_number("_FillValue", fill_value);
  }

protected:
  const IceModelVec2S &model_input() {
    const auto &melt_amount = model->background_melt();

    if (m_kind == MASS) {
      m_melt_mass.copy_from(melt_amount);
      m_melt_mass.scale(m_grid->cell_area());
      return m_melt_mass;
    }

    return melt_amount;
  }

private:
  AmountKind m_kind;
  IceModelVec2S m_melt_mass;
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
    { "debms_insolation_driven_melt_flux", Diagnostic::Ptr(new DEBMSInsolationMelt(this, AMOUNT)) },
    { "debms_insolation_driven_melt_rate", Diagnostic::Ptr(new DEBMSInsolationMelt(this, MASS)) },
    { "debms_temperature_driven_melt_flux", Diagnostic::Ptr(new DEBMSTemperatureMelt(this, AMOUNT)) },
    { "debms_temperature_driven_melt_rate", Diagnostic::Ptr(new DEBMSTemperatureMelt(this, MASS)) },
    { "debms_background_melt_flux", Diagnostic::Ptr(new DEBMSBackroundMelt(this, AMOUNT)) },
    { "debms_background_melt_rate", Diagnostic::Ptr(new DEBMSBackroundMelt(this, MASS)) },
    { "air_temp_sd", Diagnostic::wrap(*m_air_temp_sd) },
    { "snow_depth", Diagnostic::wrap(m_snow_depth) },
    { "firn_depth", Diagnostic::wrap(m_firn_depth) },
    { "surface_albedo", Diagnostic::wrap(m_surface_albedo) },
    { "atmosphere_transmissivity", Diagnostic::wrap(m_transmissivity) },
    { "insolation", Diagnostic::wrap(m_insolation) }
  };

  result = pism::combine(result, SurfaceModel::diagnostics_impl());

  return result;
}

} // end of namespace surface
} // end of namespace pism
