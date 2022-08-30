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
#include "localITM.hh"
#include "localMassBalance.hh"
#include "pism/coupler/AtmosphereModel.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Time.hh"
#include "pism/util/Vars.hh"
#include "pism/util/io/File.hh"
#include "pism/util/pism_options.hh"

#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/IceModelVec2CellType.hh"
#include "pism/util/ScalarForcing.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/iceModelVec2T.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {
namespace surface {

///// PISM surface model implementing a dEBM scheme.

DEBMSimple::DEBMSimple(IceGrid::ConstPtr g, std::shared_ptr<atmosphere::AtmosphereModel> input)
    : SurfaceModel(g, input),
      m_mbscheme(*m_config, m_sys),
      m_mass_flux(m_grid, "climatic_mass_balance", WITHOUT_GHOSTS),
      m_firn_depth(m_grid, "firn_depth", WITHOUT_GHOSTS),
      m_snow_depth(m_grid, "snow_depth", WITHOUT_GHOSTS),
      m_tempmelt(m_grid, "surface_temperature_melt_flux", WITHOUT_GHOSTS),
      m_insolmelt(m_grid, "surface_insolation_melt_flux", WITHOUT_GHOSTS),
      m_cmelt(m_grid, "surface_offset_melt_flux", WITHOUT_GHOSTS),
      m_albedo(m_grid, "albedo", WITHOUT_GHOSTS),
      m_transmissivity(m_grid, "transmissivity", WITHOUT_GHOSTS),
      m_TOAinsol(m_grid, "TOAinsol", WITHOUT_GHOSTS),
      m_qinsol(m_grid, "qinsol", WITHOUT_GHOSTS) {

  m_base_pddStdDev    = m_config->get_number("surface.itm.std_dev");
  m_sd_use_param      = m_config->get_flag("surface.itm.std_dev_use_param");
  m_sd_param_a        = m_config->get_number("surface.itm.std_dev_param_a");
  m_sd_param_b        = m_config->get_number("surface.itm.std_dev_param_b");
  m_refreeze_fraction = m_config->get_number("surface.itm.refreeze");

  ForcingOptions albedo_input(*m_grid->ctx(), "surface.itm.albedo_input");
  if (not albedo_input.filename.empty()) {
    m_log->message(2, " Albedo is read in from %s...", albedo_input.filename.c_str());

    File file(m_grid->com, albedo_input.filename, PISM_GUESS, PISM_READONLY);

    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");
    m_input_albedo           = IceModelVec2T::ForcingField(m_grid, file, "albedo",
                                                           "", // no standard name
                                                           buffer_size, albedo_input.periodic, LINEAR);
  } else {
    m_input_albedo = nullptr;
  }

  // initialize the spatially-variable air temperature standard deviation

  ForcingOptions air_temp_sd(*m_grid->ctx(), "surface.itm.std_dev");
  if (not air_temp_sd.filename.empty()) {
    m_log->message(2, "  Reading standard deviation of near-surface air temperature from '%s'...\n",
                   air_temp_sd.filename.c_str());

    int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, air_temp_sd.filename, PISM_GUESS, PISM_READONLY);

    m_air_temp_sd          = IceModelVec2T::ForcingField(m_grid, file, "air_temp_sd",
                                                         "", // no standard name
                                                         buffer_size, air_temp_sd.periodic, LINEAR);
    m_use_air_temp_sd_file = true;
  } else {
    m_air_temp_sd = IceModelVec2T::Constant(m_grid, "air_temp_sd", m_base_pddStdDev);
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

    m_tempmelt.set_attrs("diagnostic", "surface temp melt", "kg m-2", "kg m-2", "", 0);
    m_tempmelt.set(0.0);

    m_insolmelt.set_attrs("diagnostic", "surface insol melt", "kg m-2", "kg m-2", "", 0);
    m_insolmelt.set(0.0);

    m_cmelt.set_attrs("diagnostic", "surface c melt", "kg m-2", "kg m-2", "", 0);
    m_cmelt.set(0.0);
  }

  m_snow_depth.set_attrs("diagnostic", "snow cover depth (set to zero once a year)", "m", "m", "", 0);
  m_snow_depth.set(0.0);

  m_firn_depth.set_attrs("diagnostic", "firn cover depth", "m", "m", "", 0);
  m_firn_depth.metadata().set_number("valid_min", 0.0);
  m_firn_depth.set(0.0);

  m_albedo.set_attrs("diagnostic", "albedo", "", "", "", 0);
  m_albedo.set(0.0);

  m_transmissivity.set_attrs("diagnostic", "transmissivity", "", "", "", 0);
  m_transmissivity.set(0.0);

  m_TOAinsol.set_attrs("diagnostic", "insolation at the top of the atmosphere", "W m-2", "W m-2", "", 0);
  m_TOAinsol.set(0.0);

  m_qinsol.set_attrs("diagnostic",
                     "insolation at the top of the atmosphere, when the sun is above the elvation angle Phi", "W m-2",
                     "W m-2", "", 0);
  m_qinsol.set(0.0);

  std::string paleo_file = m_config->get_string("surface.itm.paleo.file");

  if (not paleo_file.empty()) {
    m_use_paleo_file = true;

    m_eccentricity.reset(
        new ScalarForcing(*g->ctx(), "surface.itm.paleo", "eccentricity", "", "", "eccentricity of the earth"));

    m_obliquity.reset(
        new ScalarForcing(*g->ctx(), "surface.itm.paleo", "obliquity", "degree", "degree", "obliquity of the earth"));

    m_perihelion_longitude.reset(
        new ScalarForcing(*g->ctx(), "surface.itm.paleo", "long_peri", "degree", "degree", "longitude of the perihelion"));
  } else {
    m_use_paleo_file = false;
  }
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

  std::string firn_file = m_config->get_string("surface.itm.firn_depth_file");

  if (input.type == INIT_RESTART) {
    if (not firn_file.empty()) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "surface.itm.firn_depth_file is not allowed when"
                                                         " re-starting from a PISM output file.");
    }

    m_firn_depth.read(input.filename, input.record);
    m_snow_depth.read(input.filename, input.record);
    m_albedo.read(input.filename, input.record);
  } else if (input.type == INIT_BOOTSTRAP) {

    m_snow_depth.regrid(input.filename, OPTIONAL, 0.0);
    m_albedo.regrid(input.filename, OPTIONAL, m_config->get_number("surface.itm.albedo_snow"));

    if (firn_file.empty()) {
      m_firn_depth.regrid(input.filename, OPTIONAL, 0.0);
    } else {
      m_firn_depth.regrid(firn_file, CRITICAL);
    }
  } else {

    m_snow_depth.set(0.0);
    m_albedo.set(m_config->get_number("surface.itm.albedo_snow"));

    if (firn_file.empty()) {
      m_firn_depth.set(0.0);
    } else {
      m_firn_depth.regrid(firn_file, CRITICAL);
    }
  }

  {
    regrid("ITM surface model", m_snow_depth);
    regrid("ITM surface model", m_firn_depth);
  }
  const bool force_albedo = m_config->get_flag("surface.itm.anomaly");
  if (force_albedo)
    m_log->message(2, " Albedo forcing sets summer albedo values to lower value\n");
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


bool DEBMSimple::albedo_anomaly_true(double time) {
  // This function is only here to perform darkening experiments, where the albedo over the whole ice sheet is reduced artificially during the summer months.

  // compute the time corresponding to the beginning of the darkening
  double anomaly_start_day = m_config->get_number("surface.itm.anomaly_start_day"),
         anomaly_end_day   = m_config->get_number("surface.itm.anomaly_end_day"),
         one_day           = units::convert(m_sys, 1.0, "days", "seconds"),
         one_year          = units::convert(m_sys, 1.0, "years", "seconds"),
         year_start        = m_grid->ctx()->time()->calendar_year_start(time),
         anomaly_start     = year_start + (anomaly_start_day - 1.0) * one_day,
         anomaly_end       = year_start + (anomaly_end_day - 1.0) * one_day;

  int frequency = m_config->get_number("surface.itm.anomaly_frequency");

  double period_seconds = frequency * one_year;
  double tmp            = time - floor(time / period_seconds) * period_seconds;
  bool frequency_true   = (tmp < one_year);

  if (time >= anomaly_start and time <= anomaly_end and frequency_true) {
    return true;
  }

  return false;
}


double DEBMSimple::get_distance2(double time) {
  // get the distance between earth and sun
  double a0 = 1.000110, a1 = 0.034221, a2 = 0.000719, b0 = 0., b1 = 0.001280, b2 = 0.000077, distance2 = 1.;

  double t  = 2. * M_PI * m_grid->ctx()->time()->year_fraction(time);
  distance2 = a0 + b0 + a1 * cos(t) + b1 * sin(t) + a2 * cos(2. * t) + b2 * sin(2. * t);
  // Equation 2.2.9 from Liou (2002)
  return distance2;
}


double DEBMSimple::get_delta(double time) {
  // get the earth declination delta
  double a0 = 0.006918, a1 = -0.399912, a2 = -0.006758, a3 = -0.002697, b0 = 0., b1 = 0.070257, b2 = 0.000907,
         b3 = 0.000148, delta = 1.;

  double t = 2. * M_PI * m_grid->ctx()->time()->year_fraction(time);
  delta =
      a0 + b0 + a1 * cos(t) + b1 * sin(t) + a2 * cos(2. * t) + b2 * sin(2. * t) + a3 * cos(3. * t) + b3 * sin(3. * t);
  // Equation 2.2.10 from Liou (2002)
  return delta;
}


double DEBMSimple::get_distance2_paleo(double time) {
  // for now the orbital parameters are as config parameters, but it would be best, if I could read in a time series
  double lambda   = get_lambda_paleo(time);
  double ecc      = 0;
  double peri_deg = 0;
  if (m_use_paleo_file) {
    ecc      = m_eccentricity->value(time);
    peri_deg = m_perihelion_longitude->value(time);
  } else {
    ecc      = m_config->get_number("surface.itm.paleo.eccentricity");
    peri_deg = m_config->get_number("surface.itm.paleo.long_peri");
  }
  double distance2 = pow((1. - ecc * cos(lambda - peri_deg * M_PI / 180.)), 2) / pow((1. - ecc * ecc), 2);
  // From Equation 2.2.5 Liou (2002)
  // (a/r)^2
  return distance2;
}


double DEBMSimple::get_delta_paleo(double time) {
  // for now the orbital parameters are as config parameters, but it would be best, if I could read in a time series
  double epsilon_deg = 0;
  if (m_use_paleo_file) {
    epsilon_deg = m_obliquity->value(time);
  } else {
    epsilon_deg = m_config->get_number("surface.itm.paleo.obliquity");
  }
  double lambda = get_lambda_paleo(time);
  double delta  = sin(epsilon_deg * M_PI / 180.) * sin(lambda);
  // Equation 2.2.4 of Liou (2002)
  return delta;
}


double DEBMSimple::get_lambda_paleo(double time) {
  // estimates solar longitude at current time in the year
  // Method is using an approximation from :cite:`Berger_1978` section 3 (lambda = 0 at spring equinox).
  // for now the orbital parameters are as config parameters, but it would be best, if I could read in a time series
  double ecc = 0, peri_deg = 0;
  if (m_use_paleo_file) {
    ecc         = m_eccentricity->value(time);
    peri_deg    = m_perihelion_longitude->value(time);
  } else {
    ecc         = m_config->get_number("surface.itm.paleo.eccentricity");
    peri_deg    = m_config->get_number("surface.itm.paleo.long_peri");
  }

  double lambda_m, lambda, delta_lambda;
  delta_lambda = 2. * M_PI * (m_grid->ctx()->time()->year_fraction(time) - 80. / 365.);
  // lambda = 0 at March equinox (80th day of the year)
  double beta     = sqrt(1 - ecc * ecc);
  double peri_rad = peri_deg * M_PI / 180.;
  lambda_m        = -2. * ((ecc / 2. + (pow(ecc, 3)) / 8.) * (1. + beta) * sin(-peri_rad) -
                    (pow(ecc, 2)) / 4. * (1. / 2. + beta) * sin(-2. * peri_rad) +
                    (pow(ecc, 3)) / 8. * (1. / 3. + beta) * sin(-3. * peri_rad)) +
             delta_lambda;
  lambda = (lambda_m + (2. * ecc - (pow(ecc, 3)) / 4.) * sin(lambda_m - peri_rad) +
            (5. / 4.) * (ecc * ecc) * sin(2. * (lambda_m - peri_rad)) +
            (13. / 12.) * (pow(ecc, 3)) * sin(3. * (lambda_m - peri_rad)));
  return lambda;
}


void DEBMSimple::update_impl(const Geometry &geometry, double t, double dt) {

  // update to ensure that temperature and precipitation time series are correct:
  m_atmosphere->update(geometry, t, dt);

  // Use near-surface air temperature as the top-of-the-ice temperature:
  m_temperature->copy_from(m_atmosphere->air_temperature());

  // Set up air temperature and precipitation time series
  int N = m_mbscheme.timeseries_length(dt);

  const double dtseries = dt / N;
  std::vector<double> ts(N), T(N), S(N), P(N), Alb(N);
  ITMMassBalance::Melt ETIM_melt;
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
  const auto &latitude         = geometry.latitude;
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
      &m_tempmelt,
      &m_insolmelt,
      &m_cmelt,
      &m_albedo,
      &m_transmissivity,
      &m_TOAinsol,
      &m_qinsol,
      &latitude,
      &surface_altitude
    };

  if ((bool)m_input_albedo) {
    m_input_albedo->update(t, dt);
    m_input_albedo->init_interpolation(ts);
    list.add(*m_input_albedo.get());
  }

  double
    ice_density    = m_config->get_number("constants.ice.density"),
    sigmalapserate = m_config->get_number("surface.pdd.std_dev_lapse_lat_rate"),
    sigmabaselat   = m_config->get_number("surface.pdd.std_dev_lapse_lat_base");

  bool
    force_albedo = m_config->get_flag("surface.itm.anomaly"),
    paleo        = m_config->get_flag("surface.itm.paleo.enabled");

  m_atmosphere->init_timeseries(ts);
  m_atmosphere->begin_pointwise_access();

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
      if (m_use_air_temp_sd_file) {
        m_air_temp_sd->interp(i, j, S);
      } else {
        double tmp = (*m_air_temp_sd)(i, j);
        for (int k = 0; k < N; ++k) {
          S[k] = tmp;
        }
      }

      if ((bool)m_input_albedo) {
        m_input_albedo->interp(i, j, Alb);
      }

      // apply standard deviation lapse rate on top of prescribed values
      double lat = latitude(i, j);

      if (sigmalapserate != 0.0) {

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
            S[k] = 0.0;
          }
        }
        (*m_air_temp_sd)(i, j) = S[0]; // ensure correct SD reporting
      }

      // Use temperature time series to remove rainfall from precipitation
      m_mbscheme.get_snow_accumulation(T,  // air temperature (input)
                                       P); // precipitation rate (input-output)


      // Use degree-day factors, the number of PDDs, and the snow precipitation to get surface mass
      // balance (and diagnostics: accumulation, melt, runoff)
      {
        double next_snow_depth_reset = m_next_balance_year_start;

        // make copies of firn and snow depth values at this point to avoid accessing 2D
        // fields in the inner loop
        double
          ice        = H(i, j),
          firn       = m_firn_depth(i, j),
          snow       = m_snow_depth(i, j),
          surfelev   = surface_altitude(i, j),
          albedo_loc = m_albedo(i, j);

        double
          A   = 0.0,            // accumulation
          M   = 0.0,            // melt
          R   = 0.0,            // runoff
          SMB = 0.0,            // resulting mass balance
          Mi  = 0.0,            // insolation melt
          Mt  = 0.0,            // temperature melt
          Mc  = 0.0,            // offset melt
          Tr  = 0.0,            // transmissivity, this is just for testing
          Ti  = 0.0,            // top of the atmosphere insolation
          Qi  = 0.0,            // insolation averaged over \Delta t_Phi
          Al  = 0.0;

        // beginning of the loop over small time steps:
        for (int k = 0; k < N; ++k) {

          if (ts[k] >= next_snow_depth_reset) {
            snow = 0.0;
            while (next_snow_depth_reset <= ts[k]) {
              next_snow_depth_reset = m_grid->ctx()->time()->increment_date(next_snow_depth_reset, 1);
            }
          }

          const double accumulation = P[k] * dtseries;

          ITMMassBalance::Changes changes;

          if (force_albedo) {

            if (albedo_anomaly_true(ts[k])) {
              albedo_loc = m_config->get_number("surface.itm.anomaly_value");
            }
          }

          if ((bool)m_input_albedo) {
            albedo_loc = Alb[k];
          }

          double delta = get_delta(ts[k]);
          double distance2 = get_distance2(ts[k]);
          if (paleo) {
            delta     = get_delta_paleo(ts[k]);
            distance2 = get_distance2_paleo(ts[k]);
          }


          ETIM_melt = m_mbscheme.calculate_melt(dtseries, S[k], T[k], surfelev, delta, distance2,
                                                lat * M_PI / 180., albedo_loc);

          //  no melt over ice-free ocean
          if (mask.ice_free_ocean(i, j)) {
            ETIM_melt.T_melt   = 0.;
            ETIM_melt.I_melt   = 0.;
            ETIM_melt.c_melt   = 0.;
            ETIM_melt.ITM_melt = 0.;
          }

          changes = m_mbscheme.step(m_refreeze_fraction, ice, ETIM_melt.ITM_melt, firn, snow, accumulation);

          if (not(bool) m_input_albedo) {
            MaskValue cell_type = static_cast<MaskValue>(mask.as_int(i, j));
            albedo_loc = m_mbscheme.albedo(changes.melt, cell_type, dtseries);
          }
          if (force_albedo) {
            if (albedo_anomaly_true(ts[k])) {
              albedo_loc = m_config->get_number("surface.itm.anomaly_value");
            }
          }

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
            A += accumulation;
            M += changes.melt;
            Mt += ETIM_melt.T_melt;
            Mi += ETIM_melt.I_melt;
            Mc += ETIM_melt.c_melt;
            R += changes.runoff;
            SMB += changes.smb, Tr += ETIM_melt.transmissivity;
            Ti += ETIM_melt.TOA_insol;
            Qi += ETIM_melt.q_insol;
            Al += albedo_loc;
          }
        } // end of the time-stepping loop

        // set firn and snow depths
        m_firn_depth(i, j)     = firn;
        m_snow_depth(i, j)     = snow;
        m_albedo(i, j)         = Al / N;
        m_transmissivity(i, j) = Tr / N; //ETIM_melt.transmissivity;
        m_TOAinsol(i, j)       = Ti / N;
        m_qinsol(i, j)         = Qi / N;

        // set melt terms at this point, converting
        // from "meters, ice equivalent" to "kg / m^2"
        m_tempmelt(i, j)  = Mt * ice_density;
        m_insolmelt(i, j) = Mi * ice_density;
        m_cmelt(i, j)     = Mc * ice_density;

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

const IceModelVec2S &DEBMSimple::surface_insolation_melt() const {
  return m_insolmelt;
}

const IceModelVec2S &DEBMSimple::surface_temperature_melt() const {
  return m_tempmelt;
}

const IceModelVec2S &DEBMSimple::surface_offset_melt() const {
  return m_cmelt;
}

const IceModelVec2S &DEBMSimple::albedo() const {
  return m_albedo;
}

const IceModelVec2S &DEBMSimple::transmissivity() const {
  return m_transmissivity;
}

const IceModelVec2S &DEBMSimple::TOAinsol() const {
  return m_TOAinsol;
}

const IceModelVec2S &DEBMSimple::qinsol() const {
  return m_qinsol;
}

void DEBMSimple::define_model_state_impl(const File &output) const {
  SurfaceModel::define_model_state_impl(output);
  m_firn_depth.define(output, PISM_DOUBLE);
  m_snow_depth.define(output, PISM_DOUBLE);
  m_albedo.define(output, PISM_DOUBLE);
}

void DEBMSimple::write_model_state_impl(const File &output) const {
  SurfaceModel::write_model_state_impl(output);
  m_firn_depth.write(output);
  m_snow_depth.write(output);
  m_albedo.write(output);
}

namespace diagnostics {

enum AmountKind { AMOUNT, MASS };

/*! @brief Report surface insolation melt, averaged over the reporting interval */
class SurfaceInsolationMelt : public DiagAverageRate<DEBMSimple> {
public:
  SurfaceInsolationMelt(const DEBMSimple *m, AmountKind kind)
      : DiagAverageRate<DEBMSimple>(
            m, kind == AMOUNT ? "surface_insolation_melt_flux" : "surface_insolation_melt_rate", TOTAL_CHANGE),
        m_kind(kind),
        m_melt_mass(m_grid, "insolation_melt_mass", WITHOUT_GHOSTS) {

    std::string name          = "surface_insolation_melt_flux",
                long_name     = "surface insolation melt, averaged over the reporting interval",
                standard_name = "surface_insolation_melt_flux", accumulator_units = "kg m-2",
                internal_units = "kg m-2 second-1", external_units = "kg m-2 year-1";
    if (kind == MASS) {
      name          = "surface_insolation_melt_rate";
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
    const IceModelVec2S &melt_amount = model->surface_insolation_melt();

    if (m_kind == MASS) {
      m_melt_mass.copy_from(melt_amount);
      m_melt_mass.scale(m_grid->cell_area());
      return m_melt_mass;
    } else {
      return melt_amount;
    }
  }

private:
  AmountKind m_kind;
  IceModelVec2S m_melt_mass;
};

/*! @brief Report surface temperature melt, averaged over the reporting interval */
class SurfaceTemperatureMelt : public DiagAverageRate<DEBMSimple> {
public:
  SurfaceTemperatureMelt(const DEBMSimple *m, AmountKind kind)
      : DiagAverageRate<DEBMSimple>(
            m, kind == AMOUNT ? "surface_temperature_melt_flux" : "surface_temperature_melt_rate", TOTAL_CHANGE),
        m_kind(kind),
        m_melt_mass(m_grid, "temperature_melt_mass", WITHOUT_GHOSTS) {

    std::string name          = "surface_temperature_melt_flux",
                long_name     = "surface temperature melt, averaged over the reporting interval",
                standard_name = "surface_temperature_melt_flux", accumulator_units = "kg m-2",
                internal_units = "kg m-2 second-1", external_units = "kg m-2 year-1";
    if (kind == MASS) {
      name          = "surface_temperature_melt_rate";
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
    const IceModelVec2S &melt_amount = model->surface_temperature_melt();

    if (m_kind == MASS) {
      m_melt_mass.copy_from(melt_amount);
      m_melt_mass.scale(m_grid->cell_area());
      return m_melt_mass;
    } else {
      return melt_amount;
    }
  }

private:
  AmountKind m_kind;
  IceModelVec2S m_melt_mass;
};

/*! @brief Report surface offset melt, averaged over the reporting interval */
class SurfaceOffsetMelt : public DiagAverageRate<DEBMSimple> {
public:
  SurfaceOffsetMelt(const DEBMSimple *m, AmountKind kind)
      : DiagAverageRate<DEBMSimple>(
            m, kind == AMOUNT ? "surface_offset_melt_flux" : "surface_offset_melt_rate", TOTAL_CHANGE),
        m_kind(kind),
        m_melt_mass(m_grid, "offset_melt_mass", WITHOUT_GHOSTS) {

    std::string name          = "surface_offset_melt_flux",
                long_name     = "surface offset melt, averaged over the reporting interval",
                standard_name = "surface_offset_melt_flux", accumulator_units = "kg m-2",
                internal_units = "kg m-2 second-1", external_units = "kg m-2 year-1";
    if (kind == MASS) {
      name          = "surface_offset_melt_rate";
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
    const IceModelVec2S &melt_amount = model->surface_offset_melt();

    if (m_kind == MASS) {
      m_melt_mass.copy_from(melt_amount);
      m_melt_mass.scale(m_grid->cell_area());
      return m_melt_mass;
    } else {
      return melt_amount;
    }
  }

private:
  AmountKind m_kind;
  IceModelVec2S m_melt_mass;
};

} // end of namespace diagnostics

DiagnosticList DEBMSimple::diagnostics_impl() const {
  using namespace diagnostics;

  DiagnosticList result = {
    { "surface_insolation_melt_flux", Diagnostic::Ptr(new SurfaceInsolationMelt(this, AMOUNT)) },
    { "surface_insolation_melt_rate", Diagnostic::Ptr(new SurfaceInsolationMelt(this, MASS)) },
    { "surface_temperature_melt_flux", Diagnostic::Ptr(new SurfaceTemperatureMelt(this, AMOUNT)) },
    { "surface_temperature_melt_rate", Diagnostic::Ptr(new SurfaceTemperatureMelt(this, MASS)) },
    { "surface_offset_melt_flux", Diagnostic::Ptr(new SurfaceOffsetMelt(this, AMOUNT)) },
    { "surface_offset_melt_rate", Diagnostic::Ptr(new SurfaceOffsetMelt(this, MASS)) },
    { "air_temp_sd", Diagnostic::wrap(*m_air_temp_sd) },
    { "snow_depth", Diagnostic::wrap(m_snow_depth) },
    { "firn_depth", Diagnostic::wrap(m_firn_depth) },
    { "albedo", Diagnostic::wrap(m_albedo) },
    { "transmissivity", Diagnostic::wrap(m_transmissivity) },
    { "TOAinsol", Diagnostic::wrap(m_TOAinsol) },
    { "qinsol", Diagnostic::wrap(m_qinsol) }
  };

  result = pism::combine(result, SurfaceModel::diagnostics_impl());

  return result;
}

} // end of namespace surface
} // end of namespace pism
