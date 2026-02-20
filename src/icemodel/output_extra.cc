/* Copyright (C) 2017, 2018, 2019, 2020, 2021, 2023, 2024, 2025, 2026 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include <memory>

#include "pism/icemodel/IceModel.hh"

#include "pism/util/pism_utilities.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/io/io_helpers.hh"

namespace pism {

//! Computes the maximum time-step we can take and still hit all `-spatial_times`.
MaxTimestep IceModel::spatial_diagnostics_max_timestep(double t) {

  if (m_spatial_filename.empty() or
      (not m_config->get_flag("time_stepping.hit_spatial_times"))) {
    return MaxTimestep("reporting (-spatial_times)");
  }

  double eps = m_config->get_number("time_stepping.resolution");

  return reporting_max_timestep(m_spatial_times, t, eps,
                                "reporting (-spatial_times)");
}

static std::set<std::string> process_extra_shortcuts(const Config &config,
                                                     const std::set<std::string> &input) {
  std::set<std::string> result = input;

  // process shortcuts
  if (result.find("amount_fluxes") != result.end()) {
    result.erase("amount_fluxes");
    result.insert({ "tendency_of_ice_amount", "tendency_of_ice_amount_due_to_basal_mass_flux",
                    "tendency_of_ice_amount_due_to_conservation_error",
                    "tendency_of_ice_amount_due_to_discharge", "tendency_of_ice_amount_due_to_flow",
                    "tendency_of_ice_amount_due_to_surface_mass_flux" });
  }

  if (result.find("mass_fluxes") != result.end()) {
    result.erase("mass_fluxes");
    result.insert({ "tendency_of_ice_mass", "tendency_of_ice_mass_due_to_basal_mass_flux",
                    "tendency_of_ice_mass_due_to_conservation_error",
                    "tendency_of_ice_mass_due_to_discharge", "tendency_of_ice_mass_due_to_flow",
                    "tendency_of_ice_mass_due_to_surface_mass_flux" });
  }

  if (result.find("pdd_fluxes") != result.end()) {
    result.erase("pdd_fluxes");
    result.insert({ "surface_accumulation_flux", "surface_runoff_flux", "surface_melt_flux" });
  }

  if (result.find("pdd_rates") != result.end()) {
    result.erase("pdd_rates");
    result.insert({ "surface_accumulation_rate", "surface_runoff_rate", "surface_melt_rate" });
  }

  if (result.find("hydrology_fluxes") != result.end()) {
    result.erase("hydrology_fluxes");
    result.insert({ "tendency_of_subglacial_water_mass",
                    "tendency_of_subglacial_water_mass_due_to_input",
                    "tendency_of_subglacial_water_mass_due_to_flow",
                    "tendency_of_subglacial_water_mass_due_to_conservation_error",
                    "tendency_of_subglacial_water_mass_at_grounded_margins",
                    "tendency_of_subglacial_water_mass_at_grounding_line",
                    "tendency_of_subglacial_water_mass_at_domain_boundary" });
  }

  if (result.find("ismip6") != result.end()) {

    const char *flag_name = "output.ISMIP6";

    if (not config.get_flag(flag_name)) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Please set %s to save ISMIP6 diagnostics "
                                    "(-spatial_vars ismip6).", flag_name);
    }

    result.erase("ismip6");
    for (const auto& v : set_split(config.get_string("output.ISMIP6_spatial_variables"), ',')) {
      result.insert(v);
    }
  }

  return result;
}

//! Initialize the code saving spatially-variable diagnostic quantities.
void IceModel::init_spatial_diagnostics() {

  m_last_spatial_time = 0;               // will be set in write_extras()
  m_next_spatial_index = 0;

  m_spatial_filename   = m_config->get_string("output.spatial.file");
  std::string times  = m_config->get_string("output.spatial.times");
  bool        split  = m_config->get_flag("output.spatial.split");
  bool        append = m_config->get_flag("output.spatial.append");

  bool file_set = not m_spatial_filename.empty();
  bool times_set = not times.empty();

  if (file_set ^ times_set) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "you need to set both output.spatial.file and output.spatial.times"
                       " to save spatial time-series.");
  }

  if (not file_set and not times_set) {
    m_spatial_filename.clear();
    return;
  }

  try {
    m_spatial_times = m_time->parse_times(times);
  } catch (RuntimeError &e) {
    e.add_context("parsing the output.spatial.times argument %s", times.c_str());
    throw;
  }

  if (m_spatial_times.empty()) {
    throw RuntimeError(PISM_ERROR_LOCATION, "output.spatial.times cannot be empty");
  }

  if (append and split) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "both output.spatial.split and output.spatial.append are set.");
  }

  // initialize m_extra_vars and m_extra_file_contents
  {
    auto vars = m_config->get_string("output.spatial.vars");
    if (not vars.empty()) {
      m_spatial_vars = process_extra_shortcuts(*m_config, set_split(vars, ','));
      m_log->message(2, "variables requested: %s\n", vars.c_str());
    } else {
      m_log->message(2,
                     "PISM WARNING: output.spatial.vars was not set. Writing the model state...\n");
      m_spatial_vars = {};
    }

    if (m_spatial_vars.empty()) {
      m_spatial_file_contents = state_variables();
    } else {
      m_spatial_file_contents = diagnostic_variables(m_spatial_vars);
    }
    m_spatial_file_contents = pism::combine(m_spatial_file_contents, common_metadata());
  }

  m_spatial_file = nullptr;
  if (not split) {
    m_spatial_file = std::make_shared<OutputFile>(m_extra_writer, m_spatial_filename);

    if (append) {
      // assume that the file is ready to write to, i.e. time and time_bounds are already
      // defined
      m_spatial_file->append();

      if (m_spatial_file->time_dimension_length() > 0) {
        double time_max = m_spatial_file->last_time_value();

        while (m_next_spatial_index + 1 < m_spatial_times.size() &&
               m_spatial_times[m_next_spatial_index + 1] < time_max) {
          m_next_spatial_index++;
        }

        if (m_next_spatial_index > 0) {
          m_log->message(2, "skipping times before the last record in %s (at %s)\n",
                         m_spatial_filename.c_str(), m_time->date(time_max).c_str());
        }

        // discard requested times before the beginning of the run
        std::vector<double> tmp(m_spatial_times.size() - m_next_spatial_index);
        for (unsigned int k = 0; k < tmp.size(); ++k) {
          tmp[k] = m_spatial_times[m_next_spatial_index + k];
        }

        m_spatial_times = tmp;
        m_next_spatial_index  = 0;
      }
    } else {
      // prepare the output file
      bool with_time_bounds = true;
      define_time(*m_spatial_file, with_time_bounds);
      define_variables(*m_spatial_file, m_spatial_file_contents);
    }
  }

  if (split) {
    m_log->message(2, "saving spatial time-series to '%s+date.nc'; ",
               m_spatial_filename.c_str());
  } else {
    if (not ends_with(m_spatial_filename, ".nc")) {
      m_log->message(2,
                 "PISM WARNING: spatial time-series file name '%s' does not have the '.nc' suffix!\n",
                 m_spatial_filename.c_str());
    }
    m_log->message(2, "saving spatial time-series to '%s'; ",
               m_spatial_filename.c_str());
  }

  m_log->message(2, "times requested: %s\n", times.c_str());

  if (m_spatial_times.size() > 500) {
    m_log->message(2,
               "PISM WARNING: more than 500 times requested. This might fill your hard-drive!\n");
  }

  if (pism::netcdf_version() > 0 and pism::netcdf_version() < 473) {
    if (m_spatial_times.size() > 5000 and
        m_config->get_string("output.format") == "netcdf4_parallel") {
      throw RuntimeError(
          PISM_ERROR_LOCATION,
          "more than 5000 times requested."
          "Please use -spatial_split to avoid a crash caused by a bug in NetCDF versions older than 4.7.3.\n"
          "Alternatively\n"
          "- split this simulation into several runs and then concatenate results\n"
          "- select a different output.format value\n"
          "- upgrade NetCDF to 4.7.3");
    }
  }

}

//! Write spatially-variable diagnostic quantities.
void IceModel::write_spatial_diagnostics() {
  double saving_after = -1.0e30; // initialize to avoid compiler warning; this
                                 // value is never used, because saving_after
                                 // is only used if save_now == true, and in
                                 // this case saving_after is guaranteed to be
                                 // initialized. See the code below.
  unsigned int current_extra;
  // determine if the user set the -save_at and -save_to options
  if (m_spatial_filename.empty()) {
    return;
  }

  const double time_resolution = m_config->get_number("time_stepping.resolution");
  double current_time = m_time->current();

  // do we need to save *now*?
  if (m_next_spatial_index < m_spatial_times.size() and
      (current_time >= m_spatial_times[m_next_spatial_index] or
       fabs(current_time - m_spatial_times[m_next_spatial_index]) < time_resolution)) {
    // the condition above is "true" if we passed a requested time or got to
    // within time_resolution seconds from it

    current_extra = m_next_spatial_index;

    // update next_extra
    while (m_next_spatial_index < m_spatial_times.size() and
           (m_spatial_times[m_next_spatial_index] <= current_time or
            fabs(current_time - m_spatial_times[m_next_spatial_index]) < time_resolution)) {
      m_next_spatial_index++;
    }

    saving_after = m_spatial_times[current_extra];
  } else {
    return;
  }

  if (current_extra == 0) {
    // The first time defines the left end-point of the first reporting interval; we don't write a
    // report at this time.

    // Re-initialize last_extra (the correct value is not known at the time init_extras() is
    // called).
    m_last_spatial_time = current_time;

    // ISMIP6 runs need to save diagnostics at the beginning of the run
    if (not m_config->get_flag("output.ISMIP6")) {
      return;
    }
  }

  if (saving_after < m_time->start()) {
    // Suppose a user tells PISM to write data at times 0:1000:10000. Suppose
    // also that PISM writes a backup file at year 2500 and gets stopped.
    //
    // When restarted, PISM will decide that it's time to write data for time
    // 2000, but
    // * that record was written already and
    // * PISM will end up writing at year 2500, producing a file containing one
    //   more record than necessary.
    //
    // This check makes sure that this never happens.
    return;
  }

  bool split = m_config->get_flag("output.spatial.split");

  const Profiling &profiling = m_ctx->profiling();
  profiling.begin("io.spatial_file");
  {
    if (m_spatial_file == nullptr) {

      std::string filename = m_spatial_filename;
      if (split) {
        // each time-series record is written to a separate file
        auto date_without_spaces = replace_character(m_time->date(m_time->current()), ' ', '_');
        filename = pism::printf("%s_%s.nc", m_spatial_filename.c_str(), date_without_spaces.c_str());
      }

      m_spatial_file.reset(new OutputFile(m_extra_writer, filename));

      if (m_config->get_flag("output.spatial.append")) {
        m_spatial_file->append();
      } else {
        // prepare the output file
        bool with_time_bounds = true;
        define_time(*m_spatial_file, with_time_bounds);
        define_variables(*m_spatial_file, m_spatial_file_contents);
      }
    }

    m_log->message(3, "saving spatial time-series to %s at %s\n", m_spatial_file->name().c_str(),
                   m_time->date(m_time->current()).c_str());

    {
      io::write_config(*m_config, "pism_config", *m_spatial_file);

      // use the mid-point of the current reporting interval
      double time = 0.5 * (m_last_spatial_time + current_time);
      m_spatial_file->append_time(time);

      if (m_spatial_vars.empty()) {
        write_state(*m_spatial_file);
      } else {
        write_diagnostics(*m_spatial_file, m_spatial_vars);
      }

      // write time bounds
      {
        // Get the length of the time dimension *after* it is appended to.
        auto time_length = m_spatial_file->time_dimension_length();
        auto time_start  = time_length > 0 ? (time_length - 1) : 0;

        auto bounds_name = m_time->variable_name() + "_bounds";

        m_spatial_file->write_array(bounds_name, { time_start, 0 }, { 1, 2 },
                                  { m_last_spatial_time, current_time });
      }

      write_run_stats(*m_spatial_file);
    }

    // FIXME: evaluate if we need to "sync" extra files. We should probably do this every
    // Nth (e.g. 10th) time we write to an extra file. This would accomplish most of what
    // "sync()" is for, but at a lower cost.
    //
    // make sure all changes are written
    // m_extra_file->sync();
  }
  profiling.end("io.spatial_file");

  scalar_diagnostics_flush_buffers();

  if (split) {
    // each record is saved to a new file, so we can close this one
    m_spatial_file->close();
    m_spatial_file = nullptr;
  }

  m_last_spatial_time = current_time;

  // reset accumulators in diagnostics that compute time averaged quantities
  for (auto &d : m_available_spatial_diagnostics) {
    d.second->reset();
  }
}

} // end of namespace pism
