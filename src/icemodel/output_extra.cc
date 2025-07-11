/* Copyright (C) 2017, 2018, 2019, 2020, 2021, 2023, 2024, 2025 PISM Authors
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
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/io/io_helpers.hh"

namespace pism {

//! Computes the maximum time-step we can take and still hit all `-extra_times`.
MaxTimestep IceModel::extras_max_timestep(double t) {

  if (m_extra_filename.empty() or
      (not m_config->get_flag("time_stepping.hit_extra_times"))) {
    return MaxTimestep("reporting (-extra_times)");
  }

  double eps = m_config->get_number("time_stepping.resolution");

  return reporting_max_timestep(m_extra_times, t, eps,
                                "reporting (-extra_times)");
}

static std::set<std::string> process_extra_shortcuts(const Config &config,
                                                     const std::set<std::string> &input) {
  std::set<std::string> result = input;

  // process shortcuts
  if (result.find("amount_fluxes") != result.end()) {
    result.erase("amount_fluxes");
    result.insert("tendency_of_ice_amount");
    result.insert("tendency_of_ice_amount_due_to_basal_mass_flux");
    result.insert("tendency_of_ice_amount_due_to_conservation_error");
    result.insert("tendency_of_ice_amount_due_to_discharge");
    result.insert("tendency_of_ice_amount_due_to_flow");
    result.insert("tendency_of_ice_amount_due_to_surface_mass_flux");
  }

  if (result.find("mass_fluxes") != result.end()) {
    result.erase("mass_fluxes");
    result.insert("tendency_of_ice_mass");
    result.insert("tendency_of_ice_mass_due_to_basal_mass_flux");
    result.insert("tendency_of_ice_mass_due_to_conservation_error");
    result.insert("tendency_of_ice_mass_due_to_discharge");
    result.insert("tendency_of_ice_mass_due_to_flow");
    result.insert("tendency_of_ice_mass_due_to_surface_mass_flux");
  }

  if (result.find("pdd_fluxes") != result.end()) {
    result.erase("pdd_fluxes");
    result.insert("surface_accumulation_flux");
    result.insert("surface_runoff_flux");
    result.insert("surface_melt_flux");
  }

  if (result.find("pdd_rates") != result.end()) {
    result.erase("pdd_rates");
    result.insert("surface_accumulation_rate");
    result.insert("surface_runoff_rate");
    result.insert("surface_melt_rate");
  }

  if (result.find("hydrology_fluxes") != result.end()) {
    result.erase("hydrology_fluxes");
    result.insert("tendency_of_subglacial_water_mass");
    result.insert("tendency_of_subglacial_water_mass_due_to_input");
    result.insert("tendency_of_subglacial_water_mass_due_to_flow");
    result.insert("tendency_of_subglacial_water_mass_due_to_conservation_error");
    result.insert("tendency_of_subglacial_water_mass_at_grounded_margins");
    result.insert("tendency_of_subglacial_water_mass_at_grounding_line");
    result.insert("tendency_of_subglacial_water_mass_at_domain_boundary");
  }

  if (result.find("ismip6") != result.end()) {

    const char *flag_name = "output.ISMIP6";

    if (not config.get_flag(flag_name)) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Please set %s to save ISMIP6 diagnostics "
                                    "(-extra_vars ismip6).", flag_name);
    }

    result.erase("ismip6");
    for (const auto& v : set_split(config.get_string("output.ISMIP6_extra_variables"), ',')) {
      result.insert(v);
    }
  }

  return result;
}

//! Initialize the code saving spatially-variable diagnostic quantities.
void IceModel::init_extras() {

  m_last_extra = 0;               // will be set in write_extras()
  m_next_extra = 0;

  m_extra_filename   = m_config->get_string("output.extra.file");
  std::string times  = m_config->get_string("output.extra.times");
  std::string vars   = m_config->get_string("output.extra.vars");
  bool        split  = m_config->get_flag("output.extra.split");
  bool        append = m_config->get_flag("output.extra.append");

  bool extra_file_set = not m_extra_filename.empty();
  bool times_set = not times.empty();

  if (extra_file_set ^ times_set) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "you need to set both output.extra.file and output.extra.times"
                       " to save spatial time-series.");
  }

  if (not extra_file_set and not times_set) {
    m_extra_filename.clear();
    return;
  }

  try {
    m_extra_times = m_time->parse_times(times);
  } catch (RuntimeError &e) {
    e.add_context("parsing the output.extra.times argument %s", times.c_str());
    throw;
  }

  if (m_extra_times.empty()) {
    throw RuntimeError(PISM_ERROR_LOCATION, "output.extra.times cannot be empty");
  }

  if (append and split) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "both output.extra.split and output.extra.append are set.");
  }

  m_extra_file = nullptr;
  if (not split) {
    m_extra_file = std::make_shared<OutputFile>(m_output_writer, m_extra_filename);

    if (append) {
      // assume that the file is ready to write to, i.e. time and time_bounds are already
      // defined
      m_extra_file->append();

      if (m_extra_file->time_dimension_length() > 0) {
        double time_max = m_extra_file->last_time_value();

        while (m_next_extra + 1 < m_extra_times.size() &&
               m_extra_times[m_next_extra + 1] < time_max) {
          m_next_extra++;
        }

        if (m_next_extra > 0) {
          m_log->message(2, "skipping times before the last record in %s (at %s)\n",
                         m_extra_filename.c_str(), m_time->date(time_max).c_str());
        }

        // discard requested times before the beginning of the run
        std::vector<double> tmp(m_extra_times.size() - m_next_extra);
        for (unsigned int k = 0; k < tmp.size(); ++k) {
          tmp[k] = m_extra_times[m_next_extra + k];
        }

        m_extra_times = tmp;
        m_next_extra  = 0;
      }
    } else {
      // prepare the file
      bool with_bounds = true;
      io::define_time_dimension(*m_extra_file, m_time->metadata(), with_bounds);
      define_metadata(*m_extra_file, WRITE_MAPPING);
      define_run_stats(*m_extra_file);
    }
  }

  if (split) {
    m_log->message(2, "saving spatial time-series to '%s+date.nc'; ",
               m_extra_filename.c_str());
  } else {
    if (not ends_with(m_extra_filename, ".nc")) {
      m_log->message(2,
                 "PISM WARNING: spatial time-series file name '%s' does not have the '.nc' suffix!\n",
                 m_extra_filename.c_str());
    }
    m_log->message(2, "saving spatial time-series to '%s'; ",
               m_extra_filename.c_str());
  }

  m_log->message(2, "times requested: %s\n", times.c_str());

  if (m_extra_times.size() > 500) {
    m_log->message(2,
               "PISM WARNING: more than 500 times requested. This might fill your hard-drive!\n");
  }

  if (pism::netcdf_version() > 0 and pism::netcdf_version() < 473) {
    if (m_extra_times.size() > 5000 and
        m_config->get_string("output.format") == "netcdf4_parallel") {
      throw RuntimeError(
          PISM_ERROR_LOCATION,
          "more than 5000 times requested."
          "Please use -extra_split to avoid a crash caused by a bug in NetCDF versions older than 4.7.3.\n"
          "Alternatively\n"
          "- split this simulation into several runs and then concatenate results\n"
          "- select a different output.format value\n"
          "- upgrade NetCDF to 4.7.3");
    }
  }

  if (not vars.empty()) {
    m_extra_vars = process_extra_shortcuts(*m_config, set_split(vars, ','));
    m_log->message(2, "variables requested: %s\n", vars.c_str());
  } else {
    m_log->message(2,
                   "PISM WARNING: output.extra.vars was not set. Writing the model state...\n");
  }
}

//! Write spatially-variable diagnostic quantities.
void IceModel::write_extras() {
  double saving_after = -1.0e30; // initialize to avoid compiler warning; this
                                 // value is never used, because saving_after
                                 // is only used if save_now == true, and in
                                 // this case saving_after is guaranteed to be
                                 // initialized. See the code below.
  unsigned int current_extra;
  // determine if the user set the -save_at and -save_to options
  if (m_extra_filename.empty()) {
    return;
  }

  const double time_resolution = m_config->get_number("time_stepping.resolution");
  double current_time = m_time->current();

  // do we need to save *now*?
  if (m_next_extra < m_extra_times.size() and
      (current_time >= m_extra_times[m_next_extra] or
       fabs(current_time - m_extra_times[m_next_extra]) < time_resolution)) {
    // the condition above is "true" if we passed a requested time or got to
    // within time_resolution seconds from it

    current_extra = m_next_extra;

    // update next_extra
    while (m_next_extra < m_extra_times.size() and
           (m_extra_times[m_next_extra] <= current_time or
            fabs(current_time - m_extra_times[m_next_extra]) < time_resolution)) {
      m_next_extra++;
    }

    saving_after = m_extra_times[current_extra];
  } else {
    return;
  }

  if (current_extra == 0) {
    // The first time defines the left end-point of the first reporting interval; we don't write a
    // report at this time.

    // Re-initialize last_extra (the correct value is not known at the time init_extras() is
    // called).
    m_last_extra = current_time;

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

  bool extra_split = m_config->get_flag("output.extra.split");

  const Profiling &profiling = m_ctx->profiling();
  profiling.begin("io.extra_file");
  {

    std::string time_name = m_time->variable_name();

    VariableMetadata time_bounds("time_bounds", m_sys);

    if (m_extra_file == nullptr) {

      std::string filename = m_extra_filename;
      if (extra_split) {
        // each time-series record is written to a separate file
        auto date_without_spaces = replace_character(m_time->date(m_time->current()), ' ', '_');
        filename = pism::printf("%s_%s.nc", m_extra_filename.c_str(), date_without_spaces.c_str());
      }

      m_extra_file.reset(new OutputFile(m_output_writer, filename));

      if (m_config->get_flag("output.extra.append")) {
        m_extra_file->append();
      } else {
        // Prepare the file:
        bool with_bounds = true;
        io::define_time_dimension(*m_extra_file, m_time->metadata(), with_bounds);
        define_metadata(*m_extra_file, WRITE_MAPPING);
        define_run_stats(*m_extra_file);
      }
    }

    m_log->message(3, "saving spatial time-series to %s at %s\n", m_extra_file->name().c_str(),
                   m_time->date(m_time->current()).c_str());

    write_metadata(*m_extra_file);
    // use the mid-point of the current reporting interval
    double time = 0.5 * (m_last_extra + current_time);
    write_variables(*m_extra_file, m_extra_vars.empty() ? INCLUDE_MODEL_STATE : JUST_DIAGNOSTICS,
                    m_extra_vars, time);

    // Get the length of the time dimension *after* it is appended to.
    auto time_length = m_extra_file->time_dimension_length();
    auto time_start = time_length > 0 ? (time_length - 1) : 0;

    // write time bounds
    m_extra_file->write_array(time_bounds, { time_start, 0 }, { 1, 2 },
                              { m_last_extra, current_time });
    // make sure all changes are written
    m_extra_file->sync();
  }
  profiling.end("io.extra_file");

  flush_timeseries();

  if (extra_split) {
    // each record is saved to a new file, so we can close this one
    m_extra_file->close();
    m_extra_file = nullptr;
  }

  m_last_extra = current_time;

  // reset accumulators in diagnostics that compute time averaged quantities
  reset_diagnostics();
}

} // end of namespace pism
