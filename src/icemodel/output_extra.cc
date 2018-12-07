/* Copyright (C) 2017, 2018 PISM Authors
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

#include "IceModel.hh"

#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Profiling.hh"

namespace pism {

//! Computes the maximum time-step we can take and still hit all `-extra_times`.
MaxTimestep IceModel::extras_max_timestep(double my_t) {

  if ((not m_save_extra) or
      (not m_config->get_boolean("time_stepping.hit_extra_times"))) {
    return MaxTimestep("reporting (-extra_times)");
  }

  return reporting_max_timestep(m_extra_times, my_t, "reporting (-extra_times)");
}

static std::set<std::string> process_extra_shortcuts(const std::set<std::string> &input) {
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

  return result;
}

//! Initialize the code saving spatially-variable diagnostic quantities.
void IceModel::init_extras() {

  m_last_extra = 0;               // will be set in write_extras()
  m_next_extra = 0;

  m_extra_filename   = m_config->get_string("output.extra.file");
  std::string times  = m_config->get_string("output.extra.times");
  std::string vars   = m_config->get_string("output.extra.vars");
  bool        split  = m_config->get_boolean("output.extra.split");
  bool        append = m_config->get_boolean("output.extra.append");

  bool extra_file_set = not m_extra_filename.empty();
  bool times_set = not times.empty();

  if (extra_file_set ^ times_set) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "you need to set both output.extra.file and output.extra.times"
                       " to save spatial time-series.");
  }

  if (not extra_file_set and not times_set) {
    m_save_extra = false;
    return;
  }

  try {
    m_extra_times = m_time->parse_times(times);
  } catch (RuntimeError &e) {
    e.add_context("parsing the output.extra.times argument %s", times.c_str());
    throw;
  }

  if (m_extra_times.size() == 0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "output.extra.times cannot be empty");
  }

  if (append and split) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "both output.extra.split and output.extra.append are set.");
  }

  if (append) {
    PIO file(m_grid->com, m_config->get_string("output.format"), m_extra_filename, PISM_READONLY);

    std::string time_name = m_config->get_string("time.dimension_name");
    if (file.inq_var(time_name)) {
      double time_max;

      file.inq_dim_limits(time_name, NULL, &time_max);

      while (m_next_extra + 1 < m_extra_times.size() && m_extra_times[m_next_extra + 1] < time_max) {
        m_next_extra++;
      }

      if (m_next_extra > 0) {
        m_log->message(2,
                   "skipping times before the last record in %s (at %s)\n",
                   m_extra_filename.c_str(), m_time->date(time_max).c_str());
      }

      // discard requested times before the beginning of the run
      std::vector<double> tmp(m_extra_times.size() - m_next_extra);
      for (unsigned int k = 0; k < tmp.size(); ++k) {
        tmp[k] = m_extra_times[m_next_extra + k];
      }

      m_extra_times = tmp;
      m_next_extra = 0;
    }
    file.close();
  }

  m_save_extra          = true;
  m_extra_file_is_ready = false;
  m_split_extra         = false;

  if (split) {
    m_split_extra = true;
    m_log->message(2, "saving spatial time-series to '%s+year.nc'; ",
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

  if (not vars.empty()) {
    m_extra_vars = process_extra_shortcuts(set_split(vars, ','));
    m_log->message(2, "variables requested: %s\n", vars.c_str());
  } else {
    m_log->message(2,
                   "PISM WARNING: output.extra.vars was not set. Writing the model state...\n");
  } // end of the else clause after "if (extra_vars_set)"
}

//! Write spatially-variable diagnostic quantities.
void IceModel::write_extras() {
  double saving_after = -1.0e30; // initialize to avoid compiler warning; this
                                 // value is never used, because saving_after
                                 // is only used if save_now == true, and in
                                 // this case saving_after is guaranteed to be
                                 // initialized. See the code below.
  char filename[PETSC_MAX_PATH_LEN];
  unsigned int current_extra;
  // determine if the user set the -save_at and -save_to options
  if (not m_save_extra) {
    return;
  }

  double current_time = m_time->current();

  // do we need to save *now*?
  if (m_next_extra < m_extra_times.size() and
      (current_time >= m_extra_times[m_next_extra] or
       fabs(current_time - m_extra_times[m_next_extra]) < 1.0)) {
    // the condition above is "true" if we passed a requested time or got to
    // within 1 second from it

    current_extra = m_next_extra;

    // update next_extra
    while (m_next_extra < m_extra_times.size() and
           (m_extra_times[m_next_extra] <= current_time or
            fabs(current_time - m_extra_times[m_next_extra]) < 1.0)) {
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

    return;
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

  if (m_split_extra) {
    m_extra_file_is_ready = false;        // each time-series record is written to a separate file
    snprintf(filename, PETSC_MAX_PATH_LEN, "%s_%s.nc",
             m_extra_filename.c_str(), m_time->date().c_str());
  } else {
    strncpy(filename, m_extra_filename.c_str(), PETSC_MAX_PATH_LEN);
  }

  m_log->message(3,
                 "saving spatial time-series to %s at %s\n",
                 filename, m_time->date().c_str());

  // default behavior is to move the file aside if it exists already; option allows appending
  bool append = m_config->get_boolean("output.extra.append");
  IO_Mode mode = m_extra_file_is_ready or append ? PISM_READWRITE : PISM_READWRITE_MOVE;

  const Profiling &profiling = m_ctx->profiling();
  profiling.begin("io.extra_file");
  {
    PIO file(m_grid->com, m_config->get_string("output.format"), filename, mode);
    std::string time_name = m_config->get_string("time.dimension_name");

    if (not m_extra_file_is_ready) {
      // Prepare the file:
      io::define_time(file, *m_ctx);
      file.put_att_text(time_name, "bounds", "time_bounds");

      write_metadata(file, WRITE_MAPPING, PREPEND_HISTORY);

      m_extra_file_is_ready = true;
    }

    write_run_stats(file);

    save_variables(file,
                   m_extra_vars.empty() ? INCLUDE_MODEL_STATE : JUST_DIAGNOSTICS,
                   m_extra_vars, PISM_FLOAT);

    // Get the length of the time dimension *after* it is appended to.
    unsigned int time_length = file.inq_dimlen(time_name);
    size_t time_start = time_length > 0 ? static_cast<size_t>(time_length - 1) : 0;

    io::write_time_bounds(file, m_extra_bounds, time_start, {m_last_extra, current_time});
  }
  profiling.end("io.extra_file");

  flush_timeseries();

  m_last_extra = current_time;

  // reset accumulators in diagnostics that compute time averaged quantities
  reset_diagnostics();
}

} // end of namespace pism
