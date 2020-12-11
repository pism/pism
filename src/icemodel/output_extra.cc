/* Copyright (C) 2017, 2018, 2019, 2020, 2021 PISM Authors
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

#include <netcdf.h>
#ifdef NC_HAVE_META_H
#include <netcdf_meta.h>
#endif

#include "IceModel.hh"

#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Profiling.hh"

#include <mpi.h>
extern "C"{
#include "cdipio.h"
#include "cdi.h"
#include "yaxt.h"
}

namespace pism {

//! Computes the maximum time-step we can take and still hit all `-extra_times`.
MaxTimestep IceModel::extras_max_timestep(double my_t) {

  if ((not m_save_extra) or
      (not m_config->get_flag("time_stepping.hit_extra_times"))) {
    return MaxTimestep("reporting (-extra_times)");
  }

  return reporting_max_timestep(m_extra_times, my_t, "reporting (-extra_times)");
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
    for (auto v : set_split(config.get_string("output.ISMIP6_extra_variables"), ',')) {
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
    File file(m_grid->com, m_extra_filename, PISM_NETCDF3, PISM_READONLY);

    std::string time_name = m_config->get_string("time.dimension_name");
    if (file.find_variable(time_name)) {
      double time_max = vector_max(file.read_dimension(time_name));

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

#ifdef NC_HAVE_META_H
  {
    if (100 * NC_VERSION_MAJOR + 10 * NC_VERSION_MINOR + NC_VERSION_PATCH < 473) {
      if (m_extra_times.size() > 5000 and m_config->get_string("output.format") == "netcdf4_parallel") {
        throw RuntimeError(PISM_ERROR_LOCATION,
                           "more than 5000 times requested."
                           "Please use -extra_split to avoid a crash caused by a bug in NetCDF versions older than 4.7.3.\n"
                           "Alternatively\n"
                           "- split this simulation into several runs and then concatenate results\n"
                           "- select a different output.format value\n"
                           "- upgrade NetCDF to 4.7.3");
      }
    }
  }
#endif

  if (not vars.empty()) {
    m_extra_vars = process_extra_shortcuts(*m_config, set_split(vars, ','));
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
  std::string filename;
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

  int fileID = -1;
  if (m_split_extra) {
    m_extra_file_is_ready = false;        // each time-series record is written to a separate file
    filename = pism::printf("%s_%s.nc",
                            m_extra_filename.c_str(),
                            m_time->date().c_str());
  } else {
    filename = m_extra_filename;
    if (streamIDs.count(filename) > 0) {
      fileID = streamIDs[filename];
    }
  }

  m_log->message(3,
                 "saving spatial time-series to %s at %s\n",
                 filename.c_str(), m_time->date().c_str());

  // default behavior is to move the file aside if it exists already; option allows appending
  bool append = m_config->get_flag("output.extra.append");
  IO_Mode mode = m_extra_file_is_ready or append ? PISM_READWRITE : PISM_READWRITE_MOVE;

  const Profiling &profiling = m_ctx->profiling();
  profiling.begin("io.extra_file");
  {
    if (not m_extra_file) {
      m_extra_file.reset(new File(m_grid->com,
                                  filename,
                                  string_to_backend(m_config->get_string("output.format")),
//                                  string_to_backend("pio_pnetcdf"),
                                  mode,
                                  m_ctx->pio_iosys_id(), ExtraMap, gridIDs, fileID));
    }

    std::string time_name = m_config->get_string("time.dimension_name");

    if (not m_extra_file_is_ready) {
      // Prepare the file:
      io::define_time(*m_extra_file, *m_ctx);
      m_extra_file->write_attribute(time_name, "bounds", "time_bounds");

      io::define_time_bounds(m_extra_bounds,
                             time_name, "nv", *m_extra_file);

      write_metadata(*m_extra_file, WRITE_MAPPING, PREPEND_HISTORY);

      m_extra_file_is_ready = true;
    }

    write_run_stats(*m_extra_file);

    save_variables(*m_extra_file,
                   m_extra_vars.empty() ? INCLUDE_MODEL_STATE : JUST_DIAGNOSTICS,
                   m_extra_vars,
                   0.5 * (m_last_extra + current_time), // use the mid-point of the
                                                        // current reporting interval
                   PISM_FLOAT);

    // Get the length of the time dimension *after* it is appended to.
    unsigned int time_length = m_extra_file->dimension_length(time_name);
    size_t time_start = time_length > 0 ? static_cast<size_t>(time_length - 1) : 0;

    io::write_time_bounds(*m_extra_file, m_extra_bounds,
                          time_start, {m_last_extra, current_time});
    // make sure all changes are written
    m_extra_file->sync();
    if (m_extra_file->backend() == PISM_CDI) {
      streamIDs[filename] = m_extra_file->get_streamID();
      vlistIDs[filename] = m_extra_file->get_vlistID();
      if (gridIDs.size()==0) gridIDs = m_extra_file->get_gridIDs();
      ExtraMap = m_extra_file->get_variables_map();
    }
  }
  if (current_extra < m_extra_times.size()-1) m_sthwritten = true;
  profiling.end("io.extra_file");

  flush_timeseries();

  if (m_split_extra) {
    // each record is saved to a new file, so we can close this one
    m_extra_file.reset(nullptr);
  }

  m_last_extra = current_time;

  // reset accumulators in diagnostics that compute time averaged quantities
  reset_diagnostics();
}

} // end of namespace pism
