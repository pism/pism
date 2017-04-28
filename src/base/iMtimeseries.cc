// Copyright (C) 2009-2017 Constantine Khroulev
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

#include <gsl/gsl_interp.h>     // gsl_interp_bsearch()

#include "iceModel.hh"

#include "base/energy/EnergyModel.hh"

#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMDiagnostic.hh"
#include "base/util/PISMTime.hh"
#include "base/util/error_handling.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/MaxTimestep.hh"
#include "base/util/Profiling.hh"
#include "base/util/pism_utilities.hh"

namespace pism {

//! Initializes the code writing scalar time-series.
void IceModel::init_timeseries() {

  m_ts_filename = m_config->get_string("output.timeseries.filename");

  options::String times("-ts_times", "Specifies a MATLAB-style range or a list of requested times");

  if (times.is_set() xor not m_ts_filename.empty()) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "you need to specity both -ts_file and -ts_times"
                       " to save scalar diagnostic time-series.");
  }

  if (m_ts_filename.empty()) {
    m_ts_diagnostics.clear();
    return;
  }

  try {
    m_time->parse_times(times, *m_ts_times);
  } catch (RuntimeError &e) {
    e.add_context("parsing the -ts_times argument %s", times->c_str());
    throw;
  }

  m_log->message(2, "  saving scalar time-series to '%s'\n", m_ts_filename.c_str());
  m_log->message(2, "  times requested: %s\n", times->c_str());

  std::set<std::string> vars = set_split(m_config->get_string("output.timeseries.variables"), ',');
  if (not vars.empty()) {
    m_log->message(2, "variables requested: %s\n", set_join(vars, ",").c_str());

    std::map<std::string, TSDiagnostic::Ptr> diagnostics;
    for (auto v : vars) {
      if (m_ts_diagnostics.find(v) != m_ts_diagnostics.end()) {
        diagnostics[v] = m_ts_diagnostics[v];
      } else {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "scalar diagnostic '%s' is not available",
                                      v.c_str());
      }
    }
    // replace m_ts_diagnostics with requested diagnostics, de-allocating the rest
    m_ts_diagnostics = diagnostics;
  } else {
    // use all diagnostics in m_ts_diagnostics
  }

  // prepare the output file
  {
    // default behavior is to move the file aside if it exists already; option allows appending
    bool append = m_config->get_boolean("output.timeseries.append");
    IO_Mode mode = append ? PISM_READWRITE : PISM_READWRITE_MOVE;
    PIO file(m_grid->com, "netcdf3", m_ts_filename, mode);      // Use NetCDF-3 to write time-series.
    // add the last saved time to the list of requested times so that the first time is interpreted
    // as the end of a reporting time step
    if (append and file.inq_dimlen("time") > 0) {
      double
        epsilon = 1.0,          // one second
        t       = 0.0;
      file.inq_dim_limits("time", NULL, &t);
      // add this time only if it is strictly before the first requested one
      if (t + epsilon < m_ts_times->front()) {
        m_ts_times->insert(m_ts_times->begin(), t);
      }
    }

    write_metadata(file, SKIP_MAPPING, PREPEND_HISTORY);
    write_run_stats(file);

    // initialize scalar diagnostics
    for (auto d : m_ts_diagnostics) {
      d.second->init(file, m_ts_times);
    }
  }
}

//! Initialize the code saving spatially-variable diagnostic quantities.
void IceModel::init_extras() {

  m_last_extra = 0;               // will be set in write_extras()
  m_next_extra = 0;

  options::String extra_file("-extra_file", "Specifies the output file");
  m_extra_filename = extra_file;

  options::String times("-extra_times", "Specifies times to save at");

  options::StringSet vars("-extra_vars",
                          "Specifies a comma-separated list of variables to save", "");

  bool split  = options::Bool("-extra_split", "Specifies whether to save to separate files");
  bool append = options::Bool("-extra_append", "append spatial diagnostics");

  if (extra_file.is_set() ^ times.is_set()) {
    throw RuntimeError(PISM_ERROR_LOCATION, "you need to specify both -extra_file and -extra_times to save spatial time-series.");
  }

  if (not extra_file.is_set() && not times.is_set()) {
    m_save_extra = false;
    return;
  }

  try {
    m_time->parse_times(times, m_extra_times);
  } catch (RuntimeError &e) {
    e.add_context("parsing the -extra_times argument %s", times->c_str());
    throw;
  }

  if (m_extra_times.size() == 0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "no argument for -extra_times option.");
  }

  if (append and split) {
    throw RuntimeError(PISM_ERROR_LOCATION, "both -extra_split and -extra_append are set.");
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

  m_log->message(2, "times requested: %s\n", times->c_str());

  if (m_extra_times.size() > 500) {
    m_log->message(2,
               "PISM WARNING: more than 500 times requested. This might fill your hard-drive!\n");
  }

  if (vars.is_set()) {
    m_extra_vars = vars;
    m_log->message(2, "variables requested: %s\n", vars.to_string().c_str());
  } else {
    m_log->message(2, "PISM WARNING: -extra_vars was not set. Writing the model state...\n");

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
  bool append = options::Bool("-extra_append", "append -extra_file output");
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

static MaxTimestep reporting_max_timestep(const std::vector<double> &times, double t,
                                          const std::string &description) {

  const size_t N = times.size();
  if (t >= times.back()) {
    return MaxTimestep();
  }

  size_t j = 0;
  double dt = 0.0;
  if (t < times[0]) {
    j = -1;
  } else {
    j = gsl_interp_bsearch(&times[0], t, 0, N - 1);
  }

  dt = times[j + 1] - t;

  // now make sure that we don't end up taking a time-step of less than 1
  // second long
  if (dt < 1.0) {
    if (j + 2 < N) {
      return MaxTimestep(times[j + 2] - t, description);
    } else {
      return MaxTimestep(description);
    }
  } else {
    return MaxTimestep(dt, description);
  }
}

//! Computes the maximum time-step we can take and still hit all `-extra_times`.
MaxTimestep IceModel::extras_max_timestep(double my_t) {

  if ((not m_save_extra) or
      (not m_config->get_boolean("time_stepping.hit_extra_times"))) {
    return MaxTimestep("reporting (-extra_times)");
  }

  return reporting_max_timestep(m_extra_times, my_t, "reporting (-extra_times)");
}

//! Computes the maximum time-step we can take and still hit all `-save_times`.
MaxTimestep IceModel::save_max_timestep(double my_t) {

  if ((not m_save_snapshots) or
      (not m_config->get_boolean("time_stepping.hit_save_times"))) {
    return MaxTimestep("reporting (-save_times)");
  }

  return reporting_max_timestep(m_snapshot_times, my_t, "reporting (-save_times)");
}

//! Computes the maximum time-step we can take and still hit all `-ts_times`.
MaxTimestep IceModel::ts_max_timestep(double my_t) {

  if ((not m_config->get_boolean("time_stepping.hit_ts_times")) or
      m_ts_diagnostics.empty()) {
    return MaxTimestep("reporting (-ts_times)");
  }

  return reporting_max_timestep(*m_ts_times, my_t, "reporting (-ts_times)");
}

//! Flush scalar time-series.
void IceModel::flush_timeseries() {
  // flush all the time-series buffers:
  for (auto d : m_ts_diagnostics) {
    d.second->flush();
  }

  // update run_stats in the time series output file
  if (not m_ts_diagnostics.empty()) {
    PIO file(m_grid->com, "netcdf3", m_ts_filename, PISM_READWRITE);
    write_run_stats(file);
  }
}

} // end of namespace pism
