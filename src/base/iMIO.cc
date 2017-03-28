// Copyright (C) 2004-2017 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cstring>              // strncpy
#include <cstdio>               // snprintf

#include <algorithm>
#include <set>

#include "iceModel.hh"

#include "base/util/IceGrid.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMDiagnostic.hh"
#include "base/util/PISMTime.hh"
#include "base/util/error_handling.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"

#include "base/util/PISMVars.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/Profiling.hh"
#include "base/util/pism_utilities.hh"
#include "base/util/projection.hh"
#include "base/util/PISMComponent.hh"
#include "base/energy/utilities.hh"

namespace pism {

//! Save model state in NetCDF format.
/*!
Calls save_variables() to do the actual work.
 */
void IceModel::save_results() {
  {
    update_run_stats();

    // build and put string into global attribute "history"
    char str[TEMPORARY_STRING_LENGTH];

    snprintf(str, TEMPORARY_STRING_LENGTH,
             "PISM done.  Performance stats: %.4f wall clock hours, %.4f proc.-hours, %.4f model years per proc.-hour, PETSc MFlops = %.2f.",
             m_run_stats.get_double("wall_clock_hours"),
             m_run_stats.get_double("processor_hours"),
             m_run_stats.get_double("model_years_per_processor_hour"),
             m_run_stats.get_double("PETSc_MFlops"));

    prepend_history(str);

  }

  std::string filename = m_config->get_string("output.file_name");

  if (filename.empty()) {
    m_log->message(2, "WARNING: output.file_name is empty. Using unnamed.nc instead.\n");
    filename = "unnamed.nc";
  }

  if (not ends_with(filename, ".nc")) {
    m_log->message(2,
                   "PISM WARNING: output file name does not have the '.nc' suffix!\n");
  }

  if (m_config->get_string("output.size") != "none") {
    m_log->message(2, "Writing model state to file `%s'\n", filename.c_str());
    save_variables(filename, m_output_vars);
  }
}


void IceModel::write_mapping(const PIO &file) {
  // only write mapping if it is set.
  const VariableMetadata &mapping = m_grid->get_mapping_info().mapping;
  if (mapping.has_attributes()) {
    if (not file.inq_var(mapping.get_name())) {
      file.redef();
      file.def_var(mapping.get_name(), PISM_DOUBLE, std::vector<std::string>());
    }
    io::write_attributes(file, mapping, PISM_DOUBLE, false);
  }
}

void IceModel::write_run_stats(const PIO &file) {
  update_run_stats();
  if (not file.inq_var(m_run_stats.get_name())) {
    file.redef();
    file.def_var(m_run_stats.get_name(), PISM_DOUBLE,
               std::vector<std::string>());
  }
  io::write_attributes(file, m_run_stats, PISM_DOUBLE, false);
}

void IceModel::write_global_attributes(const PIO &file) {
  io::write_global_attributes(file, m_output_global_attributes);
}

void IceModel::write_config(const PIO &file) {
  m_config->write(file);
}

void IceModel::save_variables(const std::string &filename,
                          const std::set<std::string> &variables) {
  const Profiling &profiling = m_ctx->profiling();
  profiling.begin("model state dump");

  {
    PIO nc(m_grid->com, m_config->get_string("output.format"), filename, PISM_READWRITE_MOVE);

    // Write metadata *before* everything else:
    write_mapping(nc);
    write_run_stats(nc);
    write_global_attributes(nc);
    write_config(nc);

    std::string time_name = m_config->get_string("time.dimension_name");
    io::define_time(nc, time_name, m_time->calendar(),
                    m_time->CF_units_string(),
                    m_sys);

    io::append_time(nc, time_name, m_time->current());

    define_model_state(nc);
    define_diagnostics(nc, variables, PISM_DOUBLE);

    write_model_state(nc);
    write_diagnostics(nc, variables);
  }

  profiling.end("model state dump");
}

void IceModel::define_diagnostics(const PIO &file, const std::set<std::string> &variables,
                                  IO_Type nctype) {
  // Define all the variables:
  for (auto var : variables) {
    if (m_grid->variables().is_available(var)) {
      const IceModelVec *v = m_grid->variables().get(var);
      // It has dedicated storage.
      if (var == "mask") {
        v->define(file, PISM_BYTE); // use the default data type
      } else {
        v->define(file, nctype);
      }
    } else {
      // It might be a diagnostic quantity
      Diagnostic::Ptr diag = m_diagnostics[var];

      if (diag) {
        diag->define(file);
      }
    }
  }
}

//! \brief Writes variables listed in vars to filename, using nctype to write
//! fields stored in dedicated IceModelVecs.
void IceModel::write_diagnostics(const PIO &file, const std::set<std::string> &variables) {
  for (auto var : variables) {
    if (m_grid->variables().is_available(var)) {
      m_grid->variables().get(var)->write(file);
    } else {
      Diagnostic::Ptr diag = m_diagnostics[var];

      if (diag) {
        IceModelVec::Ptr v_diagnostic = diag->compute();

        v_diagnostic->write_in_glaciological_units = true;
        v_diagnostic->write(file);
      }
    }
  }
}

void IceModel::define_model_state(const PIO &file) {
  std::set<std::string> variables = output_variables("small");

  // define
  define_diagnostics(file, variables, PISM_DOUBLE);

  for (auto m : m_submodels) {
    m.second->define_model_state(file);
  }
}

void IceModel::write_model_state(const PIO &file) {
  std::set<std::string> variables = output_variables("small");

  write_diagnostics(file, variables);

  for (auto m : m_submodels) {
    m.second->write_model_state(file);
  }
}


//! Manage regridding based on user options.
/*!
  For each variable selected by option `-regrid_vars`, we regrid it onto the current grid from
  the NetCDF file specified by `-regrid_file`.

  This `dimensions` argument can be 2 (regrid 2D variables only), 3 (3D
  only) and 0 (everything).
*/
void IceModel::regrid(int dimensions) {

  if (not (dimensions == 0 ||
           dimensions == 2 ||
           dimensions == 3)) {
    throw RuntimeError(PISM_ERROR_LOCATION, "dimensions can only be 0 (all), 2 or 3");
  }

  options::String regrid_filename("-regrid_file", "Specifies the file to regrid from");

  options::StringSet regrid_vars("-regrid_vars",
                                 "Specifies the list of variables to regrid",
                                 "");


  // Return if no regridding is requested:
  if (not regrid_filename.is_set()) {
     return;
  }

  if (dimensions != 0) {
    m_log->message(2, "regridding %dD variables from file %s ...\n",
               dimensions, regrid_filename->c_str());
  } else {
    m_log->message(2, "regridding from file %s ...\n",regrid_filename->c_str());
  }

  {
    PIO regrid_file(m_grid->com, "guess_mode", regrid_filename, PISM_READONLY);

    if (dimensions == 0) {
      regrid_variables(regrid_file, regrid_vars, 2);
      regrid_variables(regrid_file, regrid_vars, 3);
    } else {
      regrid_variables(regrid_file, regrid_vars, dimensions);
    }
  }
}

void IceModel::regrid_variables(const PIO &regrid_file, const std::set<std::string> &vars,
                                unsigned int ndims) {

  for (auto var : vars) {

    if (not m_grid->variables().is_available(var)) {
      continue;
    }

    // FIXME: remove const_cast. This is bad.
    IceModelVec *v = const_cast<IceModelVec*>(m_grid->variables().get(var));
    SpatialVariableMetadata &m = v->metadata();

    if (v->get_ndims() != ndims) {
      continue;
    }

    std::string pism_intent = m.get_string("pism_intent");
    if (pism_intent != "model_state") {
      m_log->message(2, "  WARNING: skipping '%s' (only model_state variables can be regridded)...\n",
                     var.c_str());
      continue;
    }

    v->regrid(regrid_file, CRITICAL);
  }

  // Check the range of the ice thickness.
  {
    double max_thickness = m_geometry.ice_thickness.range().max,
      Lz = m_grid->Lz();

    if (max_thickness >= Lz + 1e-6) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Maximum ice thickness (%f meters)\n"
                                    "exceeds the height of the computational domain (%f meters).",
                                    max_thickness, Lz);
    }
  }
}

//! Initializes the snapshot-saving mechanism.
void IceModel::init_snapshots() {
  m_current_snapshot = 0;

  options::String save_file("-save_file", "Specifies a snapshot filename",
                            m_snapshots_filename);
  m_snapshots_filename = save_file;

  options::String save_times("-save_times",
                             "Gives a list or a MATLAB-style range of times to save snapshots at");

  bool split = options::Bool("-save_split", "Specifies whether to save snapshots to separate files");

  m_snapshot_vars = output_variables(m_config->get_string("output.save_size"));

  if (save_file.is_set() ^ save_times.is_set()) {
    throw RuntimeError(PISM_ERROR_LOCATION, "you need to specify both -save_file and -save_times to save snapshots.");
  }

  if (!save_file.is_set() && !save_times.is_set()) {
    m_save_snapshots = false;
    return;
  }

  try {
    m_time->parse_times(save_times, m_snapshot_times);
  } catch (RuntimeError &e) {
    e.add_context("parsing the -save_times argument");
    throw;
  }

  if (m_snapshot_times.size() == 0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "no argument for -save_times option.");
  }

  m_save_snapshots = true;
  m_snapshots_file_is_ready = false;
  m_split_snapshots = false;

  if (split) {
    m_split_snapshots = true;
  } else if (!ends_with(m_snapshots_filename, ".nc")) {
    m_log->message(2,
               "PISM WARNING: snapshots file name does not have the '.nc' suffix!\n");
  }

  if (split) {
    m_log->message(2, "saving snapshots to '%s+year.nc'; ",
               m_snapshots_filename.c_str());
  } else {
    m_log->message(2, "saving snapshots to '%s'; ",
               m_snapshots_filename.c_str());
  }

  m_log->message(2, "times requested: %s\n", save_times->c_str());
}

  //! Writes a snapshot of the model state (if necessary)
void IceModel::write_snapshot() {
  double saving_after = -1.0e30; // initialize to avoid compiler warning; this
  // value is never used, because saving_after
  // is only used if save_now == true, and in
  // this case saving_after is guaranteed to be
  // initialized. See the code below.
  char filename[PETSC_MAX_PATH_LEN];

  // determine if the user set the -save_times and -save_file options
  if (!m_save_snapshots) {
    return;
  }

  // do we need to save *now*?
  if ((m_time->current() >= m_snapshot_times[m_current_snapshot]) and
      (m_current_snapshot < m_snapshot_times.size())) {
    saving_after = m_snapshot_times[m_current_snapshot];

    while ((m_current_snapshot < m_snapshot_times.size()) and
           (m_snapshot_times[m_current_snapshot] <= m_time->current())) {
      m_current_snapshot++;
    }
  } else {
    // we don't need to save now, so just return
    return;
  }

  // flush time-series buffers
  flush_timeseries();

  if (m_split_snapshots) {
    m_snapshots_file_is_ready = false;    // each snapshot is written to a separate file
    snprintf(filename, PETSC_MAX_PATH_LEN, "%s_%s.nc",
             m_snapshots_filename.c_str(), m_time->date(saving_after).c_str());
  } else {
    strncpy(filename, m_snapshots_filename.c_str(), PETSC_MAX_PATH_LEN);
  }

  m_log->message(2,
             "\nsaving snapshot to %s at %s, for time-step goal %s\n\n",
             filename, m_time->date().c_str(),
             m_time->date(saving_after).c_str());

  IO_Mode mode = m_snapshots_file_is_ready ? PISM_READWRITE : PISM_READWRITE_MOVE;
  PIO nc(m_grid->com, m_config->get_string("output.format"), filename, mode);

  if (not m_snapshots_file_is_ready) {
    // Prepare the snapshots file:
    io::define_time(nc, m_config->get_string("time.dimension_name"),
                    m_time->calendar(),
                    m_time->CF_units_string(),
                    m_sys);

    m_snapshots_file_is_ready = true;
  }

  // write metadata to the file *every time* we update it
  write_mapping(nc);
  write_run_stats(nc);
  write_global_attributes(nc);
  write_config(nc);

  io::append_time(nc, m_config->get_string("time.dimension_name"), m_time->current());

  define_model_state(nc);
  define_diagnostics(nc, m_snapshot_vars, PISM_DOUBLE);

  write_model_state(nc);
  write_diagnostics(nc, m_snapshot_vars);

  {
    // find out how much time passed since the beginning of the run
    double wall_clock_hours = pism::wall_clock_hours(m_grid->com, m_start_time);

    // Get time length now, i.e. after writing variables. This forces PISM to call PIO::enddef(), so
    // that the length of the time dimension is up to date.
    unsigned int time_length = nc.inq_dimlen(m_config->get_string("time.dimension_name"));

    // make sure that time_start is valid even if time_length is zero
    size_t time_start = 0;
    if (time_length > 0) {
      time_start = static_cast<size_t>(time_length - 1);
    } else {
      time_start = 0;
    }

    io::write_timeseries(nc, m_timestamp, time_start, wall_clock_hours);
  }

  nc.close();
}

//! Initialize the backup (snapshot-on-wallclock-time) mechanism.
void IceModel::init_backups() {

  m_backup_interval = m_config->get_double("output.backup_interval");

  std::string backup_file = m_config->get_string("output.file_name");
  if (not backup_file.empty()) {
    m_backup_filename = pism_filename_add_suffix(backup_file, "_backup", "");
  } else {
    m_backup_filename = "pism_backup.nc";
  }

  m_backup_interval = options::Real("-backup_interval",
                                  "Automatic backup interval, hours", m_backup_interval);

  m_backup_vars = output_variables(m_config->get_string("output.backup_size"));

  m_last_backup_time = 0.0;
}

  //! Write a backup (i.e. an intermediate result of a run).
void IceModel::write_backup() {
  double wall_clock_hours = pism::wall_clock_hours(m_grid->com, m_start_time);

  if (wall_clock_hours - m_last_backup_time < m_backup_interval) {
    return;
  }

  const Profiling &profiling = m_ctx->profiling();
  profiling.begin("backup");

  m_last_backup_time = wall_clock_hours;

  // create a history string:
  std::string date_str = pism_timestamp(m_grid->com);
  char tmp[TEMPORARY_STRING_LENGTH];
  snprintf(tmp, TEMPORARY_STRING_LENGTH,
           "PISM automatic backup at %s, %3.3f hours after the beginning of the run\n",
           m_time->date().c_str(), wall_clock_hours);

  PetscLogDouble backup_start_time = GetTime();

  m_log->message(2,
                 "  [%s] Saving an automatic backup to '%s' (%1.3f hours after the beginning of the run)\n",
                 pism_timestamp(m_grid->com).c_str(), m_backup_filename.c_str(), wall_clock_hours);

  prepend_history(tmp);

  PIO file(m_grid->com, m_config->get_string("output.format"),
         m_backup_filename, PISM_READWRITE_MOVE);

  // write metadata:
  io::define_time(file, m_config->get_string("time.dimension_name"),
              m_time->calendar(),
              m_time->CF_units_string(),
              m_sys);
  io::append_time(file, m_config->get_string("time.dimension_name"), m_time->current());

  // Write metadata *before* variables:
  write_mapping(file);
  write_run_stats(file);
  write_global_attributes(file);
  write_config(file);

  define_model_state(file);
  define_diagnostics(file, m_backup_vars, PISM_DOUBLE);

  write_model_state(file);
  write_diagnostics(file, m_backup_vars);

  // Also flush time-series:
  flush_timeseries();

  PetscLogDouble backup_end_time = GetTime();
  m_log->message(2,
                 "  [%s] Done saving an automatic backup in %f seconds (%f minutes).\n",
                 pism_timestamp(m_grid->com).c_str(),
                 backup_end_time - backup_start_time,
                 (backup_end_time - backup_start_time) / 60.0);

  profiling.end("backup");
}


} // end of namespace pism
