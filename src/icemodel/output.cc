// Copyright (C) 2004-2018 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "IceModel.hh"

#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Diagnostic.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_options.hh"

#include "pism/util/Vars.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/projection.hh"
#include "pism/util/Component.hh"
#include "pism/energy/utilities.hh"

namespace pism {

MaxTimestep reporting_max_timestep(const std::vector<double> &times, double t,
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

//! Write time-independent metadata to a file.
void IceModel::write_metadata(const PIO &file, MappingTreatment mapping_flag,
                              HistoryTreatment history_flag) {
  if (mapping_flag == WRITE_MAPPING) {
    write_mapping(file);
  }

  m_config->write(file);

  if (history_flag == PREPEND_HISTORY) {
    VariableMetadata tmp = m_output_global_attributes;

    std::string old_history = file.get_att_text("PISM_GLOBAL", "history");

    tmp.set_name("PISM_GLOBAL");
    tmp.set_string("history", tmp.get_string("history") + old_history);

    io::write_attributes(file, tmp, PISM_DOUBLE);
  } else {
    io::write_attributes(file, m_output_global_attributes, PISM_DOUBLE);
  }
}

//! Save model state in NetCDF format.
/*!
Calls save_variables() to do the actual work.
 */
void IceModel::save_results() {
  {
    update_run_stats();

    auto str = pism::printf(
      "PISM done. Performance stats: %.4f wall clock hours, %.4f proc.-hours, %.4f model years per proc.-hour.",
      m_run_stats.get_double("wall_clock_hours"),
      m_run_stats.get_double("processor_hours"),
      m_run_stats.get_double("model_years_per_processor_hour"));

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

  const Profiling &profiling = m_ctx->profiling();

  profiling.begin("io.model_state");
  if (m_config->get_string("output.size") != "none") {
    m_log->message(2, "Writing model state to file `%s'...\n", filename.c_str());
    PIO file(m_grid->com, m_config->get_string("output.format"), filename, PISM_READWRITE_MOVE);

    write_metadata(file, WRITE_MAPPING, PREPEND_HISTORY);

    write_run_stats(file);

    save_variables(file, INCLUDE_MODEL_STATE, m_output_vars);
  }
  profiling.end("io.model_state");
}

void IceModel::write_mapping(const PIO &file) {
  // only write mapping if it is set.
  const VariableMetadata &mapping = m_grid->get_mapping_info().mapping;
  std::string name = mapping.get_name();
  if (mapping.has_attributes()) {
    if (not file.inq_var(name)) {
      file.redef();
      file.def_var(name, PISM_DOUBLE, {});
    }
    io::write_attributes(file, mapping, PISM_DOUBLE);

    // Write the PROJ.4 string to mapping:proj4_params (for CDO).
    std::string proj4 = m_grid->get_mapping_info().proj4;
    if (not proj4.empty()) {
      file.put_att_text(name, "proj4_params", proj4);
    }
  }
}

void IceModel::write_run_stats(const PIO &file) {
  update_run_stats();
  if (not file.inq_var(m_run_stats.get_name())) {
    file.redef();
    file.def_var(m_run_stats.get_name(), PISM_DOUBLE, {});
  }
  io::write_attributes(file, m_run_stats, PISM_DOUBLE);
}

void IceModel::save_variables(const PIO &file,
                              OutputKind kind,
                              const std::set<std::string> &variables,
                              IO_Type default_diagnostics_type) {

  // define the time dimension if necessary (no-op if it is already defined)
  io::define_time(file, *m_grid->ctx());
  // define the "timestamp" (wall clock time since the beginning of the run)
  // Note: it is time-dependent, so we need to define time first.
  io::define_timeseries(m_timestamp, file, PISM_FLOAT);
  // append to the time dimension
  io::append_time(file, *m_config, m_grid->ctx()->time()->current());

  // Write metadata *before* everything else:
  //
  // FIXME: we should write this to variables instead of attributes because NetCDF-4 crashes after
  // about 2^16 attribute modifications per variable. :-(
  write_run_stats(file);

  if (kind == INCLUDE_MODEL_STATE) {
    define_model_state(file);
  }
  define_diagnostics(file, variables, default_diagnostics_type);

  if (kind == INCLUDE_MODEL_STATE) {
    write_model_state(file);
  }
  write_diagnostics(file, variables);

  // find out how much time passed since the beginning of the run and save it to the output file
  {
    unsigned int time_length = file.inq_dimlen(m_config->get_string("time.dimension_name"));
    size_t start = time_length > 0 ? static_cast<size_t>(time_length - 1) : 0;
    io::write_timeseries(file, m_timestamp, start,
                         wall_clock_hours(m_grid->com, m_start_time));
  }
}

void IceModel::define_diagnostics(const PIO &file, const std::set<std::string> &variables,
                                  IO_Type default_type) {
  for (auto variable : variables) {
    auto diag = m_diagnostics.find(variable);

    if (diag != m_diagnostics.end()) {
      diag->second->define(file, default_type);
    }
  }
}

//! \brief Writes variables listed in vars to filename, using nctype to write
//! fields stored in dedicated IceModelVecs.
void IceModel::write_diagnostics(const PIO &file, const std::set<std::string> &variables) {
  for (auto variable : variables) {
    auto diag = m_diagnostics.find(variable);

    if (diag != m_diagnostics.end()) {
      diag->second->compute()->write(file);
    }
  }
}

void IceModel::define_model_state(const PIO &file) {
  for (auto v : m_model_state) {
    v->define(file);
  }

  for (auto m : m_submodels) {
    m.second->define_model_state(file);
  }

  for (auto d : m_diagnostics) {
    d.second->define_state(file);
  }
}

void IceModel::write_model_state(const PIO &file) {
  for (auto v : m_model_state) {
    v->write(file);
  }

  for (auto m : m_submodels) {
    m.second->write_model_state(file);
  }

  for (auto d : m_diagnostics) {
    d.second->write_state(file);
  }
}

} // end of namespace pism
