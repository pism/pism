// Copyright (C) 2004-2024 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <algorithm>
#include <set>

#include "pism/icemodel/IceModel.hh"

#include "pism/util/Grid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Diagnostic.hh"
#include "pism/util/Time.hh"
#include "pism/util/io/File.hh"

#include "pism/util/io/io_helpers.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/projection.hh"
#include "pism/util/Component.hh"


namespace pism {

MaxTimestep reporting_max_timestep(const std::vector<double> &times, double t,
                                   double resolution,
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
    j = gsl_interp_bsearch(times.data(), t, 0, N - 1);
  }

  dt = times[j + 1] - t;

  // now make sure that we don't end up taking a time-step of less than "resolution"
  // second long
  if (dt < resolution) {
    if (j + 2 < N) {
      return MaxTimestep(times[j + 2] - t, description);
    }
    return MaxTimestep(description);
  }
  return MaxTimestep(dt, description);
}

//! Write time-independent metadata to a file.
void IceModel::write_metadata(const File &file, MappingTreatment mapping_flag,
                              HistoryTreatment history_flag) const {
  if (mapping_flag == WRITE_MAPPING) {
    write_mapping(file, m_grid->get_mapping_info());
  }

  m_config->write(file);

  if (history_flag == PREPEND_HISTORY) {
    VariableMetadata tmp = m_output_global_attributes;

    std::string old_history = file.read_text_attribute("PISM_GLOBAL", "history");

    tmp.set_name("PISM_GLOBAL");
    tmp["history"] = std::string(tmp["history"]) + old_history;

    io::write_attributes(file, tmp, io::PISM_DOUBLE);
  } else {
    io::write_attributes(file, m_output_global_attributes, io::PISM_DOUBLE);
  }
}

//! Save model state in NetCDF format.
/*!
Calls save_variables() to do the actual work.
 */
void IceModel::save_results() {
  {
    auto stats = run_stats();

    auto str = pism::printf(
      "PISM done. Performance stats: %.4f wall clock hours, %.4f proc.-hours, %.4f model years per proc.-hour.",
      (double)stats["wall_clock_hours"],
      (double)stats["processor_hours"],
      (double)stats["model_years_per_processor_hour"]);

    prepend_history(str);
  }

  std::string filename = m_config->get_string("output.file");

  if (filename.empty()) {
    m_log->message(2, "WARNING: output.file is empty. Using unnamed.nc instead.\n");
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
    File file(m_grid->com,
              filename,
              string_to_backend(m_config->get_string("output.format")),
              io::PISM_READWRITE_MOVE,
              m_ctx->pio_iosys_id());

    write_metadata(file, WRITE_MAPPING, PREPEND_HISTORY);

    write_run_stats(file, run_stats());

    save_variables(file, INCLUDE_MODEL_STATE, m_output_vars,
                   m_time->current());
  }
  profiling.end("io.model_state");
}

void write_mapping(const File &file, const pism::MappingInfo &info) {

  const auto &mapping = info.mapping;
  std::string name = mapping.get_name();
  if (mapping.has_attributes()) {
    if (not file.variable_exists(name)) {
      file.define_variable(name, io::PISM_DOUBLE, {});
    }
    io::write_attributes(file, mapping, io::PISM_DOUBLE);

    // Write the PROJ string to mapping:proj_params (for CDO).
    std::string proj = info.proj;
    if (not proj.empty()) {
      file.write_attribute(name, "proj_params", proj);
    }
  }
}

void write_run_stats(const File &file, const pism::VariableMetadata &stats) {
  if (not file.variable_exists(stats.get_name())) {
    file.define_variable(stats.get_name(), io::PISM_DOUBLE, {});
  }
  io::write_attributes(file, stats, io::PISM_DOUBLE);
}

void IceModel::save_variables(const File &file,
                              OutputKind kind,
                              const std::set<std::string> &variables,
                              double time,
                              io::Type default_diagnostics_type) const {

  // Compress 2D and 3D variables if output.compression_level > 0 and the output.format
  // supports it.
  file.set_compression_level(m_config->get_number("output.compression_level"));

  // define the time dimension if necessary (no-op if it is already defined)
  io::define_time(file, *m_grid->ctx());

  // define the "timestamp" (wall clock time since the beginning of the run)
  // Note: it is time-dependent, so we need to define time first.
  VariableMetadata timestamp("timestamp", m_sys);
  timestamp.long_name("wall-clock time since the beginning of the run").units("hours");
  io::define_timeseries(timestamp,
                        m_config->get_string("time.dimension_name"),
                        file, io::PISM_FLOAT);
  // append to the time dimension
  io::append_time(file, *m_config, time);

  // Write metadata *before* everything else:
  //
  // FIXME: we should write this to variables instead of attributes because NetCDF-4 crashes after
  // about 2^16 attribute modifications per variable. :-(
  write_run_stats(file, run_stats());

  if (kind == INCLUDE_MODEL_STATE) {
    define_model_state(file);
  }
  define_diagnostics(file, variables, default_diagnostics_type);

  // Done defining variables

  {
    // Note: we don't use "variables" (an argument of this method) here because it
    // contains PISM's names of diagnostic quantities which (in some cases) map to more
    // than one NetCDF variable. Moreover, here we're concerned with file contents, not
    // the list of requested variables.
    std::set<std::string> var_names;
    unsigned int n_vars = file.nvariables();
    for (unsigned int k = 0; k < n_vars; ++k) {
      var_names.insert(file.variable_name(k));
    }

    // If this output file contains variables lat and lon...
    if (member("lat", var_names) and member("lon", var_names)) {

      // add the coordinates attribute to all variables that use x and y dimensions
      for (const auto& v : var_names) {
        std::set<std::string> dims;
        for (const auto& d : file.dimensions(v)) {
          dims.insert(d);
        }

        if (not member(v, {"lat", "lon", "lat_bnds", "lon_bnds"}) and
            member("x", dims) and member("y", dims)) {
          file.write_attribute(v, "coordinates", "lat lon");
        }
      }

      // and if it also contains lat_bnds and lon_bnds, add the bounds attribute to lat
      // and lon.
      if (member("lat_bnds", var_names) and member("lon_bnds", var_names)) {
        file.write_attribute("lat", "bounds", "lat_bnds");
        file.write_attribute("lon", "bounds", "lon_bnds");
      }
    }
  }

  if (kind == INCLUDE_MODEL_STATE) {
    write_model_state(file);
  }
  write_diagnostics(file, variables);

  // find out how much time passed since the beginning of the run and save it to the output file
  {
    unsigned int time_length = file.dimension_length(m_config->get_string("time.dimension_name"));
    size_t start = time_length > 0 ? static_cast<size_t>(time_length - 1) : 0;
    io::write_timeseries(file, timestamp, start,
                         {wall_clock_hours(m_grid->com, m_start_time)});
  }
}

void IceModel::define_diagnostics(const File &file, const std::set<std::string> &variables,
                                  io::Type default_type) const {
  for (const auto& variable : variables) {
    auto diag = m_diagnostics.find(variable);

    if (diag != m_diagnostics.end()) {
      diag->second->define(file, default_type);
    }
  }
}

//! \brief Writes variables listed in vars to filename, using nctype to write
//! fields stored in dedicated Arrays.
void IceModel::write_diagnostics(const File &file, const std::set<std::string> &variables) const {
  for (const auto& variable : variables) {
    auto diag = m_diagnostics.find(variable);

    if (diag != m_diagnostics.end()) {
      diag->second->compute()->write(file);
    }
  }
}

void IceModel::define_model_state(const File &file) const {
  for (auto *v : m_model_state) {
    v->define(file, io::PISM_DOUBLE);
  }

  for (const auto& m : m_submodels) {
    m.second->define_model_state(file);
  }

  for (const auto& d : m_diagnostics) {
    d.second->define_state(file);
  }
}

void IceModel::write_model_state(const File &file) const {
  for (auto *v : m_model_state) {
    v->write(file);
  }

  for (const auto& m : m_submodels) {
    m.second->write_model_state(file);
  }

  for (const auto& d : m_diagnostics) {
    d.second->write_state(file);
  }
}

} // end of namespace pism
