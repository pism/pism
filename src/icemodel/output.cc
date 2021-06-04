// Copyright (C) 2004-2021 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "IceModel.hh"

#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Diagnostic.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/pism_options.hh"

#include "pism/util/Vars.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/projection.hh"
#include "pism/util/Component.hh"
#include "pism/energy/utilities.hh"

#if (Pism_USE_CDIPIO==1)
#include <mpi.h>
extern "C"{
#include "cdipio.h"
#include "cdi.h"
#include "yaxt.h"
}
#endif

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
void IceModel::write_metadata(const File &file, MappingTreatment mapping_flag,
                              HistoryTreatment history_flag) {
  if (mapping_flag == WRITE_MAPPING) {
    write_mapping(file);
  }

  m_config->write(file);

  if (history_flag == PREPEND_HISTORY) {
    VariableMetadata tmp = m_output_global_attributes;

    std::string old_history = file.read_text_attribute("PISM_GLOBAL", "history");

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
      m_run_stats.get_number("wall_clock_hours"),
      m_run_stats.get_number("processor_hours"),
      m_run_stats.get_number("model_years_per_processor_hour"));

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
  int fileID = -1;
  IO_Mode mode = PISM_READWRITE_MOVE;
  if (streamIDs.count(filename) > 0) {
    fileID = streamIDs[filename];
    mode = PISM_READWRITE;
  }

  profiling.begin("io.model_state");
  if (m_config->get_string("output.size") != "none") {
    m_log->message(2, "Writing model state to file `%s'...\n", filename.c_str());
    profiling.begin("io.open");
    File file(m_grid->com,
              filename,
              string_to_backend(m_config->get_string("output.format")),
              mode,
              m_ctx->pio_iosys_id(),
              fileID,
              DimOutMap);
    file.set_split(false);
    profiling.end("io.open");

    write_metadata(file, WRITE_MAPPING, PREPEND_HISTORY);
    write_run_stats(file);

    file.file_calendar(m_time->year_length(), m_time->calendar());
    save_variables(file, INCLUDE_MODEL_STATE, m_output_vars,
                   m_time->current());
    m_sthwritten = true;
    expose_windows();

    if (file.backend() == PISM_CDI) {
      streamIDs[filename] = file.get_streamID();
    }
  }
  profiling.end("io.model_state");
}

void IceModel::write_mapping(const File &file) {
  // only write mapping if it is set.
  const VariableMetadata &mapping = m_grid->get_mapping_info().mapping;
  std::string name = mapping.get_name();
  if (mapping.has_attributes()) {
    if (not file.find_variable(name)) {
      file.define_variable(name, PISM_DOUBLE, {});
    }
    io::write_attributes(file, mapping, PISM_DOUBLE);

    // Write the PROJ string to mapping:proj_params (for CDO).
    std::string proj = m_grid->get_mapping_info().proj;
    if (not proj.empty()) {
      file.write_attribute(name, "proj_params", proj);
    }
  }
}

void IceModel::write_run_stats(const File &file) {
  update_run_stats();
  if (not file.find_variable(m_run_stats.get_name())) {
    file.define_variable(m_run_stats.get_name(), PISM_DOUBLE, {});
  }
  io::write_attributes(file, m_run_stats, PISM_DOUBLE);
}

void IceModel::save_variables(const File &file,
                              OutputKind kind,
                              const std::set<std::string> &variables,
                              double time,
                              IO_Type default_diagnostics_type,
                              bool realsave) {

  // Compress 2D and 3D variables if output.compression_level > 0 and the output.format
  // supports it.
  file.set_compression_level(m_config->get_number("output.compression_level"));

  // define the time dimension if necessary (no-op if it is already defined)
  io::define_time(file, *m_grid->ctx());
  // define the "timestamp" (wall clock time since the beginning of the run)
  // Note: it is time-dependent, so we need to define time first.
  io::define_timeseries(m_timestamp,
                        m_config->get_string("time.dimension_name"),
                        file, PISM_FLOAT);

  // Write metadata *before* everything else:
  //
  // FIXME: we should write this to variables instead of attributes because NetCDF-4 crashes after
  // about 2^16 attribute modifications per variable. :-(
  write_run_stats(file);

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
      std::string vname = file.variable_name(k);
      var_names.insert(vname);
    }

    // If this output file contains variables lat and lon...
    if (member("lat", var_names) and member("lon", var_names)) {

      // add the coordinates attribute to all variables that use x and y dimensions
      for (auto v : var_names) {
        std::set<std::string> dims;
        for (auto d : file.dimensions(v)) {
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
  file.define_vlist(); // vlist object is now immutable - call inquire_vlist to change it
  if (not realsave) return;

  io::append_time(file, *m_config, time);
  file.set_dimatt();
  file.send_diagnostics(variables);
  if (kind == INCLUDE_MODEL_STATE) {
    write_model_state(file);
  }
  file.set_beforediag(false);
  write_diagnostics(file, variables);
  // find out how much time passed since the beginning of the run and save it to the output file
  {
    unsigned int time_length = file.dimension_length(m_config->get_string("time.dimension_name"));
    size_t start = time_length > 0 ? static_cast<size_t>(time_length - 1) : 0;
    io::write_timeseries(file, m_timestamp, start,
                         {wall_clock_hours(m_grid->com, m_start_time)});
  }
}

void IceModel::open_files() {
#if (Pism_USE_CDIPIO==1)
if (string_to_backend(m_config->get_string("output.format")) == PISM_CDI) {
  if (not m_opened) {
    m_opened = true;
    int filetype = m_config->get_number("output.cdi_pio.filetype");
    // Open snap file/s
    if (m_save_snapshots) {
      std::string filename;
      int nsnap;
      if (not m_split_snapshots) {
        nsnap = 1;
      } else {
        nsnap = m_snapshot_times.size();
      }

      for (int sn = 0; sn < nsnap; sn++) {
        if (not m_split_snapshots) {
          filename = m_snapshots_filename;
        } else {
	  filename = pism::printf("%s_%s.nc",
                                  m_snapshots_filename.c_str(),
                                  m_time->date(m_snapshot_times[sn]).c_str());
        }

        int fileID = -1;
        IO_Mode mode = PISM_READWRITE_MOVE;
        {
          File file(m_grid->com,
                    filename,
                    string_to_backend(m_config->get_string("output.format")),
                    mode,
                    m_ctx->pio_iosys_id(),
                    fileID,
                    DimSnapMap,
                    filetype);
          streamIDs[filename] = file.get_streamID();

          write_metadata(file, WRITE_MAPPING, PREPEND_HISTORY);
          m_snapshots_file_is_ready = true;
          write_run_stats(file);
          file.file_calendar(-1.0, m_time->calendar());
          save_variables(file, INCLUDE_MODEL_STATE, m_snapshot_vars, m_time->current(), PISM_FLOAT, false);
          vlistIDs[filename] = file.get_vlistID();
          DimSnapMap = file.get_dimensions_map();
        }
        m_snapshots_file_is_ready = false;
      }
      m_snapshots_file_is_ready = true;
    }

    // Open Extra file/s
    if (m_save_extra) {
      std::string filename;
      int nsnap;
      if (not m_split_extra) {
        nsnap = 1;
      } else {
        if (not m_config->get_flag("output.ISMIP6")) {
          nsnap = m_extra_times.size() - 1; // First one is never saved
        } else {
          nsnap = m_extra_times.size();
        }
      }

      for (int sn = 0; sn < nsnap; sn++) {
        if (not m_split_extra) {
          filename = m_extra_filename;
        } else {
          double tt;
          if (not m_config->get_flag("output.ISMIP6")) {
            tt = m_extra_times[sn+1];
          } else {
            tt = m_extra_times[sn];
          }
          filename = pism::printf("%s_%s.nc",
                                  m_extra_filename.c_str(),
                                  m_time->date(tt).c_str());
        }
        int fileID = -1;
        IO_Mode mode = PISM_READWRITE_MOVE;
        {
          File file(m_grid->com,
                    filename,
                    string_to_backend(m_config->get_string("output.format")),
                    mode,
                    m_ctx->pio_iosys_id(),
                    fileID,
                    DimExtraMap,
                    filetype);
          streamIDs[filename] = file.get_streamID();

          std::string time_name = m_config->get_string("time.dimension_name");
          io::define_time(file, *m_ctx);
          file.write_attribute(time_name, "bounds", "time_bounds");
          io::define_time_bounds(m_extra_bounds, file);
          write_metadata(file, WRITE_MAPPING, PREPEND_HISTORY);
          m_extra_file_is_ready = true;
          write_run_stats(file);
          file.file_calendar(-1.0, m_time->calendar());
          save_variables(file,
                         m_extra_vars.empty() ? INCLUDE_MODEL_STATE : JUST_DIAGNOSTICS,
                         m_extra_vars,
                         0,
                         PISM_FLOAT,
                         false);

          if (file.backend() == PISM_CDI) {
            vlistIDs[filename] = file.get_vlistID();
            DimExtraMap = file.get_dimensions_map();
          }
        }
        m_extra_file_is_ready = false;
      }
      m_extra_file_is_ready = true;
    }

    // Open Output file
    if (m_config->get_string("output.size") != "none")
    {
      std::string filename = m_config->get_string("output.file_name");
      if (filename.empty()) {
        m_log->message(2, "WARNING: output.file_name is empty. Using unnamed.nc instead.\n");
        filename = "unnamed.nc";
      }
      if (not ends_with(filename, ".nc")) {
        m_log->message(2,
                       "PISM WARNING: output file name does not have the '.nc' suffix!\n");
      }
      int fileID = -1;
      IO_Mode mode = PISM_READWRITE_MOVE;
      File file(m_grid->com,
                filename,
                string_to_backend(m_config->get_string("output.format")),
                mode,
                m_ctx->pio_iosys_id(),
                fileID,
                DimOutMap,
                filetype);
      streamIDs[filename] = file.get_streamID();

      write_metadata(file, WRITE_MAPPING, PREPEND_HISTORY);

      write_run_stats(file);
      file.file_calendar(-1.0, m_time->calendar());
      save_variables(file, INCLUDE_MODEL_STATE, m_output_vars,
                     0, PISM_FLOAT, false);
      if (file.backend() == PISM_CDI) {
        vlistIDs[filename] = file.get_vlistID();
        DimOutMap = file.get_dimensions_map();
      }
    }
  }
}
#endif
}

void IceModel::expose_windows() {
#if (Pism_USE_CDIPIO==1)
if (string_to_backend(m_config->get_string("output.format")) == PISM_CDI) {
    if (m_sthwritten) {
      pioWriteTimestep();
      m_sthwritten = false;
    }
}
#endif
}

void IceModel::close_files() {
#if (Pism_USE_CDIPIO==1)
if (string_to_backend(m_config->get_string("output.format")) == PISM_CDI) {
  const Profiling &profiling = m_ctx->profiling();
  profiling.begin("io.close_streams");
  // Close all files
  for (auto const& streamID : streamIDs) {
    streamClose(streamID.second);
  }
  // Destroy all streams
  for (auto const& vlistID : vlistIDs) {
    vlistDestroy(vlistID.second);
  }
  profiling.end("io.close_streams");
}
#endif
}

void IceModel::define_diagnostics(const File &file, const std::set<std::string> &variables,
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
void IceModel::write_diagnostics(const File &file, const std::set<std::string> &variables) {
  for (auto variable : variables) {
    auto diag = m_diagnostics.find(variable);

    if (diag != m_diagnostics.end()) {
      diag->second->compute()->write(file);
    }
  }
}

void IceModel::define_model_state(const File &file) {
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

void IceModel::write_model_state(const File &file) {
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
