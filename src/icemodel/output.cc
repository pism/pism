// Copyright (C) 2004-2025 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "pism/util/Config.hh"
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
void IceModel::write_metadata(const OutputFile &file, MappingTreatment mapping_flag) const {

  if (mapping_flag == WRITE_MAPPING) {
    auto info = m_grid->get_mapping_info();

    auto mapping = info.cf_mapping;
    if (not info.proj_string.empty()) {
      // Write the PROJ string to mapping:proj_params (for CDO).
      mapping["proj_params"] = info.proj_string;
    }

    if (mapping.has_attributes()) {
      file.define_variable(mapping, {});
    }
  }

  m_config->write(file);

  {
    file.append_history(m_output_history);
    file.set_global_attributes(m_output_global_attributes.all_strings(),
                               m_output_global_attributes.all_doubles());
  }
}

//! Save model state in NetCDF format.
/*!
Calls save_variables() to do the actual work.
 */
void IceModel::save_results() {
  append_history("PISM done");

  std::string filename = m_config->get_string("output.file");

  if (filename.empty()) {
    m_log->message(2, "WARNING: output.file is empty. Using unnamed.nc instead.\n");
    filename = "unnamed.nc";
  }

  if (not ends_with(filename, ".nc")) {
    m_log->message(2, "PISM WARNING: output file name does not have the '.nc' suffix!\n");
  }

  const Profiling &profiling = m_ctx->profiling();

  profiling.begin("io.model_state");
  if (m_config->get_string("output.size") != "none") {
    m_log->message(2, "Writing model state to file `%s'...\n", filename.c_str());
    OutputFile file(m_output_writer, filename);

    write_metadata(file, WRITE_MAPPING);

    save_variables(file, INCLUDE_MODEL_STATE, m_output_vars, m_time->current());
  }
  profiling.end("io.model_state");
}

void IceModel::define_run_stats(const OutputFile &file) const {

  auto time_name = m_time->variable_name();

  // define the "timestamp" (wall clock time since the beginning of the run)
  // Note: it is time-dependent, so we need to define time first.
  VariableMetadata wall_clock("wall_clock_time", m_sys);
  wall_clock.long_name("wall-clock time since the beginning of the run")
      .units("hours")
      .set_output_type(io::PISM_FLOAT);

  VariableMetadata myph("model_years_per_processor_hour", m_sys);
  myph.long_name("average number of model years per processor hour, since the beginning of the run")
      .units("years / hour")
      .set_output_type(io::PISM_FLOAT);

  VariableMetadata step_counter("step_counter", m_sys);
  step_counter.long_name("number of time steps since the beginning of the run")
      .units("")
      .set_output_type(io::PISM_INT);

  file.define_variable(wall_clock, { time_name });
  file.define_variable(myph, { time_name });
  file.define_variable(step_counter, { time_name });
}

void IceModel::write_run_stats(const OutputFile &file) const {

  auto time_name = m_time->variable_name();

  auto t_length = file.time_dimension_length();
  auto t_start = t_length > 0 ? t_length - 1 : 0;

  double wall_clock_hours = pism::wall_clock_hours(m_grid->com, m_start_time),
         proc_hours       = m_grid->size() * wall_clock_hours,
         model_years = m_time->convert_time_interval(m_time->current() - m_time->start(), "years");

  file.write_array({ "wall_clock_time", m_sys }, { t_start }, { 1 }, { wall_clock_hours });
  file.write_array({ "model_years_per_processor_hour", m_sys }, { t_start }, { 1 },
                   { model_years / proc_hours });
  file.write_array({ "step_counter", m_sys }, { t_start }, { 1 }, { (double)m_step_counter });
}

void IceModel::save_variables(const OutputFile &file, OutputKind kind,
                              const std::set<std::string> &variables, double time_seconds) const {

  auto time_name = m_time->variable_name();

  // define the time dimension if necessary (no-op if it is already defined)
  {
    auto time      = m_time->metadata();
    time["bounds"] = "time_bounds";

    file.define_dimension(time_name, io::PISM_UNLIMITED);
    file.define_variable(time, { time_name });
  }

  define_run_stats(file);

  // "lon" and "lat" are a part of the "model state", so we need to add "coordinates" when
  // saving the model state even if "lon" and "lat" were not requested explicitly:
  if ((member("lat", variables) and member("lon", variables)) or
      kind == INCLUDE_MODEL_STATE) {
    file.add_extra_attributes({ { "coordinates", "lat lon" } });
  }

  const auto &mapping = m_grid->get_mapping_info();
  if (not mapping.proj_string.empty() or mapping.cf_mapping.has_attributes()) {
    file.add_extra_attributes({ { "grid_mapping", mapping.cf_mapping.get_name() } });
  }

  if (kind == INCLUDE_MODEL_STATE) {
    define_model_state(file);
  }

  define_diagnostics(file, variables);

  // Done defining variables

  // append to the time dimension
  file.append_time(time_seconds);

  if (kind == INCLUDE_MODEL_STATE) {
    write_model_state(file);
  }
  write_diagnostics(file, variables);

  write_run_stats(file);
}

void IceModel::define_diagnostics(const OutputFile &file, const std::set<std::string> &variables) const {
  for (const auto& variable : variables) {
    auto diag = m_diagnostics.find(variable);

    if (diag != m_diagnostics.end()) {
      auto &D = *diag->second;
      for (unsigned int n = 0; n < D.n_variables(); ++n) {
        auto var = D.metadata(n);

        if (var.get_name() == "lat" and member("lat_bnds", variables)) {
          var["bounds"] = "lat_bnds";
        }
        if (var.get_name() == "lon" and member("lon_bnds", variables)) {
          var["bounds"] = "lon_bnds";
        }

        file.define_spatial_variable(var, m_grid->info());
      }
    }
  }
}

//! \brief Writes variables listed in vars to filename, using nctype to write
//! fields stored in dedicated Arrays.
void IceModel::write_diagnostics(const OutputFile &file, const std::set<std::string> &variables) const {
  for (const auto& variable : variables) {
    auto diag = m_diagnostics.find(variable);

    if (diag != m_diagnostics.end()) {
      diag->second->compute()->write(file);
    }
  }
}

void IceModel::define_model_state(const OutputFile &file) const {
  for (auto *v : m_model_state) {
    v->define(file);
  }

  for (const auto& m : m_submodels) {
    m.second->define_model_state(file);
  }

  for (const auto& d : m_diagnostics) {
    d.second->define_state(file);
  }
}

void IceModel::write_model_state(const OutputFile &file) const {
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
