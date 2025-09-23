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
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "pism/icemodel/IceModel.hh"

#include "pism/util/Grid.hh"
#include "pism/util/Config.hh"
#include "pism/util/Diagnostic.hh"
#include "pism/util/Time.hh"
#include "pism/util/io/File.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/projection.hh"
#include "pism/util/Component.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/io/io_helpers.hh"

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

std::set<VariableMetadata> IceModel::metadata(MappingTreatment mapping_flag,
                                              RunStatsTreatment run_stats) const {
  std::set<VariableMetadata> result{};

  // mapping
  if (mapping_flag == WRITE_MAPPING) {
    auto info = m_grid->get_mapping_info();

    auto mapping = info.cf_mapping;
    mapping.set_output_type(io::PISM_INT);

    if (not info.proj_string.empty()) {
      // Write the PROJ string to mapping:proj_params (for CDO).
      mapping["proj_params"] = info.proj_string;
    }

    if (mapping.has_attributes()) {
      result.insert(mapping);
    }
  }

  // config
  result.insert(config_metadata(*m_config));

  // global attributes
  {
    auto global       = m_output_global_attributes;
    global["history"] = m_output_history;
    result.insert(global);
  }

  // run stats
  if (run_stats == WRITE_RUN_STATS) {
    VariableMetadata wall_clock("wall_clock_time", m_sys);
    wall_clock.long_name("wall-clock time since the beginning of the run")
        .units("hours")
        .set_time_dependent(true)
        .set_output_type(io::PISM_FLOAT);

    VariableMetadata myph("model_years_per_processor_hour", m_sys);
    myph.long_name(
            "average number of model years per processor hour, since the beginning of the run")
        .units("years / hour")
        .set_time_dependent(true)
        .set_output_type(io::PISM_FLOAT);

    VariableMetadata step_counter("step_counter", m_sys);
    step_counter.long_name("number of time steps since the beginning of the run")
        .units("")
        .set_time_dependent(true)
        .set_output_type(io::PISM_INT);

    result.insert({ wall_clock, myph, step_counter });
  }

  return result;
}

//! Define metadata variables and attributes in an output file
void IceModel::define_metadata(const OutputFile &file, MappingTreatment mapping_flag,
                               RunStatsTreatment run_stats) const {

  for (const auto &v : metadata(mapping_flag, run_stats)) {
    if (v.get_name() == "PISM_GLOBAL") {
      file.append_history(v["history"]);

      // clear the history attribute
      auto tmp = v;
      tmp["history"] = "";
      file.set_global_attributes(tmp.all_strings(), tmp.all_doubles());
    } else {
      file.define_variable(v);
    }
  }
}

void IceModel::write_metadata(const OutputFile &file) const {
  write_config(*m_config, "pism_config", file);
}

//! Save model state in NetCDF format.
/*!
Calls write_variables() to do the actual work.
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

    // define the time dimension if necessary (no-op if it is already defined)
    {
      bool with_bounds = false;
      io::define_time(file, m_time->metadata(), with_bounds);
      define_metadata(file, WRITE_MAPPING, WRITE_RUN_STATS);
    }

    define_variables(file, INCLUDE_MODEL_STATE, m_output_vars);

    write_metadata(file);
    write_variables(file, INCLUDE_MODEL_STATE, m_output_vars, m_time->current());
  }
  profiling.end("io.model_state");
}

void IceModel::write_run_stats(const OutputFile &file) const {

  auto t_length = file.time_dimension_length();
  auto t_start = t_length > 0 ? t_length - 1 : 0;

  double wall_clock_hours = pism::wall_clock_hours(m_grid->com, m_start_time),
         proc_hours       = m_grid->size() * wall_clock_hours,
         model_years = m_time->convert_time_interval(m_time->current() - m_time->start(), "years");

  std::vector<unsigned int> start = { t_start };
  std::vector<unsigned int> count = { 1 };

  if (not m_config->get_string("output.experiment_id").empty()) {
    start.insert(start.cbegin(), 0);
    count.insert(count.cbegin(), 1);
  }

  file.write_array({ "wall_clock_time", m_sys }, start, count, { wall_clock_hours });
  file.write_array({ "model_years_per_processor_hour", m_sys }, start, count,
                   { model_years / proc_hours });
  file.write_array({ "step_counter", m_sys }, start, count, { (double)m_step_counter });
}

void IceModel::define_variables(const OutputFile &file, OutputKind kind,
                                const std::set<std::string> &variables) const {

  auto time_name = m_time->variable_name();

  if (kind == INCLUDE_MODEL_STATE) {
    define_state(file);
  }

  define_diagnostics(file, variables);
}

void IceModel::write_variables(const OutputFile &file, OutputKind kind,
                               const std::set<std::string> &variables, double time_seconds) const {
  // append to the time dimension
  file.append_time(time_seconds);

  if (kind == INCLUDE_MODEL_STATE) {
    write_state(file);
  }
  write_diagnostics(file, variables);

  write_run_stats(file);
}

void define_variables(std::set<SpatialVariableMetadata> &variables, const OutputFile &file,
                      bool use_internal_units, const std::string &mapping_variable_name) {
  std::set<std::string> variable_names;
  for (const auto &v : variables) {
    variable_names.insert(v.get_name());
  }

  for (auto var : variables) {

    if (use_internal_units) {
      var.output_units(var["units"]);
    }

    auto var_name = var.get_name();

    if (var_name == "lat" and set_member("lat_bnds", variable_names)) {
      var["bounds"] = "lat_bnds";
    }
    if (var_name == "lon" and set_member("lon_bnds", variable_names)) {
      var["bounds"] = "lon_bnds";
    }

    // add extra attributes such as "grid_mapping" and "coordinates". Variables lat, lon,
    // lat_bnds, and lon_bnds should not have these attributes to support CDO (see issue
    // #384).
    //
    // We check names of x and y dimensions to avoid setting extra attributes for variables
    // that use a different grid (e.g. viscous_bed_displacement written by the Lingle-Clark
    // bed deformation model).
    auto dim_names                = var.dimension_names();
    std::set<std::string> lat_lon = { "lat_bnds", "lon_bnds", "lat", "lon" };
    if (vector_member("x", dim_names) and vector_member("y", dim_names)) {
      if (not set_member(var_name, lat_lon) and set_member("lat", variable_names) and
          set_member("lon", variable_names)) {
        var["coordinates"] = "lat lon";
      }

      if (not mapping_variable_name.empty()) {
        var["grid_mapping"] = mapping_variable_name;
      }
    }

    file.define_variable(var);
  }
}

void IceModel::define_diagnostics(const OutputFile &file,
                                  const std::set<std::string> &variable_names) const {
  std::set<SpatialVariableMetadata> variables;

  std::string mapping_variable_name{};
  {
    const auto &mapping = m_grid->get_mapping_info();
    if (not mapping.proj_string.empty() or mapping.cf_mapping.has_attributes()) {
      mapping_variable_name = mapping.cf_mapping.get_name();
    }
  }

  for (const auto &var : variable_names) {
    auto diag = m_diagnostics.find(var);

    if (diag != m_diagnostics.end()) {
      const auto &D = diag->second;
      for (unsigned int k = 0; k < D->n_variables(); ++k) {
        variables.insert(D->metadata(k));
      }
    }
  }

  pism::define_variables(variables, file, m_config->get_flag("output.use_MKS"),
                         mapping_variable_name);
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

void IceModel::define_state(const OutputFile &file) const {
  std::string mapping_variable_name{};
  {
    const auto &mapping = m_grid->get_mapping_info();
    if (not mapping.proj_string.empty() or mapping.cf_mapping.has_attributes()) {
      mapping_variable_name = mapping.cf_mapping.get_name();
    }
  }

  std::set<SpatialVariableMetadata> state_variables{};
  for (auto *v : m_model_state) {
    for (unsigned int k = 0; k < v->ndof(); ++k) {
      state_variables.insert(v->metadata(k));
    }
  }

  pism::define_variables(state_variables, file, m_config->get_flag("output.use_MKS"),
                         mapping_variable_name);

  for (const auto& m : m_submodels) {
    m.second->define_state(file);
  }

  for (const auto& d : m_diagnostics) {
    d.second->define_state(file);
  }
}

void IceModel::write_state(const OutputFile &file) const {
  for (auto *v : m_model_state) {
    v->write(file);
  }

  for (const auto& m : m_submodels) {
    m.second->write_state(file);
  }

  for (const auto& d : m_diagnostics) {
    d.second->write_state(file);
  }
}

std::string IceModel::save_state_on_error(const std::string &suffix,
                                          const std::set<std::string> &additional_variables) {
  std::string filename = m_config->get_string("output.file");

  if (filename.empty()) {
    m_log->message(2, "WARNING: output.file is empty. Using unnamed.nc instead.");
    filename = "unnamed.nc";
  }

  filename = filename_add_suffix(filename, suffix, "");

  auto variables = output_variables("small");
  for (const auto &v : additional_variables) {
    variables.insert(v);
  }

  OutputFile file(m_output_writer, filename);

  // define the time dimension if necessary (no-op if it is already defined)
  {
    bool with_bounds = false;
    io::define_time(file, m_time->metadata(), with_bounds);
    define_metadata(file, WRITE_MAPPING, WRITE_RUN_STATS);
  }

  define_variables(file, INCLUDE_MODEL_STATE, variables);

  write_metadata(file);
  write_variables(file, INCLUDE_MODEL_STATE, variables, m_time->current());

  file.close();

  return filename;
}

} // end of namespace pism
