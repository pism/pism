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

std::set<VariableMetadata> IceModel::run_stats_metadata() const {
  VariableMetadata wall_clock("wall_clock_time", m_sys);
  wall_clock.long_name("wall-clock time since the beginning of the run")
      .units("hours")
      .set_time_dependent(true)
      .set_output_type(io::PISM_FLOAT);

  VariableMetadata myph("model_years_per_processor_hour", m_sys);
  myph.long_name("average number of model years per processor hour, since the beginning of the run")
      .units("years / hour")
      .set_time_dependent(true)
      .set_output_type(io::PISM_FLOAT);

  VariableMetadata step_counter("step_counter", m_sys);
  step_counter.long_name("number of time steps since the beginning of the run")
      .units("")
      .set_time_dependent(true)
      .set_output_type(io::PISM_INT);

  return { wall_clock, myph, step_counter };
}

std::set<VariableMetadata> IceModel::common_metadata(bool with_time_bounds) const {
  std::set<VariableMetadata> result = pism::combine(
      run_stats_metadata(), { m_output_global_attributes, config_metadata(*m_config) });

  result.insert(m_time->metadata(with_time_bounds));
  if (with_time_bounds) {
    result.insert(m_time->bounds_metadata());
  }

  {
    auto mapping = m_grid->get_mapping_info();

    if (mapping.has_attributes()) {
      result.insert(mapping);
    }
  }

  return result;
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

    {
      auto variables = pism::combine(common_metadata(), state_variables());
      variables      = pism::combine(variables, diagnostic_variables(m_output_vars));

      define_variables(file, variables);
    }

    {
      write_config(*m_config, "pism_config", file);
      file.append_time(m_time->current());
      write_state(file);
      write_diagnostics(file, m_output_vars);
      write_run_stats(file);
    }
  }
  profiling.end("io.model_state");
}

void IceModel::write_run_stats(const OutputFile &file) const {

  auto t_length = file.time_dimension_length();
  auto t_start = t_length > 0 ? t_length - 1 : 0;

  double wall_clock_hours = pism::wall_clock_hours(m_grid->com, m_start_time),
         proc_hours       = m_grid->size() * wall_clock_hours,
         model_years = m_time->convert_time_interval(m_time->current() - m_time->start(), "years");

  std::vector<unsigned int> start{};
  std::vector<unsigned int> count{};
  if (m_config->get_string("output.experiment_id").empty()) {
    start = { t_start };
    count = { 1 };
  } else {
    start = { 0, t_start };
    count = { 1, 1 };
  }

  file.write_array({ "wall_clock_time", m_sys }, start, count, { wall_clock_hours });
  file.write_array({ "model_years_per_processor_hour", m_sys }, start, count,
                   { model_years / proc_hours });
  file.write_array({ "step_counter", m_sys }, start, count, { (double)m_step_counter });
}

void IceModel::define_variables(const OutputFile &file,
                                const std::set<VariableMetadata> &variables) const {
  io::define_variables(file, variables,
                       m_grid->get_mapping_info(),
                       m_config->get_flag("output.use_MKS"));
}

//! \brief Writes variables listed in vars to filename, using nctype to write
//! fields stored in dedicated Arrays.
void IceModel::write_diagnostics(const OutputFile &file,
                                 const std::set<std::string> &variable_names) const {
  for (const auto &variable : variable_names) {
    auto diag = m_diagnostics.find(variable);

    if (diag != m_diagnostics.end()) {
      diag->second->compute()->write(file);
    }
  }
}

std::set<VariableMetadata> IceModel::state_variables() const {
  std::set<VariableMetadata> result{};
  {
    // IceModel's state variables:
    for (auto *v : m_model_state) {
      for (unsigned int k = 0; k < v->ndof(); ++k) {
        result.insert(v->metadata(k));
      }
    }

    // state variables from sub-models:
    for (const auto &m : m_submodels) {
      result = pism::combine(result, m.second->state());
    }

    // state variables from diagnostics:
    for (const auto &d : m_diagnostics) {
      result = pism::combine(result, d.second->state());
    }
  }

  return result;
}

std::set<VariableMetadata>
IceModel::diagnostic_variables(const std::set<std::string> &variable_names) const {
  std::set<VariableMetadata> result{};
  {
    for (const auto &var : variable_names) {
      auto diag = m_diagnostics.find(var);

      if (diag != m_diagnostics.end()) {
        const auto &D = diag->second;
        for (unsigned int k = 0; k < D->n_variables(); ++k) {
          result.insert(D->metadata(k));
        }
      }
    }
  }
  return result;
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

  auto variable_names = output_variables("small");
  for (const auto &v : additional_variables) {
    variable_names.insert(v);
  }

  OutputFile file(m_output_writer, filename);

  {
    auto variables = pism::combine(common_metadata(), state_variables());
    variables = pism::combine(variables, diagnostic_variables(variable_names));

    define_variables(file, variables);
  }

  {
    write_config(*m_config, "pism_config", file);
    file.append_time(m_time->current());
    write_state(file);
    write_diagnostics(file, variable_names);
    write_run_stats(file);
  }

  file.close();

  return filename;
}

} // end of namespace pism
