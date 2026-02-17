// Copyright (C) 2004-2026 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include "pism/util/Component.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/io/SynchronousOutputWriter.hh"

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

void IceModel::define_time(const OutputFile &file, bool with_bounds) const {
  file.define_variable(m_time->metadata(with_bounds));

  if (with_bounds) {
    file.define_variable(m_time->bounds_metadata());
  }
}

std::set<VariableMetadata> IceModel::common_metadata() const {
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

  std::set<VariableMetadata> result = { wall_clock, myph, step_counter, m_output_global_attributes,
                                        config_metadata(*m_config) };

  {
    auto mapping = m_grid->get_mapping_info();

    if (mapping.has_attributes()) {
      result.insert(mapping);
    }
  }

  return result;
}

void IceModel::init_final_output() {
  m_output_filename = m_config->get_string("output.file");

  if (m_output_filename.empty()) {
    m_output_filename = "unnamed.nc";
    m_log->message(2, "WARNING: output.file is empty. Using '%s' instead.\n",
                   m_output_filename.c_str());
  }

  if (not ends_with(m_output_filename, ".nc")) {
    m_log->message(2, "PISM WARNING: output file name does not have the '.nc' suffix!\n");
  }

  m_output_vars = output_variables(m_config->get_string("output.size"));

#if (Pism_USE_PROJ == 1)
  {
    std::string proj_string = m_grid->get_mapping_info()["proj_params"];
    if (not proj_string.empty()) {
      for (std::string v : {"lon", "lat"}) {
        if (set_member(v, m_output_vars)) {
          m_output_vars.insert(v + "_bnds");
        }
      }
    }
  }
#endif

  m_output_file_contents = pism::combine(common_metadata(), state_variables());
  m_output_file_contents =
      pism::combine(m_output_file_contents, state_variables_diagnostics(m_output_vars));
  m_output_file_contents =
      pism::combine(m_output_file_contents, state_variables_diagnostics(m_extra_vars));
  m_output_file_contents =
      pism::combine(m_output_file_contents, diagnostic_variables(m_output_vars));
}

//! Save model state in NetCDF format.
/*!
Calls write_variables() to do the actual work.
 */
void IceModel::write_final_output() {
  append_history("PISM done");

  const Profiling &profiling = m_ctx->profiling();

  profiling.begin("io.model_state");
  if (m_config->get_string("output.size") != "none") {
    m_log->message(2, "Writing model state to file `%s'...\n", m_output_filename.c_str());
    OutputFile file(m_output_writer, m_output_filename);

    {
      define_time(file);
      define_variables(file, m_output_file_contents);
    }

    {
      io::write_config(*m_config, "pism_config", file);
      file.append_time(m_time->current());
      write_state(file);
      write_state_diagnostics(file, m_output_vars);
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

  file.write_array("wall_clock_time", start, count, { wall_clock_hours });
  file.write_array("model_years_per_processor_hour", start, count, { model_years / proc_hours });
  file.write_array("step_counter", start, count, { (double)m_step_counter });
}

void IceModel::define_variables(const OutputFile &file,
                                const std::set<VariableMetadata> &variables) const {
  io::define_variables(file, variables,
                       m_grid->get_mapping_info(),
                       m_config->get_flag("output.use_MKS"));
}

/*!
 * Return the set of state variables from IceModel and all its sub-models. Does not
 * include state variables from requested diagnostics.
 */
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
  }

  return result;
}

/*!
 * Write state variables of IceModel and all its sub-models. Does not include state
 * variables corresponding to requested diagnostic variables.
 */
void IceModel::write_state(const OutputFile &file) const {
  for (auto *v : m_model_state) {
    v->write(file);
  }

  for (const auto& m : m_submodels) {
    m.second->write_state(file);
  }
}

std::string IceModel::save_state_on_error(const std::string &suffix,
                                          const std::set<std::string> &additional_variables) {

  auto filename = filename_add_suffix(m_output_filename, suffix, "");

  std::shared_ptr<OutputWriter> writer =
      std::make_shared<SynchronousOutputWriter>(m_grid->com, *m_config);

  auto variables = pism::combine(m_output_file_contents, diagnostic_variables(additional_variables));

  std::set<std::string> variable_names;
  for (const auto &v : variables) {
    variable_names.insert(v.get_name());
  }

  writer->initialize(variables);

  OutputFile file(writer, filename);

  {
    define_time(file);
    define_variables(file, variables);
  }

  {
    io::write_config(*m_config, "pism_config", file);
    file.append_time(m_time->current());
    write_state(file);
    write_state_diagnostics(file, variable_names);
    write_diagnostics(file, variable_names);
    write_run_stats(file);
  }

  file.close();

  return filename;
}

} // end of namespace pism
