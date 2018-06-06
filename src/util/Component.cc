// Copyright (C) 2008-2018 Ed Bueler and Constantine Khroulev
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

#include <cassert>

#include "Component.hh"
#include "pism/util/io/PIO.hh"
#include "IceGrid.hh"
#include "pism_utilities.hh"
#include "VariableMetadata.hh"
#include "iceModelVec.hh"
#include "pism_options.hh"
#include "error_handling.hh"
#include "ConfigInterface.hh"
#include "MaxTimestep.hh"
#include "pism/util/Time.hh"

namespace pism {

InputOptions::InputOptions(InitializationType t, const std::string &file, unsigned int index) {
  type     = t;
  filename = file;
  record   = index;
}

/*! Process command-line options -i and -bootstrap.
 *
 */
InputOptions process_input_options(MPI_Comm com, Config::ConstPtr config) {
  InitializationType type   = INIT_OTHER;
  unsigned int       record = 0;

  std::string input_filename = config->get_string("input.file");

  bool bootstrap = config->get_boolean("input.bootstrap") and (not input_filename.empty());
  bool restart   = (not config->get_boolean("input.bootstrap")) and (not input_filename.empty());

  if (restart) {
    // re-start a run by initializing from an input file
    type = INIT_RESTART;
  } else if (bootstrap) {
    // initialize from an input file using bootstrapping heuristics
    type = INIT_BOOTSTRAP;
  } else {
    // other types of initialization (usually from formulas)
    type = INIT_OTHER;
  }

  // get the index of the last record in the input file
  if (not input_filename.empty()) {
    PIO input_file(com, "guess_mode", input_filename, PISM_READONLY);

    // Find the index of the last record in the input file.
    unsigned int last_record = input_file.inq_nrecords();
    if (last_record > 0) {
      last_record -= 1;
    }

    record = last_record;
  } else {
    record = 0;
  }

  return InputOptions(type, input_filename, record);
}

Component::Component(IceGrid::ConstPtr g)
  : m_grid(g), m_config(g->ctx()->config()), m_sys(g->ctx()->unit_system()),
    m_log(g->ctx()->log()) {
  // empty
}

Component::~Component() {
  // empty
}

DiagnosticList Component::diagnostics() const {
  return this->diagnostics_impl();
}

TSDiagnosticList Component::ts_diagnostics() const {
  return this->ts_diagnostics_impl();
}

DiagnosticList Component::diagnostics_impl() const {
  return {};
}

TSDiagnosticList Component::ts_diagnostics_impl() const {
  return {};
}

IceGrid::ConstPtr Component::grid() const {
  return m_grid;
}

/*! @brief Define model state variables in an output file. */
/*!
 * This is needed to allow defining all the variables in an output file before any data is written
 * (an optimization needed to get decent performance writing NetCDF-3).
 */
void Component::define_model_state(const PIO &output) const {
  this->define_model_state_impl(output);
}

/*! @brief Write model state variables to an output file. */
void Component::write_model_state(const PIO &output) const {
  // define variables, if needed (this is a no-op if they are already defined)
  this->define_model_state(output);

  this->write_model_state_impl(output);
}

/*! @brief The default (empty implementation). */
void Component::define_model_state_impl(const PIO &output) const {
  (void) output;
}

/*! @brief The default (empty implementation). */
void Component::write_model_state_impl(const PIO &output) const {
  (void) output;
}

/**
 * Regrid a variable by processing -regrid_file and -regrid_vars.
 *
 * @param[in] module_name Module name, used to annotate options when run with -help.
 *
 * @param[out] variable pointer to an IceModelVec; @c variable has to
 *             have metadata set for this to work.
 *
 * @param[in] flag Regridding flag. If set to
 *            REGRID_WITHOUT_REGRID_VARS, regrid this variable by
 *            default, if =-regrid_vars= was not set. Otherwise a
 *            variable is only regridded if both =-regrid_file= and
 *            =-regrid_vars= are set *and* the name of the variable is
 *            found in the set of names given with =-regrid_vars=.
 */
void Component::regrid(const std::string &module_name, IceModelVec &variable,
                       RegriddingFlag flag) {

  auto regrid_file = m_config->get_string("input.regrid.file");
  auto regrid_vars = set_split(m_config->get_string("input.regrid.vars"), ',');

  if (regrid_file.empty()) {
    return;
  }

  SpatialVariableMetadata &m = variable.metadata();

  if (((not regrid_vars.empty()) and member(m.get_string("short_name"), regrid_vars)) or
      (regrid_vars.empty() and flag == REGRID_WITHOUT_REGRID_VARS)) {

    m_log->message(2,
               "  %s: regridding '%s' from file '%s' ...\n",
               module_name.c_str(),
               m.get_string("short_name").c_str(), regrid_file.c_str());

    variable.regrid(regrid_file, CRITICAL);
  }
}

MaxTimestep Component::max_timestep(double t) const {
  return this->max_timestep_impl(t);
}

MaxTimestep Component::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep();
}


} // end of namespace pism
