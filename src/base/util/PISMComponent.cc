// Copyright (C) 2008-2016 Ed Bueler and Constantine Khroulev
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

#include <gsl/gsl_math.h>
#include <cassert>

#include "PISMComponent.hh"
#include "base/util/io/PIO.hh"
#include "IceGrid.hh"
#include "pism_const.hh"
#include "pism_utilities.hh"
#include "VariableMetadata.hh"
#include "iceModelVec.hh"
#include "pism_options.hh"
#include "error_handling.hh"
#include "PISMConfigInterface.hh"
#include "MaxTimestep.hh"

namespace pism {

/*! Process command-line options -i and -bootstrap.
 *
 */
InputOptions process_input_options(MPI_Comm com) {
  InputOptions result;

  options::String input_filename("-i", "Specifies the PISM input file");
  bool bootstrap_is_set = options::Bool("-bootstrap", "enable bootstrapping heuristics");

  const bool bootstrap = input_filename.is_set() and bootstrap_is_set;
  const bool restart   = input_filename.is_set() and not bootstrap_is_set;

  result.filename = input_filename;

  if (restart) {
    // re-start a run by initializing from an input file
    result.type = INIT_RESTART;
  } else if (bootstrap) {
    // initialize from an input file using bootstrapping heuristics
    result.type = INIT_BOOTSTRAP;
  } else {
    // other types of initialization (usually from formulas)
    result.type = INIT_OTHER;
  }

  // get the index of the last record in the input file
  if (input_filename.is_set()) {
    PIO input_file(com, "guess_mode", input_filename, PISM_READONLY);

    // Find the index of the last record in the input file.
    unsigned int last_record = input_file.inq_nrecords();
    if (last_record > 0) {
      last_record -= 1;
    }

    result.record = last_record;
  } else {
    result.record = 0;
  }

  return result;
}

Component::Component(IceGrid::ConstPtr g)
  : m_grid(g), m_config(g->ctx()->config()), m_sys(g->ctx()->unit_system()),
    m_log(g->ctx()->log()) {
  // empty
}

Component::~Component() {
  // empty
}

std::map<std::string, Diagnostic::Ptr> Component::diagnostics() const {
  return this->diagnostics_impl();
}

std::map<std::string, TSDiagnostic::Ptr> Component::ts_diagnostics() const {
  return this->ts_diagnostics_impl();
}

std::map<std::string, Diagnostic::Ptr> Component::diagnostics_impl() const {
  return {};
}

std::map<std::string, TSDiagnostic::Ptr> Component::ts_diagnostics_impl() const {
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

  options::String regrid_file("-regrid_file", "regridding file name");

  options::StringSet regrid_vars("-regrid_vars",
                                 "comma-separated list of regridding variables",
                                 "");

  if (not regrid_file.is_set()) {
    return;
  }

  SpatialVariableMetadata &m = variable.metadata();

  if ((regrid_vars.is_set() and set_contains(regrid_vars, m.get_string("short_name"))) or
      (not regrid_vars.is_set() and flag == REGRID_WITHOUT_REGRID_VARS)) {

    m_log->message(2,
               "  %s: regridding '%s' from file '%s' ...\n",
               module_name.c_str(),
               m.get_string("short_name").c_str(), regrid_file->c_str());

    variable.regrid(regrid_file, CRITICAL);
  }
}

Component_TS::Component_TS(IceGrid::ConstPtr g)
  : Component(g) {
  m_t = m_dt = GSL_NAN;
}

Component_TS::~Component_TS() {
  // empty
}

MaxTimestep Component_TS::max_timestep(double t) const {
  return this->max_timestep_impl(t);
}

void Component_TS::update(double t, double dt) {
  this->update_impl(t, dt);
}

/*!
 * Update a `model` by asking it to perform time-stepping from the current time to one year in the
 * future (or as far as the time step restriction allows).
 *
 * This is sometimes necessary during initialization, but should be avoided if possible.
 */
void init_step(Component_TS &model, const Time& time) {
  const double
    now               = time.current(),
    one_year_from_now = time.increment_date(now, 1.0);

  // Take a one year long step if we can.
  MaxTimestep max_dt(one_year_from_now - now);

  max_dt = std::min(max_dt, model.max_timestep(now));

  // Do not take time-steps shorter than 1 second
  if (max_dt.value() < 1.0) {
    max_dt = MaxTimestep(1.0);
  }

  assert(max_dt.finite() == true);

  model.update(now, max_dt.value());
}


} // end of namespace pism
