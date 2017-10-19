// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017 Constantine Khroulev
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

#include "BedDef.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Time.hh"
#include "pism/util/Vars.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"

namespace pism {
namespace bed {

BedDef::BedDef(IceGrid::ConstPtr g)
  : Component_TS(g) {

  m_t_beddef_last = GSL_NAN;

  const unsigned int WIDE_STENCIL = m_config->get_double("grid.max_stencil_width");

  m_topg.create(m_grid, "topg", WITH_GHOSTS, WIDE_STENCIL);
  m_topg.set_attrs("model_state", "bedrock surface elevation",
                   "m", "bedrock_altitude");

  m_topg_last.create(m_grid, "topg", WITH_GHOSTS, WIDE_STENCIL);
  m_topg_last.set_attrs("model_state", "bedrock surface elevation",
                        "m", "bedrock_altitude");

  m_uplift.create(m_grid, "dbdt", WITHOUT_GHOSTS);
  m_uplift.set_attrs("model_state", "bedrock uplift rate",
                     "m s-1", "tendency_of_bedrock_altitude");
  m_uplift.metadata().set_string("glaciological_units", "mm year-1");
}

BedDef::~BedDef() {
  // empty
}

const IceModelVec2S& BedDef::bed_elevation() const {
  return m_topg;
}

const IceModelVec2S& BedDef::uplift() const {
  return m_uplift;
}

void BedDef::define_model_state_impl(const PIO &output) const {
  m_uplift.define(output);
  m_topg.define(output);
}

void BedDef::write_model_state_impl(const PIO &output) const {
  m_uplift.write(output);
  m_topg.write(output);
}

std::map<std::string, Diagnostic::Ptr> BedDef::diagnostics_impl() const {
  std::map<std::string, Diagnostic::Ptr> result;
  result = {
    {"dbdt", Diagnostic::wrap(m_uplift)},
    {"topg", Diagnostic::wrap(m_topg)}
  };

  return result;
}

void BedDef::init(const InputOptions &opts) {
  this->init_impl(opts);
}

//! Initialize using provided bed elevation and uplift.
void BedDef::bootstrap(const IceModelVec2S &bed,
                       const IceModelVec2S &bed_uplift,
                       const IceModelVec2S &ice_thickness) {
  this->bootstrap_impl(bed, bed_uplift, ice_thickness);
}

void BedDef::bootstrap_impl(const IceModelVec2S &bed,
                            const IceModelVec2S &bed_uplift,
                            const IceModelVec2S &ice_thickness) {
  m_topg.copy_from(bed);
  m_uplift.copy_from(bed_uplift);
  // suppress a compiler warning:
  (void) ice_thickness;
}

void BedDef::update_impl(double t, double dt) {
  const IceModelVec2S *thk = m_grid->variables().get_2d_scalar("land_ice_thickness");
  this->update_with_thickness_impl(*thk, t, dt);
}

void BedDef::update(const IceModelVec2S &ice_thickness, double t, double dt) {
  this->update_with_thickness_impl(ice_thickness, t, dt);
}

//! Initialize from the context (input file and the "variables" database).
void BedDef::init_impl(const InputOptions &opts) {
  m_t_beddef_last = m_grid->ctx()->time()->start();

  m_t  = GSL_NAN;
  m_dt = GSL_NAN;

  switch (opts.type) {
  case INIT_RESTART:
    // read bed elevation and uplift rate from file
    m_log->message(2,
                   "    reading bed topography and uplift from %s ... \n",
                   opts.filename.c_str());
    // re-starting
    m_topg.read(opts.filename, opts.record); // fails if not found!
    m_uplift.read(opts.filename, opts.record); // fails if not found!
    break;
  case INIT_BOOTSTRAP:
    // bootstrapping
    m_topg.regrid(opts.filename, OPTIONAL,
                  m_config->get_double("bootstrapping.defaults.bed"));
    m_uplift.regrid(opts.filename, OPTIONAL,
                    m_config->get_double("bootstrapping.defaults.uplift"));
    break;
  case INIT_OTHER:
  default:
    {
      // do nothing
    }
  }

  // process -regrid_file and -regrid_vars
  regrid("bed deformation", m_topg);
  // uplift is not a part of the model state, but the user may want to take it from a -regrid_file
  // during bootstrapping
  regrid("bed deformation", m_uplift);

  std::string uplift_file = m_config->get_string("bed_deformation.bed_uplift_file");
  if (not uplift_file.empty()) {
    m_log->message(2,
                   "    reading bed uplift from %s ... \n",
                   uplift_file.c_str());
    m_uplift.regrid(uplift_file, CRITICAL);
  }

  std::string correction_file = m_config->get_string("bed_deformation.bed_topography_delta_file");
  if (not correction_file.empty()) {
    apply_topg_offset(correction_file);
  }

  // this should be the last thing we do here
  m_topg_last.copy_from(m_topg);
}

/*!
 * Apply a correction to the bed topography by reading topg_delta from filename.
 */
void BedDef::apply_topg_offset(const std::string &filename) {
  m_log->message(2, "  Adding a bed topography correction read in from %s...\n",
                 filename.c_str());

  IceModelVec2S topg_delta;
  topg_delta.create(m_grid, "topg_delta", WITHOUT_GHOSTS);
  topg_delta.set_attrs("internal", "bed topography correction", "meters", "");

  topg_delta.regrid(filename, CRITICAL);

  m_topg.add(1.0, topg_delta);
}

//! Compute bed uplift (dt_beddef is in seconds).
void BedDef::compute_uplift(double dt_beddef) {
  m_topg.add(-1, m_topg_last, m_uplift);
  //! uplift = (topg - topg_last) / dt
  m_uplift.scale(1.0 / dt_beddef);
}

} // end of namespace bed
} // end of namespace pism
