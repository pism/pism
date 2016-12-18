// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Constantine Khroulev
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

#include "PISMBedDef.hh"
#include "base/util/pism_utilities.hh"
#include "base/util/PISMTime.hh"
#include "base/util/PISMVars.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMConfigInterface.hh"

namespace pism {
namespace bed {

BedDef::BedDef(IceGrid::ConstPtr g)
  : Component_TS(g) {

  m_t_beddef_last = GSL_NAN;

  const unsigned int WIDE_STENCIL = m_config->get_double("grid.max_stencil_width");

  m_topg.create(m_grid, "topg", WITH_GHOSTS, WIDE_STENCIL);
  m_topg.set_attrs("model_state", "bedrock surface elevation",
                   "m", "bedrock_altitude");

  m_topg_initial.create(m_grid, "topg_initial", WITH_GHOSTS, WIDE_STENCIL);
  m_topg_initial.set_attrs("model_state",
                           "bedrock surface elevation (at the beginning of the run)",
                           "m", "");

  m_topg_last.create(m_grid, "topg", WITH_GHOSTS, WIDE_STENCIL);
  m_topg_last.set_attrs("model_state", "bedrock surface elevation",
                        "m", "bedrock_altitude");

  m_uplift.create(m_grid, "dbdt", WITHOUT_GHOSTS);
  m_uplift.set_attrs("model_state", "bedrock uplift rate",
                     "m s-1", "tendency_of_bedrock_altitude");
  m_uplift.metadata().set_string("glaciological_units", "m year-1");
  m_uplift.write_in_glaciological_units = true;

  // Set default values (we set them early so that pismv can override
  // them in IceCompModel::set_vars_from_options(), which is called
  // before BedDef::init()).
  m_topg.set(0.0);
  m_uplift.set(0.0);
}

BedDef::~BedDef() {
  // empty
}

void BedDef::set_elevation(const IceModelVec2S &input) {
  m_topg.copy_from(input);
  m_topg.update_ghosts();
}

void BedDef::set_uplift(const IceModelVec2S &input) {
  m_uplift.copy_from(input);
}

const IceModelVec2S& BedDef::bed_elevation() const {
  return m_topg;
}

const IceModelVec2S& BedDef::uplift() const {
  return m_uplift;
}


void BedDef::define_model_state_impl(const PIO &output) const {
  m_topg_initial.define(output);
}

void BedDef::write_model_state_impl(const PIO &output) const {
  m_topg_initial.write(output);
}

void BedDef::init() {
  this->init_impl();
}

//! Initialize using provided bed elevation and uplift.
void BedDef::init(const IceModelVec2S &bed,
                  const IceModelVec2S &bed_uplift,
                  const IceModelVec2S &ice_thickness) {
  this->init_with_inputs_impl(bed, bed_uplift, ice_thickness);
}

void BedDef::init_with_inputs_impl(const IceModelVec2S &bed,
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
void BedDef::init_impl() {
  m_t_beddef_last = m_grid->ctx()->time()->start();

  m_t  = GSL_NAN;
  m_dt = GSL_NAN;

  InputOptions opts = process_input_options(m_grid->com);

  switch (opts.type) {
  case INIT_RESTART:
    // read bed elevation and uplift rate from file
    m_log->message(2,
                   "    reading bed topography and bed uplift rate\n"
                   "    from %s ... \n",
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
      // do nothing, keeping values set elsewhere
    }
  }

  // process -regrid_file and -regrid_vars
  regrid("bed deformation", m_topg);
  regrid("bed deformation", m_uplift);

  std::string correction_file = m_config->get_string("bed_deformation.bed_topography_delta_file");
  if (not correction_file.empty()) {
    apply_topg_offset(correction_file);
  }

  // this should be the last thing we do here
  m_topg_initial.copy_from(m_topg);
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
