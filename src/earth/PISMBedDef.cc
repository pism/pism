// Copyright (C) 2010, 2011, 2012, 2013, 2014, 2015 Constantine Khroulev
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

#include "PISMBedDef.hh"
#include "PISMTime.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"
#include "PISMConfig.hh"

#include <stdexcept>

namespace pism {

BedDef::BedDef(const IceGrid &g)
  : Component_TS(g) {

  m_thk    = NULL;

  m_t_beddef_last = GSL_NAN;

  PetscErrorCode ierr = pismbeddef_allocate();
  if (ierr != 0) {
    throw std::runtime_error("BedDef allocation failed");
  }
}

PetscErrorCode BedDef::pismbeddef_allocate() {
  const unsigned int WIDE_STENCIL = m_config.get("grid_max_stencil_width");

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
  m_uplift.set_glaciological_units("m year-1");
  m_uplift.write_in_glaciological_units = true;

  return 0;
}

const IceModelVec2S& BedDef::bed_elevation() const {
  return m_topg;
}

const IceModelVec2S& BedDef::uplift() const {
  return m_uplift;
}


void BedDef::add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &result) {
  result.insert("topg_initial");
}

void BedDef::define_variables(const std::set<std::string> &vars, const PIO &nc,
                              IO_Type nctype) {
  if (set_contains(vars, "topg_initial")) {
    m_topg_initial.define(nc, nctype);
  }
}

void BedDef::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  if (set_contains(vars, "topg_initial")) {
    m_topg_initial.write(nc);
  }
}

void BedDef::init() {
  m_t_beddef_last = m_grid.time->start();

  m_t  = GSL_NAN;
  m_dt = GSL_NAN;

  m_thk = m_grid.variables().get_2d_scalar("land_ice_thickness");

  std::string input_file;
  bool do_regrid = false;
  int start = -1;
  find_pism_input(input_file, do_regrid, start);

  // read bed elevation and uplift rate from file
  verbPrintf(2, m_grid.com,
             "    reading bed topography and bed uplift rate\n"
             "    from %s ... \n",
             input_file.c_str());
  if (do_regrid) {
    // bootstrapping
    m_topg.regrid(input_file, OPTIONAL,
                  config.get("bootstrapping_bed_value_no_var"));
    m_uplift.regrid(filename, OPTIONAL,
                    config.get("bootstrapping_uplift_value_no_var"));
  } else {
    // re-starting
    m_topg.read(input_file, start); // fails if not found!
    m_uplift.read(input_file, start); // fails if not found!
  }

  // process -regrid_file and -regrid_vars
  regrid("BedDef", &m_topg);
  regrid("BedDef", &m_uplift);

  // this should be the last thing we do here
  m_topg.copy_to(m_topg_initial);
  m_topg.copy_to(m_topg_last);
}

//! Compute bed uplift (dt_beddef is in seconds).
void BedDef::compute_uplift(double dt_beddef) {
  m_topg.add(-1, m_topg_last, m_uplift);
  //! uplift = (topg - topg_last) / dt
  m_uplift.scale(1.0 / dt_beddef);
}

} // end of namespace pism
