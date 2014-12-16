// Copyright (C) 2010, 2011, 2012, 2013, 2014 Constantine Khroulev
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

BedDef::BedDef(IceGrid &g)
  : Component_TS(g) {

  thk    = NULL;
  topg   = NULL;
  uplift = NULL;

  t_beddef_last = GSL_NAN;

  PetscErrorCode ierr = pismbeddef_allocate();
  if (ierr != 0) {
    throw std::runtime_error("BedDef allocation failed");
  }
}

PetscErrorCode BedDef::pismbeddef_allocate() {
  const unsigned int WIDE_STENCIL = config.get("grid_max_stencil_width");

  topg_initial.create(grid, "topg_initial", WITH_GHOSTS, WIDE_STENCIL);
  topg_initial.set_attrs("model_state", "bedrock surface elevation (at the beginning of the run)",
                         "m", "");

  topg_last.create(grid, "topg", WITH_GHOSTS, WIDE_STENCIL);
  topg_last.set_attrs("model_state", "bedrock surface elevation",
                      "m", "bedrock_altitude");

  return 0;
}

void BedDef::add_vars_to_output(const std::string &/*keyword*/, std::set<std::string> &result) {
  result.insert("topg_initial");
}

void BedDef::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                            IO_Type nctype) {
  if (set_contains(vars, "topg_initial")) {
    topg_initial.define(nc, nctype);
  }
}

void BedDef::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  if (set_contains(vars, "topg_initial")) {
    topg_initial.write(nc);
  }
}

void BedDef::init(Vars &vars) {
  t_beddef_last = grid.time->start();

  thk    = vars.get_2d_scalar("land_ice_thickness");
  topg   = vars.get_2d_scalar("bedrock_altitude");
  uplift = vars.get_2d_scalar("tendency_of_bedrock_altitude");

  // Save the bed elevation at the beginning of the run:
  topg_initial.copy_from(*topg);
}

//! Compute bed uplift (dt_beddef is in seconds).
void BedDef::compute_uplift(double dt_beddef) {
  topg->add(-1, topg_last, *uplift);
  //! uplift = (topg - topg_last) / dt
  uplift->scale(1.0 / dt_beddef);
}

} // end of namespace pism
