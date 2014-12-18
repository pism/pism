/* Copyright (C) 2013, 2014 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "PISMOceanKill.hh"
#include "pism_options.hh"
#include "PISMVars.hh"
#include "Mask.hh"
#include <assert.h>
#include "error_handling.hh"

namespace pism {

OceanKill::OceanKill(IceGrid &g)
  : Component(g) {

  m_ocean_kill_mask.create(m_grid, "ocean_kill_mask", WITH_GHOSTS,
                           m_config.get("grid_max_stencil_width"));

  m_ocean_kill_mask.set_attrs("internal",
                              "mask specifying fixed calving front locations",
                              "", "");
  m_ocean_kill_mask.set_time_independent(true);
}

OceanKill::~OceanKill() {
  // empty
}

void OceanKill::init() {
  verbPrintf(2, m_grid.com,
             "* Initializing calving at a fixed calving front...\n");

  options::String ocean_kill_file("-ocean_kill_file",
                                  "Specifies a file to get ocean_kill thickness from");

  IceModelVec2S thickness, *tmp;

  if (ocean_kill_file.is_set()) {
    verbPrintf(2, m_grid.com,
               "  setting fixed calving front location using\n"
               "  ice thickness from '%s'.\n", ocean_kill_file.c_str());

    thickness.create(m_grid, "thk", WITHOUT_GHOSTS);
    thickness.set_attrs("temporary", "land ice thickness",
                        "m", "land_ice_thickness");
    thickness.metadata().set_double("valid_min", 0.0);

    thickness.regrid(ocean_kill_file, CRITICAL);

    tmp = &thickness;
  } else {
    throw RuntimeError("option -ocean_kill_file is required.");
  }

  IceModelVec2Int *mask = m_grid.variables().get_2d_mask("mask");
  MaskQuery m(*mask);

  IceModelVec::AccessList list;
  list.add(m_ocean_kill_mask);
  list.add(*tmp);
  list.add(*mask);

  for (Points p(m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((*tmp)(i, j) > 0 || m.grounded(i, j)) { // FIXME: use GeometryCalculator
      m_ocean_kill_mask(i, j) = 0;
    } else {
      m_ocean_kill_mask(i, j) = 1;
    }
  }

  m_ocean_kill_mask.update_ghosts();
}

// Updates mask and ice thickness, including ghosts.
void OceanKill::update(IceModelVec2Int &pism_mask, IceModelVec2S &ice_thickness) {
  IceModelVec::AccessList list;
  list.add(m_ocean_kill_mask);
  list.add(pism_mask);
  list.add(ice_thickness);

  unsigned int GHOSTS = pism_mask.get_stencil_width();
  assert(m_ocean_kill_mask.get_stencil_width() >= GHOSTS);
  assert(ice_thickness.get_stencil_width()     >= GHOSTS);

  for (PointsWithGhosts p(m_grid, GHOSTS); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_ocean_kill_mask(i, j) > 0.5) {
      pism_mask(i, j)     = MASK_ICE_FREE_OCEAN;
      ice_thickness(i, j) = 0.0;
    }
  }
}

void OceanKill::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  if (keyword == "medium" || keyword == "big") {
    result.insert(m_ocean_kill_mask.metadata().get_string("short_name"));
  }
}

void OceanKill::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                               IO_Type nctype) {

  if (set_contains(vars, m_ocean_kill_mask.metadata().get_string("short_name"))) {
    m_ocean_kill_mask.define(nc, nctype);
  }
}

void OceanKill::write_variables(const std::set<std::string> &vars, const PIO& nc) {

  if (set_contains(vars, m_ocean_kill_mask.metadata().get_string("short_name"))) {
    m_ocean_kill_mask.write(nc);
  }
}

} // end of namespace pism
