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

OceanKill::OceanKill(IceGrid &g, const Config &conf)
  : Component(g, conf) {
  PetscErrorCode ierr;

  ierr = m_ocean_kill_mask.create(grid, "ocean_kill_mask", WITH_GHOSTS,
                                  config.get("grid_max_stencil_width"));
  if (ierr != 0) {
    throw std::runtime_error("OceanKill allocation failed");
  }

  ierr = m_ocean_kill_mask.set_attrs("internal",
                                     "mask specifying fixed calving front locations",
                                     "", ""); CHKERRCONTINUE(ierr);
  m_ocean_kill_mask.set_time_independent(true);
}

OceanKill::~OceanKill() {
  // empty
}

PetscErrorCode OceanKill::init(Vars &vars) {
  PetscErrorCode ierr;
  std::string filename;
  bool flag;

  ierr = verbPrintf(2, grid.com,
                    "* Initializing calving at a fixed calving front...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Fixed calving front options", "");
  PISM_PETSC_CHK(ierr, "PetscOptionsBegin");
  {
    ierr = OptionsString("-ocean_kill_file", "Specifies a file to get ocean_kill thickness from",
                             filename, flag); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd();
  PISM_PETSC_CHK(ierr, "PetscOptionsEnd");

  IceModelVec2Int *mask = vars.get_2d_mask("mask");
  MaskQuery m(*mask);

  IceModelVec2S thickness, *tmp;

  if (flag == false) {
    throw RuntimeError("option -ocean_kill_file is required.");
  } else {
    ierr = verbPrintf(2, grid.com,
       "  setting fixed calving front location using\n"
       "  ice thickness from '%s'.\n", filename.c_str()); CHKERRQ(ierr);

    ierr = thickness.create(grid, "thk", WITHOUT_GHOSTS); CHKERRQ(ierr);
    ierr = thickness.set_attrs("temporary", "land ice thickness",
                               "m", "land_ice_thickness"); CHKERRQ(ierr);
    thickness.metadata().set_double("valid_min", 0.0);

    ierr = thickness.regrid(filename, CRITICAL); CHKERRQ(ierr);

    tmp = &thickness;
  }

  IceModelVec::AccessList list;
  list.add(m_ocean_kill_mask);
  list.add(*tmp);
  list.add(*mask);

  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if ((*tmp)(i, j) > 0 || m.grounded(i, j)) // FIXME: use GeometryCalculator
      m_ocean_kill_mask(i, j) = 0;
    else
      m_ocean_kill_mask(i, j) = 1;
  }

  ierr = m_ocean_kill_mask.update_ghosts(); CHKERRQ(ierr);

  return 0;
}

// Updates mask and ice thickness, including ghosts.
PetscErrorCode OceanKill::update(IceModelVec2Int &pism_mask, IceModelVec2S &ice_thickness) {
  IceModelVec::AccessList list;
  list.add(m_ocean_kill_mask);
  list.add(pism_mask);
  list.add(ice_thickness);

  unsigned int GHOSTS = pism_mask.get_stencil_width();
  assert(m_ocean_kill_mask.get_stencil_width() >= GHOSTS);
  assert(ice_thickness.get_stencil_width()     >= GHOSTS);

  for (PointsWithGhosts p(grid, GHOSTS); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_ocean_kill_mask(i, j) > 0.5) {
      pism_mask(i, j)     = MASK_ICE_FREE_OCEAN;
      ice_thickness(i, j) = 0.0;
    }
  }

  return 0;
}

void OceanKill::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  if (keyword == "medium" || keyword == "big")
    result.insert(m_ocean_kill_mask.metadata().get_string("short_name"));
}

PetscErrorCode OceanKill::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                               IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, m_ocean_kill_mask.metadata().get_string("short_name"))) {
    ierr = m_ocean_kill_mask.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode OceanKill::write_variables(const std::set<std::string> &vars, const PIO& nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, m_ocean_kill_mask.metadata().get_string("short_name"))) {
    ierr = m_ocean_kill_mask.write(nc); CHKERRQ(ierr);
  }

  return 0;
}

} // end of namespace pism
