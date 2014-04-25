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

namespace pism {

OceanKill::OceanKill(IceGrid &g, const Config &conf)
  : Component(g, conf) {
  PetscErrorCode ierr;

  ierr = m_ocean_kill_mask.create(grid, "ocean_kill_mask", WITH_GHOSTS,
                                  grid.max_stencil_width);
  if (ierr != 0) {
    PetscPrintf(grid.com, "PISM ERROR: failed to allocate storage for OceanKill.\n");
    PISMEnd();
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

  ierr = PetscOptionsBegin(grid.com, "", "Fixed calving front options", ""); CHKERRQ(ierr);
  {
    ierr = OptionsString("-ocean_kill_file", "Specifies a file to get ocean_kill thickness from",
                             filename, flag); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  IceModelVec2Int *mask = dynamic_cast<IceModelVec2Int*>(vars.get("mask"));
  assert(mask != NULL);
  MaskQuery m(*mask);

  IceModelVec2S thickness, *tmp;

  if (flag == false) {
    PetscPrintf(grid.com, "PISM ERROR: option -ocean_kill_file is required.\n");
    PISMEnd();
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

  ierr = m_ocean_kill_mask.begin_access(); CHKERRQ(ierr);
  ierr = tmp->begin_access();              CHKERRQ(ierr);
  ierr = mask->begin_access();             CHKERRQ(ierr);

  for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if ((*tmp)(i, j) > 0 || m.grounded(i, j) ) // FIXME: use GeometryCalculator
        m_ocean_kill_mask(i, j) = 0;
      else
        m_ocean_kill_mask(i, j) = 1;
    }
  }

  ierr = mask->end_access();             CHKERRQ(ierr);
  ierr = tmp->end_access();              CHKERRQ(ierr);
  ierr = m_ocean_kill_mask.end_access(); CHKERRQ(ierr);

  ierr = m_ocean_kill_mask.update_ghosts(); CHKERRQ(ierr);

  return 0;
}

// Updates mask and ice thickness, including ghosts.
PetscErrorCode OceanKill::update(IceModelVec2Int &pism_mask, IceModelVec2S &ice_thickness) {
  PetscErrorCode ierr;

  ierr = m_ocean_kill_mask.begin_access(); CHKERRQ(ierr);
  ierr = pism_mask.begin_access();         CHKERRQ(ierr);
  ierr = ice_thickness.begin_access();     CHKERRQ(ierr);
  int GHOSTS = grid.max_stencil_width;
  for (int   i = grid.xs - GHOSTS; i < grid.xs+grid.xm + GHOSTS; ++i) {
    for (int j = grid.ys - GHOSTS; j < grid.ys+grid.ym + GHOSTS; ++j) {
      if (m_ocean_kill_mask(i, j) > 0.5) {
        pism_mask(i, j)     = MASK_ICE_FREE_OCEAN;
        ice_thickness(i, j) = 0.0;
      }
    }
  }
  ierr = ice_thickness.end_access();     CHKERRQ(ierr);
  ierr = pism_mask.end_access();         CHKERRQ(ierr);
  ierr = m_ocean_kill_mask.end_access(); CHKERRQ(ierr);

  return 0;
}

void OceanKill::add_vars_to_output(std::string keyword, std::set<std::string> &result) {
  if (keyword == "medium" || keyword == "big")
    result.insert(m_ocean_kill_mask.metadata().get_string("short_name"));
}

PetscErrorCode OceanKill::define_variables(std::set<std::string> vars, const PIO &nc,
                                               IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, m_ocean_kill_mask.metadata().get_string("short_name"))) {
    ierr = m_ocean_kill_mask.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode OceanKill::write_variables(std::set<std::string> vars, const PIO& nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, m_ocean_kill_mask.metadata().get_string("short_name"))) {
    ierr = m_ocean_kill_mask.write(nc); CHKERRQ(ierr);
  }

  return 0;
}

} // end of namespace pism
