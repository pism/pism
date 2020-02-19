/* Copyright (C) 2013, 2014, 2015, 2016 PISM Authors
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
#include <cassert>

#include "pism/geometry/Geometry.hh"
#include "pism/util/Vars.hh"
#include "pism/util/Time.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_options.hh"

#include "LakeCC.hh"
#include "LakeLevel_ConnectedComponents.hh"

namespace pism {
namespace ocean {
namespace lake_level {

LakeCC::LakeCC(IceGrid::ConstPtr g)
  : LakeLevel(g),
    m_gc(*m_config),
    m_target_level(grid, "target_lake_level", WITHOUT_GHOSTS) {

  m_log->message(2, "  - Setting up LakeCC Module...\n");

  m_option = "-lakecc";

  //Patch
  m_patch_iter = m_config->get_number("lake_level.lakecc.max_patch_iterations");
  m_max_update_interval_years = m_config->get_number("lake_level.lakecc.max_update_interval",
                                "years");
  if (m_max_update_interval_years < 0) {
    //Ful update every timestep requested
    m_max_update_interval_years = 0;
  }

  m_next_update_time = m_grid->ctx()->time()->current();

}

LakeCC::~LakeCC() {
  // empty
}

void LakeCC::init_impl(const Geometry &geometry) {
}

void LakeCC::update_impl(const Geometry &geometry, double t, double dt) {

  const IceModelVec2S &old_sl = geometry.sea_level_elevation,
                      &new_sl = *m_grid->variables().get_2d_scalar("sea_level");

  bool full_update = false;

  //Check is a complete update is due!
  if ((t >= m_next_update_time) or (fabs(t - m_next_update_time) < 1.0)) {
    full_update = true;
  }

  if (!full_update) {
    //Full update when ocean basins have vanished.
    full_update = checkOceanBasinsVanished(geometry.bed_elevation,
                                           old_sl, new_sl);
  }

  if (!full_update) {
    //Full update when patch iteration does not finish.
    full_update = iterativelyPatchTargetLevel(geometry.bed_elevation,
                                              geometry.ice_thickness,
                                              new_sl, m_target_level);
  }

  if (full_update) {
    //Update Target lake level using LakeCC model!
    updateTargetLevel();

    //Do update after the 'fake' initialization step
    if (t != m_grid->ctx()->time()->start()) {
      m_next_update_time = m_grid->ctx()->time()->increment_date(t, m_max_update_interval_years);
    }
  }

}

MaxTimestep LakeCC::max_timestep_impl(double t) const {

  if (m_max_update_interval_years > 0.0) {
    double dt = m_next_update_time - t;
    if (dt < 1.0) {
      double update_time_after_next = m_grid->ctx()->time()->increment_date(m_next_update_time,
                                                                            m_max_update_interval_years);
      dt = update_time_after_next - m_next_update_time;
      assert(dt > 0.0);
    }

    MaxTimestep lakecc_dt(dt, "lake level forcing");
    return lakecc_dt;
  } else {
    MaxTimestep lakecc_dt("lake level forcing");
    return lakecc_dt;
  }

}

bool LakeCC::expandMargins_impl() const {
  return true;
}













bool LakeCC::checkOceanBasinsVanished(const IceModelVec2S &bed,
                                      const IceModelVec2S &old_sl,
                                      const IceModelVec2S &new_sl) {
  IceModelVec2Int mask(m_grid, "mask", WITHOUT_GHOSTS);
  IceModelVec2S min_basin(m_grid, "min_basin", WITHOUT_GHOSTS),
                max_wl(m_grid, "max_wl", WITHOUT_GHOSTS);

  { //Check which ocean cells are newly added
    ParallelSection ParSec(m_grid->com);
    try {
      FilterExpansionCC FExCC(m_grid, m_fill_value, bed, old_sl);
      FExCC.filter_ext2(old_sl, new_sl, mask, min_basin, max_wl);
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }

  bool basinVanished = false;

  IceModelVec::AccessList list{ &mask };
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    if (mask(i, j) == -1) {
      basinVanished = true;
      break;
    }
  }

  return (GlobalOr(m_grid->com, basinVanished));
}

bool LakeCC::iterativelyPatchTargetLevel(const IceModelVec2S &bed,
                                         const IceModelVec2S &thk,
                                         const IceModelVec2S &sl,
                                         IceModelVec2S &target_level) {

  for (unsigned int n = 0; n < m_patch_iter; ++n) {
    const int unsigned local_patch_result = patch_lake_levels(bed, thk, sl, target_level),
                       patch_result = GlobalMax(m_grid->com, local_patch_result);
    if (patch_result == 0) {
      //No further iteration needed
      return false;
    }
    if ((patch_result == 2) or (n >= (m_patch_iter -1))) {
      //update needed
      return true;
    }
  }
}

unsigned int LakeCC::patch_lake_levels(const IceModelVec2S &bed,
                                       const IceModelVec2S &thk,
                                       const IceModelVec2S &sea_level,
                                       IceModelVec2S &lake_level) {
  const Direction dirs[] = { North, East, South, West };

  IceModelVec2S lake_level_old(m_grid, "ll", WITH_GHOSTS, 1);
  lake_level_old.copy_from(lake_level);

  IceModelVec::AccessList list{ &lake_level, &thk, &bed, &sea_level, &lake_level_old };

  unsigned int return_value = 0;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double bed_ij = bed(i, j),
                 thk_ij = thk(i, j),
                 sl_ij  = sea_level(i, j);
    if (not m_gc.islake(lake_level_old(i, j))) {
      StarStencil<double> ll_star = lake_level_old.star(i, j);

      double Level        = m_fill_value;
      unsigned int nLevel = 0;
      bool becomesLake    = false;

      for (int n = 0; n < 4; ++n) {
        const Direction direction = dirs[n];
        if (m_gc.islake(ll_star[direction])) {
          if (Level != ll_star[direction]) {
            nLevel++;
            Level = ll_star[direction];
            if (mask::ocean(m_gc.mask(m_fill_value, bed_ij, thk_ij, Level))) {
              // Cell below Level -> mark to be filled
              becomesLake = true;
            }
            if (becomesLake and (nLevel > 1)) {
              // Lake cell adjoins several lakes -> full Update
              return 2;
            }
          }
        }
      }
      if ((nLevel > 0) and mask::ocean(m_gc.mask(sl_ij, bed_ij, thk_ij))) {
        // Cell is now ocean and is next to lakes -> recalculate
        return 2;
      }
      if (becomesLake) {
        // Set cell to become a lake
        lake_level(i, j) = Level;
        return_value = 1;
      }
    } else { // cell was lake
      if (mask::ocean(m_gc.mask(sl_ij, bed_ij, thk_ij))) {
        // Was lake and is now Ocean -> full Update
        return 2;
      } else {
        if (not mask::ocean(m_gc.mask(m_fill_value, bed_ij, thk_ij, lake_level_old(i, j)))) {
          // If cell that was previously lake is not lake anymore -> remove
          // label
          lake_level(i, j) = m_fill_value;
          return_value = 1;
        }
      }
    }
  }
  // If function gets here all cells have been checked and patched
  return return_value;
}

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism

