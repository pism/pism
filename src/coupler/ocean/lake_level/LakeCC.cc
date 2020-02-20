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
    m_target_level(m_grid, "target_lake_level", WITHOUT_GHOSTS) {

  m_log->message(2, "  - Setting up LakeCC Model...\n");

  m_option = "-lakecc";

  //Patch
  m_patch_iter = m_config->get_number("lake_level.lakecc.max_patch_iterations");
  m_max_update_interval_years = m_config->get_number("lake_level.lakecc.max_update_interval",
                                "years");
  if (m_max_update_interval_years < 0) {
    //Ful update every timestep requested
    m_max_update_interval_years = 0;
  }

  //LakeCC
  m_icefree_thickness = m_config->get_number("geometry.ice_free_thickness_standard");

  const double ice_density        = m_config->get_number("constants.ice.density"),
               freshwater_density = m_config->get_number("constants.fresh_water.density");
  m_drho = ice_density / freshwater_density;

  m_lake_level_min = m_config->get_number("lake_level.lakecc.zmin");
  m_lake_level_max = m_config->get_number("lake_level.lakecc.zmax");
  m_lake_level_dh  = m_config->get_number("lake_level.lakecc.dz");

  m_filter_size = m_config->get_number("lake_level.lakecc.filter_size");
  if (m_filter_size < 0) {
    //No filtering at all
    m_filter_size = 0;
  }

  m_check_sl_diagonal = m_config->get_flag("lake_level.lakecc.check_sl_diagonal");
  m_keep_existing_lakes = m_config->get_flag("lake_level.lakecc.keep_existing_lakes");

  m_topg_overlay.create(m_grid, "topg_overlay", WITHOUT_GHOSTS);
  m_topg_overlay.set_attrs("internal",
                           "topography overlay",
                           "meter", "meter", "", 0);

}

LakeCC::~LakeCC() {
  // empty
}

void LakeCC::init_impl(const Geometry &geometry) {

  m_log->message(2, "  *Initializing LakeCC model.\n");

  m_log->message(2, "  LakeCC: lake levels between %g and %gm, with %gm spacing\n",
                 m_lake_level_min, m_lake_level_max, m_lake_level_dh);

  if (m_filter_size > 0) {
    m_log->message(3, "  *LakeCC: Filter size: %i \n", m_filter_size);
  }

 {
    std::string overlay_file = m_config->get_string("lake_level.lakecc.topg_overlay_file");

    if (not overlay_file.empty()) {
      m_topg_overlay.regrid(overlay_file, OPTIONAL, 0.0);
      m_use_topg_overlay = true;
    } else {
      m_topg_overlay.set(0.0);
      m_use_topg_overlay = false;
    }
  }

  InputOptions opts = process_input_options(m_grid->com, m_config);

  if (opts.type == INIT_RESTART) {

    m_log->message(2, "* Reading lake level forcing from '%s' for re-starting...\n",
                    opts.filename.c_str());

    File file(m_grid->com, opts.filename, PISM_GUESS, PISM_READONLY);
    const unsigned int time_length = file.nrecords(),
                      last_record = time_length > 0 ? time_length - 1 : 0;

    m_lake_level.read(file, last_record);

    file.close();
  } else if (opts.type == INIT_BOOTSTRAP) {
    try {
      //effective_lake_level might be available in input file
      m_lake_level.regrid(opts.filename, CRITICAL);
    } catch (...) {
      //if it was not found...
      m_lake_level.set(m_fill_value);

      bool init_filled = m_config->get_flag("lake_level.lakecc.init_filled");
      if (init_filled) {
        m_log->message(2, "  *LakeCC: init with filled basins requested. Running LakeCC model.\n");

        const IceModelVec2S &sea_level = *m_grid->variables().get_2d_scalar("sea_level");
        updateLakeCC(geometry.bed_elevation,
                     geometry.ice_thickness,
                     sea_level,
                     m_lake_level,
                     m_lake_level);
      }
    }
  }

  m_log->message(3, "  *LakeCC: number of iterations used by patch-algorithm: %d \n", m_patch_iter);
  m_log->message(3, "  *LakeCC: maximum interval of full update: %d years \n", m_max_update_interval_years);

  //Full update in first timestep
  m_next_update_time = m_grid->ctx()->time()->current();
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
    updateLakeCC(geometry.bed_elevation,
                 geometry.ice_thickness,
                 new_sl,
                 geometry.lake_level_elevation,
                 m_target_level);

    if (t != m_grid->ctx()->time()->start()) {
      //Reset next update time.
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

void LakeCC::updateLakeCC(const IceModelVec2S& bed,
                          const IceModelVec2S& thk,
                          const IceModelVec2S& sea_level,
                          const IceModelVec2S& eff_lake_level,
                          IceModelVec2S& lake_level) {

  //Add topography overlay onto PISM bed if available
  IceModelVec2S lakecc_bed(m_grid, "filtered_bed_elevation", WITHOUT_GHOSTS);
  if (m_use_topg_overlay) {
    m_topg_overlay.add(1.0, bed, lakecc_bed);
  } else {
    lakecc_bed.copy_from(bed);
  }



  //do lake update
  {
    //Create PISM's mask
    IceModelVec2Int pism_mask(m_grid, "pism_mask", WITHOUT_GHOSTS);
    m_gc.compute_mask(sea_level, lakecc_bed, thk, pism_mask);



    if (m_check_sl_diagonal) {
      IceModelVec2Int mask_wide(m_grid, "mask_wide", WITH_GHOSTS);
      mask_wide.copy_from(pism_mask);

      IceModelVec2S sl_wide(m_grid, "sl_wide", WITH_GHOSTS);
      sl_wide.copy_from(sea_level);

      std::vector<int>i_diagonals = {1, 1,-1,-1},
                      j_diagonals = {1,-1,-1, 1};

      IceModelVec::AccessList list({ &mask_wide, &pism_mask, &sl_wide, &lakecc_bed, &thk });

      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        if (sl_wide(i, j) == m_fill_value) {
          for (int k = 0; k < 4; k++) {
            int i_diag = i + i_diagonals[k],
                j_diag = j + j_diagonals[k];
            if (mask::ocean(mask_wide(i_diag, j_diag))) {
              pism_mask(i, j) = m_gc.mask(sl_wide(i_diag, j_diag), lakecc_bed(i, j), thk(i, j));
              break;
            }
          }
        }
      }
      pism_mask.update_ghosts();
    }



    //Prepare validity mask. It marks areas invalid that are completely surrounded by ice.
    IceModelVec2Int validity_mask(m_grid, "pism_mask", WITHOUT_GHOSTS);
    {
      IsolationCC IsoCC(m_grid, thk, m_icefree_thickness);
      IsoCC.find_isolated_spots(validity_mask);

      if (m_keep_existing_lakes) {
        IceModelVec::AccessList list{ &validity_mask, &bed, &thk, &eff_lake_level };

        if (eff_lake_level.state_counter() > 0) {
          for (Points p(*m_grid); p; p.next()) {
            const int i = p.i(), j = p.j();

            //Set valid where a lake already exists
            if ( mask::ocean(m_gc.mask(m_fill_value, bed(i, j), thk(i, j), eff_lake_level(i, j))) ) {
              validity_mask(i, j) = 1;
            }
          }
        }
      }
      validity_mask.update_ghosts();
    }



    m_log->message(2, "->LakeCC: Update of Lake Levels! \n");

    ParallelSection ParSec(m_grid->com);
    try {
      // Initialze LakeCC Model
      LakeLevelCC LM(m_grid, m_drho, bed, thk, pism_mask, m_fill_value, validity_mask);
      LM.computeLakeLevel(m_lake_level_min, m_lake_level_max, m_lake_level_dh, lake_level);
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();

    m_log->message(2, "          Done!\n");
  }



  //Transfer back onto PISM topography
  if (m_use_topg_overlay) {
    IceModelVec::AccessList list({ &lake_level, &bed, &thk });

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      const double ll_ij = lake_level(i, j);
      // If cell is labled as lake on higher "resolved" bed,
      // but is not valid on topography used by PISM -> set invalid
      if ( m_gc.islake(ll_ij) and not
          mask::ocean(m_gc.mask(m_fill_value, bed(i, j), thk(i, j), ll_ij)) ) {
        lake_level(i, j) = m_fill_value;
      }
    }
  }



  //Filter lakes
  if (m_filter_size > 0) {
    ParallelSection ParSec(m_grid->com);
    try {
      FilterLakesCC FL(m_grid, m_fill_value);
      FL.filter_map(m_filter_size, lake_level);
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }

  lake_level.update_ghosts();
}

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism

