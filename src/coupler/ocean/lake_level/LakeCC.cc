/* Copyright (C) 2013, 2014, 2015, 2016, 2023 PISM Authors
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
#include "IsolationCC.hh"
#include "FilterLakesCC.hh"

namespace pism {
namespace ocean {
namespace lake_level {

LakeCC::LakeCC(IceGrid::ConstPtr g)
  : LakeLevel(g),
    m_target_level(m_grid, "target_lake_level", WITHOUT_GHOSTS),
    m_fill_rate(m_grid, "lake_fill_rate", WITHOUT_GHOSTS),
    m_gc(*m_config),
    m_topg_overlay(m_grid, "topg_overlay", WITHOUT_GHOSTS) {

  m_log->message(2, "  - Setting up LakeCC Model...\n");

  m_option = "lake_level.lakecc";

  m_target_level.set_attrs("model_state", "target lake level",
                           "meter", "meter", "target_level", 0);
  m_target_level.metadata().set_number("_FillValue", m_fill_value);

  m_fill_rate.set_attrs("model_state", "lake fill rate",
                        "meter second-1", "meter year-1", "lake_fill_rate", 0);
  m_fill_rate.metadata().set_number("_FillValue", m_fill_value);

  //Patch (deprecated - should not be used anymore)
  //max_update_interval should be set to 0 (default)
  m_patch_iter = m_config->get_number(m_option + ".max_patch_iterations");
  m_max_update_interval_years = m_config->get_number(m_option + ".max_update_interval",
                                "years");
  if (m_max_update_interval_years < 0) {
    //Ful update every timestep requested
    m_max_update_interval_years = 0;
  }

  //LakeCC
  m_icefree_thickness = m_config->get_number(m_option + ".ice_free_thickness");

  const double ice_density        = m_config->get_number("constants.ice.density"),
               freshwater_density = m_config->get_number("constants.fresh_water.density");
  m_drho = ice_density / freshwater_density;

  m_lake_level_min = m_config->get_number(m_option + ".zmin");
  m_lake_level_max = m_config->get_number(m_option + ".zmax");
  m_lake_level_dh  = m_config->get_number(m_option + ".dz");

  m_filter_size = m_config->get_number(m_option + ".filter_size");
  if (m_filter_size < 0) {
    //No filtering at all
    m_filter_size = 0;
  }

  //check_sl_diagonal (deprecated - only use default value *true*)
  m_check_sl_diagonal = m_config->get_flag(m_option + ".check_sl_diagonal");
  m_keep_existing_lakes = m_config->get_flag(m_option + ".keep_existing_lakes");

  m_topg_overlay.set_attrs("internal",
                           "topography overlay",
                           "meter", "meter", "", 0);
  m_topg_overlay.set_time_independent(true);

  //Gradual
  m_max_lake_fill_rate = m_config->get_number(m_option + ".max_fill_rate", "meter second-1");

  //use_const_fill_rate (deprecated - only use default option *true*)
  m_use_const_fill_rate = m_config->get_flag(m_option + ".use_constant_fill_rate");

  if (m_use_const_fill_rate) {
    m_fill_rate.set(m_max_lake_fill_rate);
  }
}

LakeCC::~LakeCC() {
  // empty
}

void LakeCC::init_impl(const Geometry &geometry) {

  m_log->message(2, "  *Initializing LakeCC model.\n");

  m_log->message(3, "  LakeCC: lake levels between %g and %gm, with %gm spacing\n",
                 m_lake_level_min, m_lake_level_max, m_lake_level_dh);

  if (m_filter_size > 0) {
    m_log->message(3, "  LakeCC: Filter size: %i \n", m_filter_size);
  }

 {
    std::string overlay_file = m_config->get_string(m_option + ".topg_overlay_file");

    if (not overlay_file.empty()) {
      m_topg_overlay.regrid(overlay_file, OPTIONAL, 0.0);
      m_use_topg_overlay = true;
    } else {
      m_topg_overlay.set(0.0);
      m_use_topg_overlay = false;
    }
  }

  {
    InputOptions opts = process_input_options(m_grid->com, m_config);
    IceModelVec2S tmp(m_grid, "effective_lake_level_elevation", WITHOUT_GHOSTS);
    tmp.set_attrs("diagnostic",
                  "lake level elevation, relative to the geoid",
                  "meter", "meter", "", 0);
    tmp.metadata().set_number("_FillValue", m_fill_value);

    if (opts.type == INIT_RESTART) {

      m_log->message(3, "* Reading lake level forcing from '%s' for re-starting...\n",
                     opts.filename.c_str());

      File file(m_grid->com, opts.filename, PISM_GUESS, PISM_READONLY);
      const unsigned int time_length = file.nrecords(),
                         last_record = time_length > 0 ? time_length - 1 : 0;

      tmp.read(file, last_record);

      file.close();
    } else if (opts.type == INIT_BOOTSTRAP) {
      try {
        //effective_lake_level might be available in input file
        tmp.regrid(opts.filename, CRITICAL);
      } catch (...) {
        //if it was not found...
        tmp.set(m_fill_value);

        bool init_filled = m_config->get_flag(m_option + ".init_filled");
        if (init_filled) {
          m_log->message(3, "  LakeCC: init with filled basins requested. Running LakeCC model.\n");

          const IceModelVec2S &sea_level = *m_grid->variables().get_2d_scalar("sea_level");
          updateLakeCC(geometry.bed_elevation,
                       geometry.ice_thickness,
                       sea_level,
                       tmp,
                       tmp);
        }
      }
    }
    m_lake_level.copy_from(tmp);
  }

  if (m_max_update_interval_years > 0) {
    m_log->message(1, "  LakeCC: number of iterations used by patch-algorithm: %d \n",
                  m_patch_iter);
    m_log->message(1, "  LakeCC: maximum interval of full update: %d years \n",
                  m_max_update_interval_years);
    m_log->message(1, "  -> Warning! This option is deprecated; only the default value should be used: 0 years \n");
  }

  if (not m_use_const_fill_rate) {
    m_log->message(1, "  LakeCC: use of constant fill rate disabled. \n");
    m_log->message(1, "  -> Warning! This option is deprecated; only the default value should be used: true \n");
  }

  double max_fill_rate_m_y = m_config->get_number(m_option + ".max_fill_rate", "meter year-1");
  m_log->message(3, "  LakeCC: maximum fill rate: %g meter/year \n",
                 max_fill_rate_m_y);

  //Full update in first timestep
  m_next_update_time = m_grid->ctx()->time()->current();
}

void LakeCC::update_impl(const Geometry &geometry, double t, double dt) {

  const IceModelVec2S &old_sl = geometry.sea_level_elevation,
                      &new_sl = *m_grid->variables().get_2d_scalar("sea_level"),
                      &bed    = geometry.bed_elevation,
                      &thk    = geometry.ice_thickness,
                      &old_ll = geometry.lake_level_elevation;

  bool full_update = false;

  if (m_max_update_interval_years <= 0) {
    m_next_update_time = t;
  }

  //Check is a complete update is due!
  //By defualt this is alway true (upadte every timestep)
  if ((t >= m_next_update_time) or (fabs(t - m_next_update_time) < 1.0)) {
    full_update = true;
  }

  //By default this is skipped
  if (!full_update) {
    //Full update when ocean basins have vanished.
    full_update = checkOceanBasinsVanished(bed,
                                           old_sl,
                                           new_sl);
  }

  //By default this is skipped
  if (!full_update) {
    //Full update when patch iteration does not finish.
    full_update = iterativelyPatchTargetLevel(bed,
                                              thk,
                                              new_sl,
                                              m_target_level);
  }


  if (full_update) {
    //Update Target lake level using LakeCC model!
    updateLakeCC(bed,
                 thk,
                 new_sl,
                 geometry.lake_level_elevation,
                 m_target_level);

    if (t != m_grid->ctx()->time()->start()) {
      //Reset next update time.
      m_next_update_time = m_grid->ctx()->time()->increment_date(t, m_max_update_interval_years);
    }
  }


  //Gradually fill
  {

    //By default this is skipped
    if (not m_use_const_fill_rate) {
      const IceModelVec2S &bmb               = *m_grid->variables().get_2d_scalar("effective_BMB"),
                          &tc_calving        = *m_grid->variables().get_2d_scalar("thickness_change_due_to_calving"),
                          &tc_frontal_melt   = *m_grid->variables().get_2d_scalar("thickness_change_due_to_frontal_melt"),
                          &tc_forced_retreat = *m_grid->variables().get_2d_scalar("thickness_change_due_to_forced_retreat");

      compute_fill_rate(dt,
                        m_target_level,
                        bmb,
                        tc_calving,
                        tc_frontal_melt,
                        tc_forced_retreat,
                        m_fill_rate);
    }

    IceModelVec2S min_level(m_grid, "min_level", WITHOUT_GHOSTS),
                  max_level(m_grid, "max_level", WITHOUT_GHOSTS),
                  min_basin(m_grid, "min_basin", WITHOUT_GHOSTS);

    updateLakeLevelMinMax(m_lake_level,
                          m_target_level,
                          min_level,
                          max_level);

    const bool LakeLevelChanged = prepareLakeLevel(m_target_level,
                                                   bed,
                                                   thk,
                                                   min_level,
                                                   old_ll,
                                                   old_sl,
                                                   min_basin,
                                                   m_lake_level);

    if (LakeLevelChanged) {
      //if a new lake basin was added we need to update the min and max lake level
      updateLakeLevelMinMax(m_lake_level,
                            m_target_level,
                            min_level,
                            max_level);
    }

    gradually_fill(dt,
                   m_max_lake_fill_rate,
                   m_target_level,
                   bed,
                   thk,
                   new_sl,
                   min_level,
                   max_level,
                   min_basin,
                   m_fill_rate,
                   m_lake_level);
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

    return {dt, "lake level forcing"};
  }
  return {"lake level forcing"};
}

bool LakeCC::expandMargins_impl() const {
  return true;
}

// Write diagnostic variables to extra files if requested
DiagnosticList LakeCC::diagnostics_impl() const {

  DiagnosticList result = {
    { "lakecc_gradual_target",       Diagnostic::wrap(m_target_level) },
    { "lakecc_gradual_fill_rate",    Diagnostic::wrap(m_fill_rate) },
  };

  return result;
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

  for (int n = 0; n < m_patch_iter; ++n) {
    const int unsigned local_patch_result = patch_lake_levels(bed, thk, sl, target_level),
                       patch_result = GlobalMax(m_grid->com, local_patch_result);
    if (patch_result == 0) {
      //No further iteration needed
      return false;
    }
    if ((patch_result == 2) or (n >= (m_patch_iter - 1))) {
      //update needed
      return true;
    }
  }

  // FIXME: we don't know if this is correct
  return false;
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



    m_log->message(3, "->LakeCC: Update of Lake Levels! \n");

    ParallelSection ParSec(m_grid->com);
    try {
      // Initialze LakeCC Model
      LakeLevelCC LM(m_grid, m_drho, bed, thk, pism_mask, m_fill_value, validity_mask);
      LM.computeLakeLevel(m_lake_level_min, m_lake_level_max, m_lake_level_dh, lake_level);
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();

    m_log->message(3, "          Done!\n");
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

void LakeCC::compute_fill_rate(const double dt,
                               const IceModelVec2S &lake_level,
                               const IceModelVec2S &bmb,
                               const IceModelVec2S &tc_calving,
                               const IceModelVec2S &tc_frontal_melt,
                               const IceModelVec2S &tc_forced_retreat,
                               IceModelVec2S &lake_fill_rate) {

  IceModelVec2S lake_area(m_grid, "lake_area", WITHOUT_GHOSTS),
                lake_mass_input_discharge(m_grid, "lake_mass_input_discharge", WITHOUT_GHOSTS),
                lake_mass_input_basal(m_grid, "lake_mass_input_basal", WITHOUT_GHOSTS);

  {
    //Initialize Lake accumulator
    LakeAccumulatorCCSerial Lacc(m_grid, m_fill_value);
    Lacc.init(lake_level);

    IceModelVec2S mass_discharge(m_grid, "mass_discharge", WITHOUT_GHOSTS),
                  mass_basal(m_grid, "mass_basal", WITHOUT_GHOSTS);

    double rho_ice = m_config->get_number("constants.ice.density"),
           cell_area = m_grid->cell_area();

    IceModelVec::AccessList list({ &tc_calving, &tc_frontal_melt, &tc_forced_retreat,
                                   &bmb, &mass_discharge, &mass_basal });

    //Update lake extend depending on exp_mask
    ParallelSection ParSec(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        //Convert from basal mass balance to basal mass loss [m].
        double bml_ij = -bmb(i, j);
        bml_ij = (bml_ij >= 0.0) ? bml_ij : 0.0;

        double discharge_ij = -1.0 * (tc_calving(i, j) + tc_frontal_melt(i, j) + tc_forced_retreat(i, j));
        discharge_ij = (discharge_ij >= 0.0) ? discharge_ij : 0.0;

        const double C = cell_area * rho_ice / dt;

        //Convert to mass flux [kg/s]
        mass_discharge(i, j) = C * discharge_ij;
        mass_basal(i, j)     = C * bml_ij;
      }
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();

    //Lake surface area
    IceModelVec2S cell_area2D(m_grid, "ceall_area", WITHOUT_GHOSTS);
    cell_area2D.set(cell_area);
    Lacc.accumulate(cell_area2D, lake_area);

    Lacc.accumulate(mass_discharge, lake_mass_input_discharge);
    Lacc.accumulate(mass_basal, lake_mass_input_basal);
  }

  {
    double rho_fresh_water = m_config->get_number("constants.fresh_water.density");

    IceModelVec::AccessList list({ &lake_area, &lake_mass_input_discharge,
                                   &lake_mass_input_basal, &lake_fill_rate });

    //Update lake extend depending on exp_mask
    ParallelSection ParSec(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();
        const double lake_area_ij = lake_area(i, j);
        if (m_gc.islake(lake_area_ij)){
          const double lake_mass_input_total_ij = lake_mass_input_discharge(i, j) + lake_mass_input_basal(i, j);
          lake_fill_rate(i, j) = lake_mass_input_total_ij / (rho_fresh_water * lake_area_ij);
        } else {
          lake_fill_rate(i, j) = m_fill_value;
        }
      }
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }
  lake_fill_rate.update_ghosts();
}

void LakeCC::updateLakeLevelMinMax(const IceModelVec2S &lake_level,
                                   const IceModelVec2S &target_level,
                                   IceModelVec2S &min_level,
                                   IceModelVec2S &max_level) {
  //Compute min max level
  ParallelSection ParSec(m_grid->com);
  try {
    // Initialze LakeProperties Model
    LakePropertiesCC LpCC(m_grid, m_fill_value, target_level, lake_level);
    LpCC.getLakeProperties(min_level, max_level);
  } catch (...) {
    ParSec.failed();
  }
  ParSec.check();
}

bool LakeCC::prepareLakeLevel(const IceModelVec2S &target_level,
                              const IceModelVec2S &bed,
                              const IceModelVec2S &thk,
                              const IceModelVec2S &min_level,
                              const IceModelVec2S &old_ll,
                              const IceModelVec2S &sea_level,
                              IceModelVec2S &min_basin,
                              IceModelVec2S &lake_level) {

  //Calculate basins that were filled by the ocean
  IceModelVec2S old_sl_basins(m_grid, "sl_basins", WITHOUT_GHOSTS);
  {
    IceModelVec::AccessList list({ &old_sl_basins, &sea_level, &bed, &thk });

    ParallelSection ParSec(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();
        if (mask::ocean(m_gc.mask(sea_level(i, j), bed(i, j), thk(i, j)))) {
          old_sl_basins(i, j) = sea_level(i, j);
        } else {
          old_sl_basins(i, j) = m_fill_value;
        }
      }
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();

    old_sl_basins.update_ghosts();
  }

  IceModelVec2Int exp_mask(m_grid, "exp_mask", WITHOUT_GHOSTS);
  IceModelVec2S max_sl_basin(m_grid, "max_sl_basin", WITHOUT_GHOSTS);

  { //Check which lake cells are newly added
    ParallelSection ParSec(m_grid->com);
    try {
      FilterExpansionCC FExCC(m_grid, m_fill_value, bed, old_sl_basins);
      FExCC.filter_ext(lake_level, target_level, exp_mask, min_basin, max_sl_basin);
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }

  bool MinMaxChanged = false;
  {
    IceModelVec::AccessList list({ &lake_level, &min_level, &min_basin,
                                  &exp_mask, &old_ll, &max_sl_basin });

    //Update lake extend depending on exp_mask
    ParallelSection ParSec(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        const int mask_ij = exp_mask.as_int(i, j);
        const bool new_lake = ( (mask_ij > 0) and not (m_gc.islake(min_level(i, j))) );
        if ( mask_ij == 1 or new_lake ) {
          //new lake or new basin added to existing lake
          if (max_sl_basin(i, j) != m_fill_value) {
            //basin was previously occupied by the ocean
            lake_level(i, j) = max_sl_basin(i, j);
          } else {
            //it was "dry" before
            if (new_lake) {
              //Initialize new basin with its lowest basin elevation
              lake_level(i, j) = min_basin(i, j);
            } else {
              //New basin added to existing lake -> set to min level
              lake_level(i, j) = std::min(min_level(i, j), min_basin(i, j));
            }
          }
          MinMaxChanged = true;
        } else if (mask_ij == 2) {
          //Extend existing lake by new cells
          if ( m_gc.islake(old_ll(i, j)) ) {
            lake_level(i, j) = old_ll(i, j);
          } else {
            lake_level(i, j) = min_level(i, j);
          }
        }
      }
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }

  lake_level.update_ghosts();

  return (GlobalOr(m_grid->com, MinMaxChanged));
}

void LakeCC::gradually_fill(const double dt,
                            const double max_fill_rate,
                            const IceModelVec2S &target_level,
                            const IceModelVec2S &bed,
                            const IceModelVec2S &thk,
                            const IceModelVec2S &sea_level,
                            const IceModelVec2S &min_level,
                            const IceModelVec2S &max_level,
                            const IceModelVec2S &min_bed,
                            const IceModelVec2S &fill_rate,
                            IceModelVec2S &lake_level) {

  double dh_max = max_fill_rate * dt;

  IceModelVec::AccessList list{ &lake_level, &target_level,
                                &min_level, &max_level, &min_bed,
                                &sea_level, &bed, &thk };

  if (not m_use_const_fill_rate) {
    list.add(fill_rate);
  }

  //Update lakes
  ParallelSection ParSec(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_gc.islake(lake_level(i, j))) {
        const double min_ij = min_level(i, j),
                     max_ij = max_level(i, j),
                     current_ij = lake_level(i, j),
                     target_level_ij  = target_level(i, j);
        const bool unbalanced_level = (min_ij != max_ij),
                   disappear = not m_gc.islake(target_level_ij);

        double target_ij = target_level_ij;

        // if lakes vanish ...
        if (disappear) {
          // ... and become ocean, set to sea level, else to floatation level
          if (mask::ocean(m_gc.mask(sea_level(i, j), bed(i, j), thk(i, j)))) {
            target_ij = sea_level(i, j);
          } else {
            if (mask::ocean(m_gc.mask(sea_level(i, j), bed(i, j), thk(i, j), current_ij))) {
              // if "real" lake, gradually empty it to floatation level
              target_ij = bed(i, j) + m_drho * thk(i, j);
            } else {
              //else (if it is below floatation level) immediately get rid of it.
              target_ij = current_ij;
            }
          }
        }

        const bool rising = (current_ij < target_ij);

        if (not m_use_const_fill_rate or not rising) {
          dh_max = max_fill_rate * dt;
        }

        if (rising) {
          //Only if rising, adjust fill rate to mass loss of ice sheet
          if (not (m_use_const_fill_rate or unbalanced_level)) {
            const double fill_rate_ij = fill_rate(i, j);

            if (m_gc.islake(fill_rate_ij) and (fill_rate_ij <= max_fill_rate) and not disappear) {
              dh_max = fill_rate_ij * dt;
            } else {
              dh_max = max_fill_rate * dt;
            }
          }

          const double new_level = ( disappear ? current_ij : min_ij ) + std::min(dh_max, std::abs(target_ij - min_ij));

          if ( (new_level > current_ij) or disappear ) {
            if (disappear and (new_level >= target_ij)) {
              lake_level(i, j) = m_fill_value;
            } else {
              lake_level(i, j) = new_level;
            }
          }
        } else {
          const double new_level = ( disappear ? current_ij : max_ij ) - std::min(dh_max, std::abs(max_ij - target_ij));

          if ( (new_level < current_ij) or disappear ) {
            if (disappear and (new_level <= target_ij)) {
              lake_level(i, j) = m_fill_value;
            } else {
              lake_level(i, j) = new_level;
            }
          }
        }
      }
    }
  } catch (...) {
    ParSec.failed();
  }
  ParSec.check();
}

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism

