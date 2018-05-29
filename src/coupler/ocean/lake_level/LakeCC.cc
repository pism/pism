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
#include <vector>

#include "pism/util/Time.hh"
#include "pism/util/Vars.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/geometry/Geometry.hh"

#include "LakeCC.hh"
#include "LakeLevel_ConnectedComponents.hh"

namespace pism {
namespace ocean {
namespace lake_level {

LakeCC::LakeCC(IceGrid::ConstPtr g)
  : LakeLevel(g),
    m_gc(*m_config) {

  m_log->message(2, "  - Setting up LakeCC Module...\n");

  m_option = "-lakecc";

  m_next_update_time          = m_grid->ctx()->time()->current();
  m_min_update_interval_years = 100;

  m_lake_level_min = 0.;
  m_lake_level_max = 1000.;
  m_lake_level_dh  = 10.;

  m_update_passive  = false;
  m_update_patch    = false;
  m_update_periodic = false;

  m_patch_iter = 3;

  m_icefree_thickness = m_config->get_double("geometry.ice_free_thickness_standard");

  const double ice_density        = m_config->get_double("constants.ice.density"),
               freshwater_density = m_config->get_double("constants.fresh_water.density");
  m_drho = ice_density / freshwater_density;

  m_last_update = 0.0;

}

LakeCC::~LakeCC() {
  // empty
}

void LakeCC::init_impl(const Geometry &geometry) {

  m_log->message(2, "* Initializing the lakeCC ocean model modifier...\n");

  m_lake_level.set(m_fill_value);

  process_options();

  m_next_update_time = m_grid->ctx()->time()->current();

}

void LakeCC::process_options() {
  std::string default_update_scheme = "periodic";
  options::Keyword update_scheme(m_option + "_update_scheme", "Specify the scheme how the lakecc model is updated",
                                 "passive,patch,periodic", default_update_scheme);
  if (update_scheme == "patch") {
    m_update_patch    = true;
    m_update_periodic = true;
  } else if (update_scheme == "periodic") {
    m_update_periodic = true;
  } else if (update_scheme == "passive") {
    m_update_passive  = true;
  }
  m_log->message(2, "  LakeCC: use update scheme '%s'\n", update_scheme.value().c_str());
  m_update = true;

  int update_interval = m_min_update_interval_years;
  update_interval     = options::Integer(m_option + "_update_interval",
                                         "Interval (in years) between updates of the lake mask.",
                                         update_interval);
  m_min_update_interval_years = update_interval;
  if (m_update_periodic) {
    m_log->message(2, "  LakeCC: min update interval %d years \n", update_interval);
  }

  std::vector<double> default_levels;
  default_levels.push_back(m_lake_level_min);
  default_levels.push_back(m_lake_level_max);
  default_levels.push_back(m_lake_level_dh);
  options::RealList levels(m_option + "_level", "Comma-seperated list of three numbers: z_min, z_max, dz",
                           default_levels);
  m_lake_level_min = levels[0];
  m_lake_level_max = levels[1];
  m_lake_level_dh  = levels[2];

  m_log->message(2, "  LakeCC: lake levels between %g and %gm, with %gm spacing\n",
                 m_lake_level_min, m_lake_level_max, m_lake_level_dh);

  int patch_iter = m_patch_iter;
  patch_iter     = options::Integer(m_option + "_patch_iter",
                                    "Number of iterations used every time-step by patch-algorithm.", patch_iter);
  m_patch_iter   = patch_iter;
  if (m_update_patch) {
    m_log->message(2, "  LakeCC: number of iterations used by patch-algorithm: %d \n", m_patch_iter);
  }

  double icefree_thickness = m_icefree_thickness;
  icefree_thickness        = options::Real(m_option + "_icefree_thickness",
                                           "Minimum ice thickness that is regarded as a closed ice cover.", icefree_thickness);
  m_icefree_thickness      = icefree_thickness;
  m_log->message(2, "  LakeCC: icefree thickness: %gm \n", m_icefree_thickness);
}

void LakeCC::update_impl(const Geometry &geometry, double my_t, double my_dt) {

  bool init = false;
  const IceModelVec2S *bed, *thk, *sl;
  IceModelVec2S tmp;
  if (geometry.bed_elevation.state_counter() < m_grid->ctx()->size()) {
    //Fake sea level timestep -> geometry.bed_elevation not available yet.
    //Try to get it from somewhere else
    init = true;

    tmp.create(m_grid, "topg", WITHOUT_GHOSTS);
    tmp.set_attrs("model_state", "bedrock surface elevation",
                  "m", "bedrock_altitude");

    //Try to initialize topg from file
    InputOptions opts = process_input_options(m_grid->com, m_config);
    if (opts.type != INIT_RESTART) {
      // bootstrapping
      tmp.regrid(opts.filename, OPTIONAL,
                    m_config->get_double("bootstrapping.defaults.bed"));
    }
    bed = &tmp;
  } else {
    bed = &(geometry.bed_elevation);
  }
  thk = &(geometry.ice_thickness);

  if (geometry.sea_level_elevation.state_counter() < m_grid->ctx()->size()) {
    sl = m_grid->variables().get_2d_scalar("sea_level");
  } else {
    sl = &(geometry.sea_level_elevation);
  }

  if (m_update_periodic) {
    if ((my_t >= m_next_update_time) or (fabs(my_t - m_next_update_time) < 1.0)) {
      m_update = true;
    }
  }
  if (m_update_patch and not m_update) {
    //Save copy of original lake levels
    IceModelVec2S ll_old;
    ll_old.create(m_grid, "lake_levels_old", WITHOUT_GHOSTS);
    ll_old.copy_from(m_lake_level);

    for (unsigned int n=0; n < m_patch_iter; ++n) {
      const int unsigned local_patch_result = patch_lake_levels(bed, thk, sl),
                         patch_result = GlobalMax(m_grid->com, local_patch_result);
      if (patch_result == 0) {
        //No further iteration needed
        break;
      }
      if ((patch_result == 2) or (n == (m_patch_iter -1))) {
        //update needed
        m_update = true;
        break;
      }
    }

    if (m_update) {
      //restore lake levels
      m_lake_level.copy_from(ll_old);
    }
  }

  if (m_update) {
    // Update Lakes
    do_lake_update(bed, thk, sl);

    m_last_update      = my_t;
    m_next_update_time = m_grid->ctx()->time()->increment_date(my_t, m_min_update_interval_years);
    if (not m_update_passive and not init) {
      m_update = false;
    }
  }

}

MaxTimestep LakeCC::max_timestep_impl(double t) const {
  if (m_update_periodic) {
    double dt = m_next_update_time - t;
    if (dt < 1.0) {
      double update_time_after_next = m_grid->ctx()->time()->increment_date(m_next_update_time,
                                                                            m_min_update_interval_years);
      dt = update_time_after_next - m_next_update_time;
      assert(dt > 0.0);
    }

    MaxTimestep lakecc_dt(dt, "lake level forcing");

    return lakecc_dt;
  } else {
    return MaxTimestep("lake level forcing");
  }
}

unsigned int LakeCC::patch_lake_levels(const IceModelVec2S *bed, const IceModelVec2S *thk, const IceModelVec2S *sea_level) {
  const Direction dirs[] = { North, East, South, West };

  IceModelVec2S lake_level_old(m_grid, "ll", WITH_GHOSTS, 1);
  IceModelVec::AccessList list{ &m_lake_level, thk, bed, sea_level, &lake_level_old };
  m_lake_level.update_ghosts();
  lake_level_old.copy_from(m_lake_level);

  unsigned int return_value = 0;

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double bed_ij = (*bed)(i, j),
                 thk_ij = (*thk)(i, j),
                 sl_ij  = (*sea_level)(i, j);
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
              // Lake cell adjoins several lakes -> recalculate
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
        m_lake_level(i, j) = Level;
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
          m_lake_level(i, j) = m_fill_value;
          return_value = 1;
        }
      }
    }
  }
  // If function gets here all cells have been checked and patched
  return return_value;
}

void LakeCC::prepare_mask_validity(const IceModelVec2S *thk, IceModelVec2Int& valid_mask) {
  IsolationCC IsoCC(m_grid, *thk, m_icefree_thickness);
  IsoCC.find_isolated_spots();
  IsoCC.isolation_mask(valid_mask);

  IceModelVec::AccessList list{ &valid_mask, &m_lake_level };
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    //Set sink, where pism_mask is ocean or at margin of computational domain
    if (m_gc.islake(m_lake_level(i, j))) {
      valid_mask(i, j) = 1;
    }
  }
  valid_mask.update_ghosts();
}

void LakeCC::do_lake_update(const IceModelVec2S *bed, const IceModelVec2S *thk, const IceModelVec2S *sea_level) {
  IceModelVec2Int pism_mask, valid_mask;
  pism_mask.create(m_grid, "mask", WITHOUT_GHOSTS);
  valid_mask.create(m_grid, "mask", WITHOUT_GHOSTS);
  m_gc.compute_mask(*sea_level, *bed, *thk, pism_mask);
  prepare_mask_validity(thk, valid_mask);

  m_log->message(2, "->LakeCC: Update of Lake Levels! \n");

  ParallelSection ParSec(m_grid->com);
  try {
    // Initialze LakeCC Model
    LakeLevelCC LM(m_grid, m_drho, m_icefree_thickness, *bed, *thk, pism_mask, m_fill_value, valid_mask);
    LM.floodMap(m_lake_level_min, m_lake_level_max, m_lake_level_dh);
    LM.lake_levels(m_lake_level);
  } catch (...) {
    ParSec.failed();
  }
  ParSec.check();

  m_log->message(2, "          Done!\n");
}

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism

