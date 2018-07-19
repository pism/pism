/* Copyright (C) 2018 PISM Authors
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

#include "pism/util/Time.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/Vars.hh"
#include "pism/geometry/Geometry.hh"

#include "Patch.hh"

namespace pism {
namespace ocean {
namespace lake_level {

Patch::Patch(IceGrid::ConstPtr grid,
             std::shared_ptr<LakeLevel> in)
  : LakeLevel(grid, in) {

  m_next_update_time          = m_grid->ctx()->time()->current();
  m_min_update_interval_years = 100;

  m_patch_iter = 3;

  m_last_update = 0.0;
}

Patch::~Patch() {
  //empty
}

void Patch::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  int patch_iter = m_patch_iter;
  patch_iter     = options::Integer("-lake_patch_iter",
                                    "Number of iterations used every time-step by patch-algorithm.", patch_iter);
  m_patch_iter   = patch_iter;
  m_log->message(2, "  LakeCC: number of iterations used by patch-algorithm: %d \n", m_patch_iter);

  int update_interval = m_min_update_interval_years;
  update_interval     = options::Integer("-lake_patch_update_interval",
                                         "Interval (in years) between updates of the lake mask.",
                                         update_interval);
  m_min_update_interval_years = update_interval;
  m_log->message(2, "  LakeCC: min update interval %d years \n", update_interval);

  m_lake_level.set(m_fill_value);

  m_next_update_time = m_grid->ctx()->time()->current();
}

void Patch::update_impl(const Geometry &geometry, double t, double dt) {

  bool full_update = false;
  if ((t >= m_next_update_time) or (fabs(t - m_next_update_time) < 1.0)) {
    full_update = true;
  }

  if (!full_update) {
    const IceModelVec2S &bed = geometry.bed_elevation,
                        &thk = geometry.ice_thickness,
                        &sl  = geometry.sea_level_elevation;

    for (unsigned int n = 0; n < m_patch_iter; ++n) {
      const int unsigned local_patch_result = patch_lake_levels(bed, thk, sl),
                         patch_result = GlobalMax(m_grid->com, local_patch_result);
      if (patch_result == 0) {
        //No further iteration needed
        break;
      }
      if ((patch_result == 2) or (n >= (m_patch_iter -1))) {
        //update needed
        full_update = true;
        break;
      }
    }
  }

  if (full_update) {
    m_input_model->update(geometry, t, dt);
    m_lake_level.copy_from(m_input_model->elevation());

    m_last_update = t;
    if (t != m_grid->ctx()->time()->start()) {
      m_next_update_time = m_grid->ctx()->time()->increment_date(t, m_min_update_interval_years);
    }
  }
}

MaxTimestep Patch::max_timestep_impl(double t) const {

  double dt = m_next_update_time - t;
  if (dt < 1.0) {
    double update_time_after_next = m_grid->ctx()->time()->increment_date(m_next_update_time,
                                                                          m_min_update_interval_years);
    dt = update_time_after_next - m_next_update_time;
    assert(dt > 0.0);
  }

  MaxTimestep lakecc_dt(dt, "lake level forcing");
  MaxTimestep input_max_timestep = m_input_model->max_timestep(t);

  if (input_max_timestep.finite()) {
    return std::min(input_max_timestep, lakecc_dt);
  } else {
    return lakecc_dt;
  }
}

unsigned int Patch::patch_lake_levels(const IceModelVec2S &bed,
                                      const IceModelVec2S &thk,
                                      const IceModelVec2S &sea_level) {
  const Direction dirs[] = { North, East, South, West };

  IceModelVec2S lake_level_old(m_grid, "ll", WITH_GHOSTS, 1);
  lake_level_old.copy_from(m_lake_level);

  IceModelVec::AccessList list{ &m_lake_level, &thk, &bed, &sea_level, &lake_level_old };

  unsigned int return_value = 0;
  GeometryCalculator gc(*m_config);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    const double bed_ij = bed(i, j),
                 thk_ij = thk(i, j),
                 sl_ij  = sea_level(i, j);
    if (not gc.islake(lake_level_old(i, j))) {
      StarStencil<double> ll_star = lake_level_old.star(i, j);

      double Level        = m_fill_value;
      unsigned int nLevel = 0;
      bool becomesLake    = false;

      for (int n = 0; n < 4; ++n) {
        const Direction direction = dirs[n];
        if (gc.islake(ll_star[direction])) {
          if (Level != ll_star[direction]) {
            nLevel++;
            Level = ll_star[direction];
            if (mask::ocean(gc.mask(m_fill_value, bed_ij, thk_ij, Level))) {
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
      if ((nLevel > 0) and mask::ocean(gc.mask(sl_ij, bed_ij, thk_ij))) {
        // Cell is now ocean and is next to lakes -> recalculate
        return 2;
      }
      if (becomesLake) {
        // Set cell to become a lake
        m_lake_level(i, j) = Level;
        return_value = 1;
      }
    } else { // cell was lake
      if (mask::ocean(gc.mask(sl_ij, bed_ij, thk_ij))) {
        // Was lake and is now Ocean -> full Update
        return 2;
      } else {
        if (not mask::ocean(gc.mask(m_fill_value, bed_ij, thk_ij, lake_level_old(i, j)))) {
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

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism
