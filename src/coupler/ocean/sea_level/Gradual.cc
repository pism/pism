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

#include "pism/util/pism_options.hh"
#include "pism/util/Time.hh"
#include "pism/util/Vars.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/geometry/Geometry.hh"

#include "Gradual.hh"
#include "../lake_level/LakeLevel_ConnectedComponents.hh"

namespace pism {
namespace ocean {
namespace sea_level {

Gradual::Gradual(IceGrid::ConstPtr g, std::shared_ptr<SeaLevel> in)
  : SeaLevel(g, in) {

  m_log->message(2,
                 "  - Setting up SL Gradual Module...\n");

  m_max_fill_rate = units::convert(m_sys, 1.0, "m/years", "m/seconds");

  m_target_level.create(m_grid, "target_sea_level", WITHOUT_GHOSTS);
  m_target_level.set_attrs("model_state", "target sea level",
                           "m", "m", "target_sea_level", 0);
  m_target_level.metadata().set_number("_FillValue", m_fill_value);

  m_min_basin.create(m_grid, "min_sl_bed", WITHOUT_GHOSTS);
  m_min_basin.set_attrs("model_state", "min sl bed",
                      "m", "m", "min_sl_bed", 0);
  m_min_basin.metadata().set_number("_FillValue", m_fill_value);

  m_max_ll_basin.create(m_grid, "max_ll_basin", WITHOUT_GHOSTS);
  m_max_ll_basin.set_attrs("model_state", "max ll basin",
                           "m", "max_ll_basin");
  m_max_ll_basin.metadata().set_number("_FillValue", m_fill_value);

  m_expansion_mask.create(m_grid, "sl_expansion_mask", WITHOUT_GHOSTS);
  m_expansion_mask.metadata().set_number("_FillValue", m_fill_value);

  m_bootstrap = false;
}

Gradual::~Gradual() {
  // empty
}


void Gradual::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);
  m_log->message(2,
                 "* Initializing the Gradual ocean model modifier...\n");

  bool regridded = false;
  bool restarted = false;
  {
    IceModelVec2S tmp;
    tmp.create(m_grid, "effective_sea_level_elevation", WITHOUT_GHOSTS);
    tmp.set_attrs("diagnostic",
                  "sea level elevation, relative to the geoid",
                  "meter", "m", "", 0);
    tmp.metadata().set_number("_FillValue", m_fill_value);

    InputOptions opts = process_input_options(m_grid->com, m_config);

    if (opts.type == INIT_RESTART) {

      m_log->message(2, "* Reading sea level forcing from '%s' for re-starting...\n",
                    opts.filename.c_str());

      File file(m_grid->com, opts.filename, PISM_GUESS, PISM_READONLY);
      const unsigned int time_length = file.nrecords(),
                         last_record = time_length > 0 ? time_length - 1 : 0;

      tmp.read(file, last_record);

      file.close();

      restarted = true;
    } else {
      tmp.set(m_fill_value);
    }

    // Support regridding. This is needed to ensure that initialization using "-i" is
    // equivalent to "-i ... -bootstrap -regrid_file ..."
    {
      regrid("sea level gradual filling modifier", tmp,
             REGRID_WITHOUT_REGRID_VARS);
    }

    if (tmp.state_counter() == 2) {
      regridded = true;
    }

    m_sea_level.copy_from(tmp);
  }

  if (not restarted and not regridded) {
    m_bootstrap = true;
  }

  double max_fill_rate_m_y = units::convert(m_sys, m_max_fill_rate,
                                            "m/seconds", "m/year");
  max_fill_rate_m_y = options::Real("-sl_gradual_max_fill_rate",
                                    "Maximum rate at which SL is filled (m/year)",
                                    max_fill_rate_m_y);
  m_max_fill_rate = units::convert(m_sys, max_fill_rate_m_y, "m/year", "m/seconds");
}

void Gradual::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_target_level.copy_from(m_input_model->elevation());

  if (m_bootstrap) {
    m_sea_level.copy_from(m_target_level);
    m_bootstrap = false;
  }

  const IceModelVec2S &bed = geometry.bed_elevation,
                      &thk = geometry.ice_thickness,
                      &ll  = geometry.lake_level_elevation;

  prepareSeaLevel(m_target_level, bed, ll, m_expansion_mask, m_min_basin, m_max_ll_basin, m_sea_level);

  gradually_fill(dt, m_max_fill_rate, m_target_level, bed, thk, m_min_basin, m_sea_level);
}


void Gradual::prepareSeaLevel(const IceModelVec2S &target_level,
                              const IceModelVec2S &bed,
                              const IceModelVec2S &lake_level,
                              IceModelVec2Int &mask,
                              IceModelVec2S &min_basin,
                              IceModelVec2S &max_ll_basin,
                              IceModelVec2S &sea_level) {

  { //Check which lake cells are newly added
    ParallelSection ParSec(m_grid->com);
    try {
      FilterExpansionCC FExCC(m_grid, m_fill_value, bed, lake_level);
      FExCC.filter_ext2(sea_level, target_level, mask, min_basin, max_ll_basin);
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }

  {
    IceModelVec::AccessList list{ &sea_level, &min_basin, &mask, &target_level, &max_ll_basin };

    //Update lake extend depending on exp_mask
    ParallelSection ParSec(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        const int mask_ij = mask.as_int(i, j);
        if (mask_ij > 0) {
          //New basin
          if (max_ll_basin(i, j) == m_fill_value) {
            sea_level(i, j) = std::min(min_basin(i, j), target_level(i, j));
          } else {
            sea_level(i, j) = max_ll_basin(i, j);
          }
        } else if (mask_ij < 0) {
          sea_level(i, j) = m_fill_value;
        }
      }
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }

  sea_level.update_ghosts();
}


void Gradual::gradually_fill(const double dt,
                             const double max_fill_rate,
                             const IceModelVec2S &target_level,
                             const IceModelVec2S &bed,
                             const IceModelVec2S &thk,
                             const IceModelVec2S &min_bed,
                             IceModelVec2S &sea_level) {

  double dh_max = max_fill_rate * dt;

  GeometryCalculator gc(*m_config);

  IceModelVec::AccessList list{ &target_level, &min_bed,
                                &sea_level, &bed, &thk };

  //Update SL
  ParallelSection ParSec(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      //If SL is not m_fill_value
      if (gc.islake(sea_level(i, j))) {
        const double current_ij = sea_level(i, j),
                     target_ij  = target_level(i, j);

        const bool rising = ((current_ij < target_ij) and gc.islake(target_ij));
        if (rising) {
          const double new_level = current_ij + std::min(dh_max, (target_ij - current_ij));
          if (new_level > current_ij) {
            sea_level(i, j) = new_level;
          }
        } else {
          const double dh_ij = gc.islake(target_ij) ? (current_ij - target_ij) : dh_max,
                       new_level = current_ij - std::min(dh_max, dh_ij);

          if ( new_level < bed(i, j) and not gc.islake(target_ij) ) {
            sea_level(i, j) = m_fill_value;
          } else {
            if (new_level < current_ij) {
              sea_level(i, j) = new_level;
            }
          }
        }
      }

    } // end of loop
  } catch (...) {
    ParSec.failed();
  }
  ParSec.check();
}


// Write diagnostic variables to extra files if requested
DiagnosticList Gradual::diagnostics_impl() const {

  DiagnosticList result = {
    { "sl_gradual_target",       Diagnostic::wrap(m_target_level) },
    { "sl_gradual_min_bed",      Diagnostic::wrap(m_min_basin) },
    { "sl_gradual_max_ll_basin", Diagnostic::wrap(m_max_ll_basin) },
    { "sl_expansion_mask",       Diagnostic::wrap(m_expansion_mask) },
  };

  return combine(result, m_input_model->diagnostics());
}

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism
