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

#include "pism/util/pism_options.hh"
#include "pism/util/Time.hh"
#include "pism/util/MaxTimestep.hh"

#include "SeaLevel2DCC.hh"
#include "SeaLevel_ConnectedComponents.hh"

namespace pism {
namespace ocean {
namespace sea_level {

SeaLevel2DCC::SeaLevel2DCC(IceGrid::ConstPtr g, std::shared_ptr<SeaLevel> in)
  : SeaLevel(g, in),
    m_target_level(m_grid, "target_sea_level", WITHOUT_GHOSTS),
    m_gc(*m_config),
    m_topg_overlay(m_grid, "topg_overlay", WITHOUT_GHOSTS),
    m_mask(m_grid, "sl_mask", WITHOUT_GHOSTS) {

  m_option = "sea_level.sl2dcc";

  m_max_update_interval_years = m_config->get_number(m_option + ".max_update_interval",
                                "years");
  if (m_max_update_interval_years < 0) {
    //Ful update every timestep requested
    m_max_update_interval_years = 0;
  }

  const double ice_density   = m_config->get_number("constants.ice.density"),
               ocean_density = m_config->get_number("constants.sea_water.density");
  m_density_ratio = ice_density/ocean_density;

  m_topg_overlay.set_attrs("internal",
                           "topography overlay",
                           "meter", "meter", "", 0);
  m_topg_overlay.set_time_independent(true);

  m_max_fill_rate = m_config->get_number(m_option + ".max_fill_rate", "meter second-1");

  m_gradual_active = true;
  std::string lake_models = m_config->get_string("lake_level.models");
  if (lake_models.find("lakecc") != std::string::npos) {
    //LakeCC model activated. -> Gradual model must be deactivated
    m_gradual_active = false;
  }

  m_sl_offset = m_config->get_number(m_option + ".sl_offset", "meter");
}

SeaLevel2DCC::~SeaLevel2DCC() {
  // empty
}

void SeaLevel2DCC::init_impl(const Geometry &geometry) {

  m_input_model->init(geometry);

  m_log->message(2, "  *Initializing SeaLevel2DCC model.\n");

  m_mask.set(0);

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

  m_log->message(3, "  SeaLevel2DCC: maximum interval of full update: %d years \n",
                 m_max_update_interval_years);

  double max_fill_rate_m_y = m_config->get_number(m_option + ".max_fill_rate", "meter year-1");
  m_log->message(3, "  SeaLevel2DCC: maximum fill rate: %g meter/year \n",
                 max_fill_rate_m_y);

  {
    InputOptions opts = process_input_options(m_grid->com, m_config);
    IceModelVec2S tmp(m_grid, "effective_sea_level_elevation", WITHOUT_GHOSTS);
    tmp.set_attrs("diagnostic",
                  "sea level elevation, relative to the geoid",
                  "meter", "meter", "", 0);
    tmp.metadata().set_number("_FillValue", m_fill_value);

    if (opts.type == INIT_RESTART) {

      m_log->message(3, "* Reading sea level forcing from '%s' for re-starting...\n",
                     opts.filename.c_str());

      File file(m_grid->com, opts.filename, PISM_GUESS, PISM_READONLY);
      const unsigned int time_length = file.nrecords(),
                         last_record = time_length > 0 ? time_length - 1 : 0;

      tmp.read(file, last_record);

      file.close();
    } else if (opts.type == INIT_BOOTSTRAP) {
      try {
        //effective_sea_level might be available in input file
        tmp.regrid(opts.filename, CRITICAL);
      } catch (...) {
        //if it was not found...
        tmp.copy_from(m_input_model->elevation());

        update_mask(geometry.bed_elevation,
                    geometry.ice_thickness,
                    tmp,
                    m_mask);

        apply_mask(m_mask,
                   tmp);
      }
    }
    m_sea_level.copy_from(tmp);
  }

  //Full update in first timestep
  m_next_update_time = m_grid->ctx()->time()->current();
}

void SeaLevel2DCC::update_impl(const Geometry &geometry, double t, double dt) {

  m_input_model->update(geometry, t, dt);
  m_target_level.copy_from(m_input_model->elevation());

  if (m_max_update_interval_years <= 0) {
    m_next_update_time = t;
  }

  //Check is a update of Mask is due!
  if ((t >= m_next_update_time) or (fabs(t - m_next_update_time) < 1.0)) {
    // Update sea level mask
    update_mask(geometry.bed_elevation,
                geometry.ice_thickness,
                m_target_level,
                m_mask);

    if (t != m_grid->ctx()->time()->start()) {
      //Reset next update time.
      m_next_update_time = m_grid->ctx()->time()->increment_date(t, m_max_update_interval_years);
    }
  }

  //Apply mask
  apply_mask(m_mask,
             m_target_level);

  if (m_gradual_active) {

    prepareSeaLevel(geometry.bed_elevation,
                    m_target_level,
                    m_sea_level);

    gradually_fill(dt,
                   geometry.bed_elevation,
                   geometry.ice_thickness,
                   m_target_level,
                   m_sea_level);
  } else {
    //LakeCC is active and takes care of this.
    m_sea_level.copy_from(m_target_level);
  }

}

void SeaLevel2DCC::update_mask(const IceModelVec2S &bed,
                               const IceModelVec2S &thk,
                               const IceModelVec2S &sea_level,
                               IceModelVec2Int &mask) {

  m_log->message(3, "->SL2dCC: Update of 2D Sea Level mask! \n");

  IceModelVec2S sl2dcc_bed(m_grid, "sl2dcc_bed", WITHOUT_GHOSTS);

  if (m_use_topg_overlay) {
    m_topg_overlay.add(1.0, bed, sl2dcc_bed);
  } else {
    sl2dcc_bed.copy_from(bed);
  }

  ParallelSection ParSec(m_grid->com);
  try {
    // Initialze LakeCC Model
    SeaLevelCC SLM(m_grid, m_density_ratio, sl2dcc_bed, thk, m_fill_value);
    SLM.computeMask(sea_level, m_sl_offset, mask);
  } catch (...) {
    ParSec.failed();
  }
  ParSec.check();

  mask.update_ghosts();

  m_log->message(3, "          Done!\n");
}

void SeaLevel2DCC::apply_mask(const IceModelVec2Int &mask,
                              IceModelVec2S &sea_level) {
  //Crop values according to mask
  IceModelVec::AccessList list( { &sea_level, &mask } );
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask(i, j) > 0) {
      sea_level(i, j) = m_fill_value;
    }
  }
  sea_level.update_ghosts();
}

void SeaLevel2DCC::prepareSeaLevel(const IceModelVec2S &bed,
                                   const IceModelVec2S &target_level,
                                   IceModelVec2S &sea_level) {

  IceModelVec2Int exp_mask(m_grid, "exp_mask", WITHOUT_GHOSTS);
  IceModelVec2S max_sl_basin(m_grid, "max_sl_basin", WITHOUT_GHOSTS),
                min_basin(m_grid, "min_basin", WITHOUT_GHOSTS);

  { //Check which lake cells are newly added
    ParallelSection ParSec(m_grid->com);
    try {
      //max_sl_basin will not contain any information, as sea_level was "m_fill_value" in the new basin
      //It should not be calculated in the future...
      FilterExpansionCC FExCC(m_grid, m_fill_value, bed, sea_level);
      FExCC.filter_ext(sea_level, target_level, exp_mask, min_basin, max_sl_basin);
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }

  {
    IceModelVec::AccessList list({ &sea_level, &target_level,
                                   &min_basin, &exp_mask });

    //Update lake extend depending on exp_mask
    ParallelSection ParSec(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        const int mask_ij = exp_mask.as_int(i, j);
        if (mask_ij > 0) {
          //Set the sea_level in the new basin either to the target level, if it is below the basin
          //or to the lowest basin elevation, to gradually increase the water..
          sea_level(i, j) = std::min(min_basin(i, j), target_level(i, j));
        }
      }
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }
  sea_level.update_ghosts();
}

void SeaLevel2DCC::gradually_fill(const double dt,
                                  const IceModelVec2S &bed,
                                  const IceModelVec2S &thk,
                                  const IceModelVec2S &target_level,
                                  IceModelVec2S &sea_level) {

  double dh_max = m_max_fill_rate * dt;

  IceModelVec::AccessList list({ &target_level, &sea_level, &bed, &thk });

  //Update SL
  ParallelSection ParSec(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      //If SL is not "invalid" -> prepareSeaLavel took care of this
      if (m_gc.islake(sea_level(i, j))) {
        const double current_ij = sea_level(i, j),
                     target_ij  = target_level(i, j);

        const bool rising = ((current_ij < target_ij) and m_gc.islake(target_ij));
        if (rising) {
          const double new_level = current_ij + std::min(dh_max, (target_ij - current_ij));
          if (new_level > current_ij) {
            sea_level(i, j) = new_level;
          }
        } else {
          const double dh_ij = m_gc.islake(target_ij) ? (current_ij - target_ij) : dh_max,
                       new_level = current_ij - std::min(dh_max, dh_ij);

          if ( (new_level < (bed(i, j) + m_density_ratio * thk(i, j))) and not m_gc.islake(target_ij) ) {
            //Ocean basin vanishes
            sea_level(i, j) = m_fill_value;
          } else {
            sea_level(i, j) = new_level;
          }
        }
      }

    } // end of loop
  } catch (...) {
    ParSec.failed();
  }
  ParSec.check();

  sea_level.update_ghosts();
}

MaxTimestep SeaLevel2DCC::max_timestep_impl(double t) const {

  if (m_max_update_interval_years > 0.0) {
    double dt = m_next_update_time - t;
    if (dt < 1.0) {
      double update_time_after_next = m_grid->ctx()->time()->increment_date(m_next_update_time,
                                                                            m_max_update_interval_years);
      dt = update_time_after_next - m_next_update_time;
      assert(dt > 0.0);
    }

    MaxTimestep sl2dcc_dt(dt, "sea level forcing");
    return sl2dcc_dt;
  } else {
    MaxTimestep sl2dcc_dt("sea level forcing");
    return sl2dcc_dt;
  }
}

bool SeaLevel2DCC::expandMargins_impl() const {
  return true;
}

DiagnosticList SeaLevel2DCC::diagnostics_impl() const {

  DiagnosticList result = {
    { "sl2dcc_gradual_target",       Diagnostic::wrap(m_target_level) },
  };

  return result;
}

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism
