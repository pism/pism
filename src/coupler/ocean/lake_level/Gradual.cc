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

#include "pism/util/pism_options.hh"
#include "pism/util/Vars.hh"
#include "pism/geometry/Geometry.hh"

#include "Gradual.hh"
#include "LakeLevel_ConnectedComponents.hh"

namespace pism {
namespace ocean {
namespace lake_level {

Gradual::Gradual(IceGrid::ConstPtr grid,
                 std::shared_ptr<LakeLevel> in)
  : LakeLevel(grid, in) {

  //Set default filling rate for lakes
  m_max_lake_fill_rate = units::convert(m_sys, 10.0, "m/years", "m/seconds");
  m_target_level.create(m_grid, "target_level", WITHOUT_GHOSTS);
  m_target_level.set_attrs("model_state", "target lake level",
                           "m", "target_level");
  m_target_level.metadata().set_double("_FillValue", m_fill_value);

  m_min_level.create(m_grid, "min_lake_level", WITHOUT_GHOSTS);
  m_min_level.set_attrs("model_state", "min lake level",
                        "m", "min_level");
  m_min_level.metadata().set_double("_FillValue", m_fill_value);

  m_max_level.create(m_grid, "max_lake_level", WITHOUT_GHOSTS);
  m_max_level.set_attrs("model_state", "max lake level",
                        "m", "max_level");
  m_max_level.metadata().set_double("_FillValue", m_fill_value);

  m_min_bed.create(m_grid, "min_lake_bed", WITHOUT_GHOSTS);
  m_min_bed.set_attrs("model_state", "min lake bed",
                      "m", "min_bed");
  m_min_bed.metadata().set_double("_FillValue", m_fill_value);

  m_expansion_mask.create(m_grid, "expansion_mask", WITHOUT_GHOSTS);
  m_expansion_mask.metadata().set_double("_FillValue", m_fill_value);
}


void Gradual::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  {
    IceModelVec2S tmp;
    tmp.create(m_grid, "effective_lake_level_elevation", WITHOUT_GHOSTS);
    tmp.set_attrs("diagnostic",
                  "lake level elevation, relative to the geoid",
                  "meter", "");
    tmp.metadata().set_double("_FillValue", m_fill_value);

    InputOptions opts = process_input_options(m_grid->com, m_config);

    if (opts.type == INIT_RESTART) {

      m_log->message(2, "* Reading lake level forcing from '%s' for re-starting...\n",
                    opts.filename.c_str());

      PIO file(m_grid->com, "guess_mode", opts.filename, PISM_READONLY);
      const unsigned int time_length = file.inq_nrecords(),
                        last_record = time_length > 0 ? time_length - 1 : 0;

      tmp.read(file, last_record);

      file.close();
    } else {
      tmp.set(m_fill_value);
    }

    // Support regridding. This is needed to ensure that initialization using "-i" is
    // equivalent to "-i ... -bootstrap -regrid_file ..."
    {
      regrid("lake gradual filling modifier", tmp,
            REGRID_WITHOUT_REGRID_VARS);
    }

    m_lake_level.copy_from(tmp);
  }

  double max_lake_fill_rate_m_y = units::convert(m_sys, m_max_lake_fill_rate,
                                                 "m/seconds", "m/year");
  max_lake_fill_rate_m_y = options::Real("-lake_gradual_max_fill_rate",
                                         "Maximum rate at which lakes do fill (m/year)",
                                         max_lake_fill_rate_m_y);
  m_max_lake_fill_rate = units::convert(m_sys, max_lake_fill_rate_m_y, "m/year", "m/seconds");
}

void Gradual::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_target_level.copy_from(m_input_model->elevation());

  const IceModelVec2S &bed = geometry.bed_elevation,
                      &thk = geometry.ice_thickness,
                      &sl  = geometry.sea_level_elevation;

  { //Compute min max level and min bed...
    ParallelSection ParSec(m_grid->com);
    try {
      // Initialze LakeProperties Model
      LakePropertiesCC LpCC(m_grid, m_fill_value, m_target_level,
                            m_lake_level, bed);
      LpCC.getLakeProperties(m_min_level, m_max_level, m_min_bed);
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }

  prepareLakeLevel(m_target_level, bed, m_min_level, m_min_bed);

  gradually_fill(dt, m_target_level, bed, thk, sl, m_min_level, m_max_level, m_min_bed);
}


void Gradual::prepareLakeLevel(const IceModelVec2S &target_level,
                               const IceModelVec2S &bed,
                               const IceModelVec2S &min_level,
                               const IceModelVec2S &min_bed) {
  {
    ParallelSection ParSec(m_grid->com);
    try {
      FilterExpansionCC FExCC(m_grid, m_fill_value);
      FExCC.filter_ext(m_lake_level, target_level, m_expansion_mask);
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }

  {
    IceModelVec::AccessList list{ &m_lake_level, &min_level,
                                  &min_bed, &m_expansion_mask };

    GeometryCalculator gc(*m_config);

    //Update lake extend depending on exp_mask
    ParallelSection ParSec(m_grid->com);
    try {
      for (Points p(*m_grid); p; p.next()) {
        const int i = p.i(), j = p.j();

        const int mask_ij = m_expansion_mask.as_int(i, j);
        if (mask_ij == 1 or not gc.islake(min_level(i, j))) {
          //New lake basin
          m_lake_level(i, j) = min_bed(i, j);
        } else if (mask_ij == 2) {
          //Extend existing lake by new cells
          m_lake_level(i, j) = min_level(i, j);
        }
      }
    } catch (...) {
      ParSec.failed();
    }
    ParSec.check();
  }

  m_lake_level.update_ghosts();
}


void Gradual::gradually_fill(double dt,
                             const IceModelVec2S &target_level,
                             const IceModelVec2S &bed,
                             const IceModelVec2S &thk,
                             const IceModelVec2S &sea_level,
                             const IceModelVec2S &min_level,
                             const IceModelVec2S &max_level,
                             const IceModelVec2S &min_bed) {

  const double dh_max = m_max_lake_fill_rate * dt;

  GeometryCalculator gc(*m_config);

  IceModelVec::AccessList list{ &m_lake_level, &target_level,
                                &min_level, &max_level, &min_bed,
                                &sea_level, &bed, &thk };
  //Update lakes
  ParallelSection ParSec(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (gc.islake(m_lake_level(i, j))) {
        const double current_ij = m_lake_level(i, j),
                     target_ij  = target_level(i, j);

        const bool rising = ((current_ij < target_ij) and gc.islake(target_ij)) ? true : false;
        if (rising) {
          const double min_ij = min_level(i, j),
                       new_level = min_ij + std::min(dh_max, (target_ij - min_ij));
          if (new_level > current_ij) {
            m_lake_level(i, j) = new_level;
          }
        } else {
          const double max_ij = max_level(i, j),
                       dh_ij = gc.islake(target_ij) ? (current_ij - target_ij) : dh_max,
                       new_level = max_ij - std::min(dh_max, dh_ij);
          if (new_level > bed(i, j)
              and not mask::ocean(gc.mask(sea_level(i, j), bed(i, j), thk(i, j)))) {
            if (new_level < current_ij) {
              m_lake_level(i, j) = new_level;
            }
          } else {
            m_lake_level(i, j) = m_fill_value;
          }
        }
      }
    }
  } catch (...) {
    ParSec.failed();
  }
  ParSec.check();
}

// Write diagnostic variables to extra files if requested
DiagnosticList Gradual::diagnostics_impl() const {

  DiagnosticList result = {
    { "lake_gradual_target",    Diagnostic::wrap(m_target_level) },
    { "lake_gradual_min_level", Diagnostic::wrap(m_min_level) },
    { "lake_gradual_max_level", Diagnostic::wrap(m_max_level) },
    { "lake_gradual_min_bed",   Diagnostic::wrap(m_min_bed) },
    { "lake_expansion_mask",    Diagnostic::wrap(m_expansion_mask) },
  };

  return combine(result, m_input_model->diagnostics());
}

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism
