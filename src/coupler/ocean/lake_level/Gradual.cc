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

  IceModelVec2S target_level;
  target_level.create(m_grid, "target_level", WITHOUT_GHOSTS);
  target_level.set_attrs("model_state", "target lake level",
                         "m", "target_level");

  target_level.copy_from(m_input_model->elevation());

  const IceModelVec2S *bed, *thk, *sl;
  IceModelVec2S tmp;
  if (geometry.bed_elevation.state_counter() == 0) {
    //Fake lake level timestep -> geometry.bed_elevation not available yet.
    //Try to get it from somewhere else
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

  gradually_fill(dt, target_level, *bed, *thk, *sl);
}

void Gradual::gradually_fill(double dt, const IceModelVec2S &target_level,
                             const IceModelVec2S &bed, const IceModelVec2S &thk,
                             const IceModelVec2S &sea_level) {

  const double dh_max = m_max_lake_fill_rate * dt;

  const double ice_density        = m_config->get_double("constants.ice.density"),
               freshwater_density = m_config->get_double("constants.fresh_water.density"),
               drho               = ice_density / freshwater_density;

  IceModelVec2S floating_thresh_level(m_grid, "floating_threshold_level", WITHOUT_GHOSTS);
  floating_thresh_level.copy_from(bed);
  floating_thresh_level.add(drho, thk);

  IceModelVec2S min_level(m_grid, "min_level", WITHOUT_GHOSTS),
                max_level(m_grid, "max_level", WITHOUT_GHOSTS),
                min_float_level_lake(m_grid, "min_float_level_lake", WITHOUT_GHOSTS);

  //Compute min max ...
  ParallelSection ParSec(m_grid->com);
  try {
    // Initialze LakeProperties Model
    LakePropertiesCC LpCC(m_grid, m_fill_value, target_level,
                          m_lake_level, floating_thresh_level);
    LpCC.getLakeProperties(min_level, max_level, min_float_level_lake);
  } catch (...) {
    ParSec.failed();
  }
  ParSec.check();

  GeometryCalculator gc(*m_config);

  IceModelVec::AccessList list{ &m_lake_level, &target_level,
                                &floating_thresh_level, &min_level, &max_level,
                                &min_float_level_lake, &sea_level, &bed, &thk };
  //Update lakes
  ParallelSection ParSec2(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      const bool part_of_lake_domain = (min_float_level_lake(i, j) != m_fill_value);

      if (part_of_lake_domain) {
        const double current_ij = m_lake_level(i, j),
                     target_ij  = target_level(i, j);
        bool rising;

        if (current_ij == m_fill_value) {
          //target_ij != m_fill_value
          rising = true;
        } else if (target_ij == m_fill_value) {
          //current_ij != m_fill_value
          rising = false;
        } else {
          rising = ((target_ij - current_ij) >= 0.0) ? true : false;
        }

        if (rising) {
          double min_ij = min_level(i, j);
          if (min_ij == m_fill_value) {
            //filling lake from very bottom
            min_ij = min_float_level_lake(i, j);
          }
          double dh_ij = target_ij - min_ij;
          dh_ij = std::min(dh_max, dh_ij);
          const double new_level = min_ij + dh_ij;
          if (new_level > floating_thresh_level(i, j)) {
            if ((new_level > current_ij) or (current_ij == m_fill_value)) {
              m_lake_level(i, j) = new_level;
            }
          } else {
            m_lake_level(i, j) = m_fill_value;
          }
        } else {
          double max_ij = max_level(i, j),
                 target_level = target_ij;
          if (target_level == m_fill_value) {
            target_level = floating_thresh_level(i, j);
          }
          double dh_ij = target_level - max_ij;
          dh_ij = std::min(dh_max, std::abs(dh_ij));
          const double new_level = max_ij - dh_ij;
          if (new_level > floating_thresh_level(i, j)
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
    ParSec2.failed();
  }
  ParSec2.check();
}

} // end of namespace lake_level
} // end of namespace ocean
} // end of namespace pism
