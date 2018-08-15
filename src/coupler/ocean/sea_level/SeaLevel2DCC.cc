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

#include "SeaLevel2DCC.hh"
#include "SeaLevel_ConnectedComponents.hh"

namespace pism {
namespace ocean {
namespace sea_level {

SeaLevel2DCC::SeaLevel2DCC(IceGrid::ConstPtr g, std::shared_ptr<SeaLevel> in)
  : SeaLevel(g, in) {

  m_log->message(2,
                 "  - Setting up SeaLevel2DCC Module...\n");

  m_option_prefix = "-ocean_sl2dcc";

  m_next_update_time = m_grid->ctx()->time()->current();
  m_update_interval_years = 100;

  m_mask.create(m_grid, "sl_mask", WITH_GHOSTS, 1);
  m_mask.set(0);

  m_update_periodic = false;
  m_update_passive  = false;
  m_update_startup  = false;

  m_offset = 0.0;

  const double ice_density = m_config->get_double("constants.ice.density"),
               ocean_density = m_config->get_double("constants.sea_water.density");
  m_drho = ice_density/ocean_density;
}

SeaLevel2DCC::~SeaLevel2DCC() {
  // empty
}


void SeaLevel2DCC::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);
  m_log->message(2,
                 "* Initializing the SL2dCC ocean model modifier...\n");

  process_options();

  m_next_update_time = m_grid->ctx()->time()->current();
}

void SeaLevel2DCC::process_options() {
  std::string default_update_scheme = "passive";
  options::Keyword update_scheme(m_option_prefix + "_update_scheme",
                                 "Specify the scheme how the lakecc model is updated",
                                 "passive,periodic",
                                 default_update_scheme);
  if (update_scheme == "periodic") {
    m_update_periodic = true;
  } else if (update_scheme == "passive") {
    m_update_passive = true;
  }
  m_log->message(2,
                 "  SL2dCC: use update scheme '%s'\n",
                 update_scheme.value().c_str());
  m_update = true;

  int update_interval = m_update_interval_years;
  update_interval = options::Integer(m_option_prefix + "_update_interval",
                                     "Interval (in years) between updates of the sea level mask.",
                                     update_interval);
  m_update_interval_years = update_interval;
  if (m_update_periodic) {
    m_log->message(2,
                   "  SL2dCC: min update interval %dyears \n",
                   update_interval);
  }

  double offset = m_offset;
  offset = options::Integer(m_option_prefix + "_offset_level",
                            "Level offset (m)",
                            offset);
  m_offset = offset;
  m_log->message(2,
                 "  SL2dCC: Level offset %gm \n",
                 offset);
}

void SeaLevel2DCC::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_sea_level.copy_from(m_input_model->elevation());

  if (m_update_periodic) {
    if (t >= m_next_update_time or fabs(t - m_next_update_time) < 1.0) {
      m_update = true;
    }
  }

  if (m_update) {
    const IceModelVec2S &bed = geometry.bed_elevation,
                        &thk = geometry.ice_thickness;

    // Update sea level mask
    do_sl_mask_update(bed, thk);
    m_next_update_time = m_grid->ctx()->time()->increment_date(t, m_update_interval_years);
    if (not m_update_passive and (t != m_grid->ctx()->time()->start())) {
      m_update = false;
    }
  }

  const Direction dirs[] = { North, East, South, West };

  //Crop values according to mask
  IceModelVec::AccessList list({&m_sea_level, &m_mask});
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_mask(i, j) > 0) {
      m_sea_level(i, j) = m_fill_value;
    }
  }
}

MaxTimestep SeaLevel2DCC::max_timestep_impl(double t) const {
  if (m_update_periodic) {
    double dt = m_next_update_time - t;

    if (dt < 1.0) {
      double update_time_after_next = m_grid->ctx()->time()->increment_date(m_next_update_time,
                                                                            m_update_interval_years);
      dt = update_time_after_next - m_next_update_time;
      assert(dt > 0.0);
    }

    MaxTimestep sl2dcc_dt(dt, "ocean sl2dcc");
    MaxTimestep input_max_timestep = m_input_model->max_timestep(t);

    if (input_max_timestep.finite()) {
      return std::min(input_max_timestep, sl2dcc_dt);
    } else {
      return sl2dcc_dt;
    }
  } else {
    return m_input_model->max_timestep(t);
  }
}

void SeaLevel2DCC::do_sl_mask_update(const IceModelVec2S &bed, const IceModelVec2S &thk) {

  m_log->message(3,
                 "->SL2dCC: Update of 2D Sea Level! \n");

  ParallelSection ParSec(m_grid->com);
  try {
    // Initialze LakeCC Model
    SeaLevelCC SLM(m_grid, m_drho, bed, thk, m_fill_value);
    SLM.computeMask(m_sea_level, m_offset, m_mask);
  } catch (...) {
    ParSec.failed();
  }
  ParSec.check();

  m_mask.update_ghosts();

  m_log->message(3,
                 "          Done!\n");
}

bool SeaLevel2DCC::expandMargins_impl() const {
  return true;
}

} // end of namespace sea_level
} // end of namespace ocean
} // end of namespace pism
