// Copyright (C) 2010, 2011, 2013, 2014, 2015, 2016, 2017, 2018 Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "BedDef.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Time.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Vars.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace bed {

PointwiseIsostasy::PointwiseIsostasy(IceGrid::ConstPtr g)
  : BedDef(g) {
  m_load_last.create(m_grid, "load_last", WITH_GHOSTS, m_config->get_double("grid.max_stencil_width"));
}

PointwiseIsostasy::~PointwiseIsostasy() {
  // empty
}

void PointwiseIsostasy::init_impl(const InputOptions &opts, const IceModelVec2S &ice_thickness,
                                  const IceModelVec2S &sea_level_elevation) {

  m_log->message(2,
                 "* Initializing the pointwise isostasy bed deformation model...\n");

  BedDef::init_impl(opts, ice_thickness, sea_level_elevation);

  // store the initial load
  compute_load(m_topg, ice_thickness, sea_level_elevation, m_load_last);
}


void PointwiseIsostasy::bootstrap_impl(const IceModelVec2S &bed_elevation,
                                       const IceModelVec2S &bed_uplift,
                                       const IceModelVec2S &ice_thickness,
                                       const IceModelVec2S &sea_level_elevation) {
  BedDef::bootstrap_impl(bed_elevation, bed_uplift, ice_thickness, sea_level_elevation);

  // store initial load and bed elevation
  compute_load(bed_elevation, ice_thickness, sea_level_elevation, m_load_last);
  m_topg_last.copy_from(bed_elevation);
}

MaxTimestep PointwiseIsostasy::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("bed_def iso");
}

//! Updates the pointwise isostasy model.
void PointwiseIsostasy::update_impl(const IceModelVec2S &ice_thickness,
                                    const IceModelVec2S &sea_level_elevation,
                                    double t, double dt) {

  double t_final = t + dt;

  // Check if it's time to update:
  double dt_beddef = t_final - m_t_beddef_last; // in seconds
  if ((dt_beddef < m_config->get_double("bed_deformation.update_interval", "seconds") &&
       t_final < m_grid->ctx()->time()->end()) ||
      dt_beddef < 1e-12) {
    return;
  }

  m_t_beddef_last = t_final;

  const double
    mantle_density = m_config->get_double("bed_deformation.mantle_density"),
    load_density   = m_config->get_double("constants.ice.density"),
    ocean_density  = m_config->get_double("constants.sea_water.density"),
    ice_density    = m_config->get_double("constants.ice.density"),
    f              = load_density / mantle_density;

  //! Our goal: topg = topg_last - f*(load - load_last)

  IceModelVec::AccessList list{&m_topg, &m_topg_last,
                               &ice_thickness, &sea_level_elevation, &m_load_last};

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double load = compute_load(m_topg(i, j),
                                 ice_thickness(i, j),
                                 sea_level_elevation(i, j),
                                 ice_density, ocean_density);

      m_topg(i, j) = m_topg_last(i, j) - f * (load - m_load_last(i, j));
      m_load_last(i, j) = load;
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  //! Finally, we need to update bed uplift, topg_last and load_last.
  compute_uplift(m_topg, m_topg_last, dt_beddef, m_uplift);

  m_topg_last.copy_from(m_topg);

  //! Increment the topg state counter. SIAFD relies on this!
  m_topg.inc_state_counter();
}

} // end of namespace bed
} // end of namespace pism
