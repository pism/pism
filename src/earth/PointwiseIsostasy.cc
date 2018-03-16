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
  m_thk_last.create(m_grid, "thk_last", WITH_GHOSTS, m_config->get_double("grid.max_stencil_width"));
}

PointwiseIsostasy::~PointwiseIsostasy() {
  // empty
}

void PointwiseIsostasy::init_impl(const InputOptions &opts) {

  m_log->message(2,
             "* Initializing the pointwise isostasy bed deformation model...\n");

  BedDef::init_impl(opts);

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  m_thk_last.copy_from(*ice_thickness);
}

MaxTimestep PointwiseIsostasy::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("bed_def iso");
}

//! Updates the pointwise isostasy model.
void PointwiseIsostasy::update_impl(const IceModelVec2S &ice_thickness,
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

  const double mantle_density = m_config->get_double("bed_deformation.mantle_density"),
    ice_density = m_config->get_double("constants.ice.density"),
    f = ice_density / mantle_density;

  //! Our goal: topg = topg_last - f*(thk - thk_last)

  //! Step 1: topg = topg_last - f*thk
  m_topg_last.add(-f, ice_thickness, m_topg);
  //! Step 2: topg = topg + f*thk_last = (topg_last - f*thk) + f*thk_last = topg_last - f*(thk - thk_last)
  m_topg.add(f, m_thk_last);
  //! This code is written this way to avoid allocating temp. storage for (thk - thk_last).

  //! Finally, we need to update bed uplift, topg_last and thk_last.
  compute_uplift(m_topg, m_topg_last, dt_beddef, m_uplift);

  m_thk_last.copy_from(ice_thickness);
  m_topg_last.copy_from(m_topg);

  //! Increment the topg state counter. SIAFD relies on this!
  m_topg.inc_state_counter();
}

} // end of namespace bed
} // end of namespace pism
