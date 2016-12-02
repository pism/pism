// Copyright (C) 2010, 2011, 2013, 2014, 2015, 2016 Constantine Khroulev
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

#include "PISMBedDef.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMTime.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMVars.hh"
#include "base/util/MaxTimestep.hh"

namespace pism {
namespace bed {

PBPointwiseIsostasy::PBPointwiseIsostasy(IceGrid::ConstPtr g)
  : BedDef(g) {
  m_thk_last.create(m_grid, "thk_last", WITH_GHOSTS, m_config->get_double("grid.max_stencil_width"));
}

PBPointwiseIsostasy::~PBPointwiseIsostasy() {
  // empty
}

void PBPointwiseIsostasy::init_impl() {

  m_log->message(2,
             "* Initializing the pointwise isostasy bed deformation model...\n");

  BedDef::init_impl();

  const IceModelVec2S *ice_thickness = m_grid->variables().get_2d_scalar("land_ice_thickness");
  m_thk_last.copy_from(*ice_thickness);
}

MaxTimestep PBPointwiseIsostasy::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("bed_def iso");
}

//! Updates the pointwise isostasy model.
void PBPointwiseIsostasy::update_with_thickness_impl(const IceModelVec2S &ice_thickness,
                                                     double my_t, double my_dt) {
  if ((fabs(my_t - m_t)   < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = my_t;
  m_dt = my_dt;

  double t_final = m_t + m_dt;

  // Check if it's time to update:
  double dt_beddef = t_final - m_t_beddef_last; // in seconds
  if ((dt_beddef < m_config->get_double("bed_deformation.update_interval", "seconds") &&
       t_final < m_grid->ctx()->time()->end()) ||
      dt_beddef < 1e-12) {
    return;
  }

  m_t_beddef_last = t_final;

  const double lithosphere_density = m_config->get_double("bed_deformation.lithosphere_density"),
    ice_density = m_config->get_double("constants.ice.density"),
    f = ice_density / lithosphere_density;

  //! Our goal: topg = topg_last - f*(thk - thk_last)

  //! Step 1: topg = topg_last - f*thk
  m_topg_last.add(-f, ice_thickness, m_topg);
  //! Step 2: topg = topg + f*thk_last = (topg_last - f*thk) + f*thk_last = topg_last - f*(thk - thk_last)
  m_topg.add(f, m_thk_last);
  //! This code is written this way to avoid allocating temp. storage for (thk - thk_last).

  //! Finally, we need to update bed uplift, topg_last and thk_last.
  compute_uplift(dt_beddef);

  m_thk_last.copy_from(ice_thickness);
  m_topg_last.copy_from(m_topg);

  //! Increment the topg state counter. SIAFD relies on this!
  m_topg.inc_state_counter();
}

} // end of namespace bed
} // end of namespace pism
