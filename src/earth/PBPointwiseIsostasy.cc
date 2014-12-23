// Copyright (C) 2010, 2011, 2013, 2014 Constantine Khroulev
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
#include "IceGrid.hh"
#include "PISMTime.hh"
#include "PISMConfig.hh"

#include <stdexcept>

namespace pism {

PBPointwiseIsostasy::PBPointwiseIsostasy(const IceGrid &g)
  : BedDef(g) {
  PetscErrorCode ierr;

  ierr = allocate();
  if (ierr != 0) {
    throw std::runtime_error("PBPointwiseIsostasy allocation failed");
  }

}

PetscErrorCode PBPointwiseIsostasy::allocate() {

  thk_last.create(m_grid, "thk_last", WITH_GHOSTS, m_config.get("grid_max_stencil_width"));

  return 0;
}

void PBPointwiseIsostasy::init() {

  BedDef::init();

  verbPrintf(2, m_grid.com,
             "* Initializing the pointwise isostasy bed deformation model...\n");

  thk->copy_to(thk_last);
  topg->copy_to(topg_last);
}

//! Updates the pointwise isostasy model.
void PBPointwiseIsostasy::update(double my_t, double my_dt) {
  if ((fabs(my_t - m_t)   < 1e-12) &&
      (fabs(my_dt - m_dt) < 1e-12)) {
    return;
  }

  m_t  = my_t;
  m_dt = my_dt;

  double t_final = m_t + m_dt;

  // Check if it's time to update:
  double dt_beddef = t_final - t_beddef_last; // in seconds
  if ((dt_beddef < m_config.get("bed_def_interval_years", "years", "seconds") &&
       t_final < m_grid.time->end()) ||
      dt_beddef < 1e-12) {
    return;
  }

  t_beddef_last = t_final;

  const double lithosphere_density = m_config.get("lithosphere_density"),
    ice_density = m_config.get("ice_density"),
    f = ice_density / lithosphere_density;

  //! Our goal: topg = topg_last - f*(thk - thk_last)

  //! Step 1: topg = topg_last - f*thk
  topg_last.add(-f, *thk, *topg);
  //! Step 2: topg = topg + f*thk_last = (topg_last - f*thk) + f*thk_last = topg_last - f*(thk - thk_last)
  topg->add(f, thk_last);
  //! This code is written this way to avoid allocating temp. storage for (thk - thk_last).

  //! Finally, we need to update bed uplift, topg_last and thk_last.
  compute_uplift(dt_beddef);

  thk->copy_to(thk_last);
  topg->copy_to(topg_last);

  //! Increment the topg state counter. SIAFD relies on this!
  topg->inc_state_counter();
}

} // end of namespace pism
