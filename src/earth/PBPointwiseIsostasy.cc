// Copyright (C) 2010, 2011 Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
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

PBPointwiseIsostasy::PBPointwiseIsostasy(IceGrid &g, const NCConfigVariable &conf)
  : PISMBedDef(g, conf) {
  PetscErrorCode ierr;

  ierr = allocate();
  if (ierr != 0) {
    PetscPrintf(grid.com, "PBPointwiseIsostasy::PBPointwiseIsostasy(...): allocate() failed\n");
    PISMEnd();
  }

}

PetscErrorCode PBPointwiseIsostasy::allocate() {
  PetscErrorCode ierr;
  PetscInt WIDE_STENCIL = grid.max_stencil_width;

  ierr = thk_last.create(grid, "thk_last", true, WIDE_STENCIL); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PBPointwiseIsostasy::init(PISMVars &vars) {
  PetscErrorCode ierr;

  ierr = PISMBedDef::init(vars); CHKERRQ(ierr);

  ierr = verbPrintf(2, grid.com,
		    "* Initializing the pointwise isostasy bed deformation model...\n"); CHKERRQ(ierr);

  ierr = thk->copy_to(thk_last);   CHKERRQ(ierr);
  ierr = topg->copy_to(topg_last); CHKERRQ(ierr);

  return 0;
}

//! Updates the pointwise isostasy model.
PetscErrorCode PBPointwiseIsostasy::update(PetscReal t_years, PetscReal dt_years) {
  PetscErrorCode ierr;

  if ((fabs(t_years - t)   < 1e-12) &&
      (fabs(dt_years - dt) < 1e-12))
    return 0;

  t  = t_years;
  dt = dt_years;

  // Check if it's time to update:
  PetscScalar dt_beddef = t_years - t_beddef_last;
  if (dt_beddef < config.get("bed_def_interval_years"))
    return 0;

  t_beddef_last = t_years;

  const PetscScalar lithosphere_density = config.get("lithosphere_density"),
    ice_density = config.get("ice_density"),
    f = ice_density / lithosphere_density;

  //! Our goal: topg = topg_last - f*(thk - thk_last)

  //! Step 1: topg = topg_last - f*thk
  ierr = topg_last.add(-f, *thk, *topg); CHKERRQ(ierr);
  //! Step 2: topg = topg + f*thk_last = (topg_last - f*thk) + f*thk_last = topg_last - f*(thk - thk_last)
  ierr = topg->add(f, thk_last); CHKERRQ(ierr);
  //! This code is written this way to avoid allocating temp. storage for (thk - thk_last).

  //! Finally, we need to update bed uplift, topg_last and thk_last.
  ierr = compute_uplift(dt_beddef); CHKERRQ(ierr);

  ierr =  thk->copy_to(thk_last);  CHKERRQ(ierr);
  ierr = topg->copy_to(topg_last); CHKERRQ(ierr);

  return 0;
}
