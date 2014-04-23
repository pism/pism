// Copyright (C) 2012-2014 PISM Authors
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


#include "hydrology_diagnostics.hh"

namespace pism {

PISMHydrology_bwat::PISMHydrology_bwat(PISMHydrology *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiag<PISMHydrology>(m, g, my_vars) {
  vars[0].init_2d("bwat", grid);
  set_attrs("thickness of transportable water in subglacial layer", "", "m", "m", 0);
}

PetscErrorCode PISMHydrology_bwat::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "bwat", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;
  ierr = model->subglacial_water_thickness(*result); CHKERRQ(ierr);
  output = result;
  return 0;
}

PISMHydrology_bwp::PISMHydrology_bwp(PISMHydrology *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiag<PISMHydrology>(m, g, my_vars) {
  vars[0].init_2d("bwp", grid);
  set_attrs("pressure of transportable water in subglacial layer", "", "Pa", "Pa", 0);
}


PetscErrorCode PISMHydrology_bwp::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "bwp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;
  ierr = model->subglacial_water_pressure(*result); CHKERRQ(ierr);
  output = result;
  return 0;
}


PISMHydrology_bwprel::PISMHydrology_bwprel(PISMHydrology *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiag<PISMHydrology>(m, g, my_vars) {
  vars[0].init_2d("bwprel", grid);
  set_attrs("pressure of transportable water in subglacial layer as fraction of the overburden pressure", "",
            "", "", 0);
  vars[0].set_double("_FillValue", grid.config.get("fill_value"));
}


PetscErrorCode PISMHydrology_bwprel::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  double fill = grid.config.get("fill_value");
  IceModelVec2S *Po     = new IceModelVec2S,
                *result = new IceModelVec2S;
  ierr = result->create(grid, "bwprel", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  ierr = Po->create(grid, "Po_temporary", WITHOUT_GHOSTS); CHKERRQ(ierr);

  ierr = model->subglacial_water_pressure(*result); CHKERRQ(ierr);
  ierr = model->overburden_pressure(*Po); CHKERRQ(ierr);

  ierr = result->begin_access(); CHKERRQ(ierr);
  ierr = Po->begin_access(); CHKERRQ(ierr);
  for (int i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (int j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ((*Po)(i,j) > 0.0)
        (*result)(i,j) /= (*Po)(i,j);
      else
        (*result)(i,j) = fill;
    }
  }
  ierr = result->end_access(); CHKERRQ(ierr);
  ierr = Po->end_access(); CHKERRQ(ierr);

  output = result;
  return 0;
}


PISMHydrology_effbwp::PISMHydrology_effbwp(PISMHydrology *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiag<PISMHydrology>(m, g, my_vars) {
  vars[0].init_2d("effbwp", grid);
  set_attrs("effective pressure of transportable water in subglacial layer (overburden pressure minus water pressure)",
            "", "Pa", "Pa", 0);
}


PetscErrorCode PISMHydrology_effbwp::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *P      = new IceModelVec2S,
                *result = new IceModelVec2S;
  ierr = result->create(grid, "effbwp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  ierr = P->create(grid, "P_temporary", WITHOUT_GHOSTS); CHKERRQ(ierr);

  ierr = model->subglacial_water_pressure(*P); CHKERRQ(ierr);
  ierr = model->overburden_pressure(*result); CHKERRQ(ierr);
  ierr = result->add(-1.0,*P); CHKERRQ(ierr);  // result <-- result + (-1.0) P = Po - P

  output = result;
  return 0;
}


PISMHydrology_hydroinput::PISMHydrology_hydroinput(PISMHydrology *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiag<PISMHydrology>(m, g, my_vars) {
  vars[0].init_2d("hydroinput", grid);
  set_attrs("total water input into subglacial hydrology layer",
            "", "m s-1", "m/year", 0);
}


PetscErrorCode PISMHydrology_hydroinput::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "hydroinput", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;
  // the value reported diagnostically is merely the last value filled
  ierr = (model->total_input).copy_to(*result); CHKERRQ(ierr);
  output = result;
  return 0;
}


PISMHydrology_wallmelt::PISMHydrology_wallmelt(PISMHydrology *m, IceGrid &g, PISMVars &my_vars)
    : PISMDiag<PISMHydrology>(m, g, my_vars) {
  vars[0].init_2d("wallmelt", grid);
  set_attrs("wall melt into subglacial hydrology layer from (turbulent) dissipation of energy in transportable water",
            "", "m s-1", "m/year", 0);
}


PetscErrorCode PISMHydrology_wallmelt::compute(IceModelVec* &output) {
  PetscErrorCode ierr;
  IceModelVec2S *result = new IceModelVec2S;
  ierr = result->create(grid, "wallmelt", WITHOUT_GHOSTS); CHKERRQ(ierr);
  result->metadata() = vars[0];
  result->write_in_glaciological_units = true;
  ierr = model->wall_melt(*result); CHKERRQ(ierr);
  output = result;
  return 0;
}


} // end of namespace pism
