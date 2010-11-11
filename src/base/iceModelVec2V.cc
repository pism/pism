// Copyright (C) 2009, 2010 Constantine Khroulev
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

#include "iceModelVec.hh"
#include "pism_const.hh"

IceModelVec2V::IceModelVec2V() : IceModelVec2() {
  dof = 2;
  vars.resize(dof);

  reset_attrs(0);
  reset_attrs(1);
  component_da = PETSC_NULL;
}

IceModelVec2V::~IceModelVec2V() {
  if (!shallow_copy) {
    destroy();
  }
}

PetscErrorCode IceModelVec2V::destroy() {
  PetscErrorCode ierr;

  ierr = IceModelVec::destroy(); CHKERRQ(ierr);

  if (component_da != PETSC_NULL) {
    ierr = DADestroy(component_da); CHKERRQ(ierr);
  }
  return 0;
}

PetscErrorCode  IceModelVec2V::create(IceGrid &my_grid, const char my_short_name[], bool local,
				      int stencil_width) {

  PetscErrorCode ierr = IceModelVec2::create(my_grid, my_short_name, local, DA_STENCIL_BOX,
					     stencil_width, dof); CHKERRQ(ierr);

  // component DA:
  ierr = DACreate2d(grid->com, DA_XYPERIODIC, DA_STENCIL_BOX,
		    grid->My, grid->Mx,
		    grid->Ny, grid->Nx,
		    1, stencil_width, // "1" is the dof of a component
                    grid->procs_y, grid->procs_x,
		    &component_da); CHKERRQ(ierr);

  string s_name = name;
  vars[0].init("u" + s_name, my_grid, GRID_2D);
  vars[1].init("v" + s_name, my_grid, GRID_2D);

  return 0;
}

PetscErrorCode IceModelVec2V::get_array(PISMVector2** &a) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  a = static_cast<PISMVector2**>(array);
  return 0;
}

PISMVector2& IceModelVec2V::operator()(int i, int j) {
  check_array_indices(i, j);
  return static_cast<PISMVector2**>(array)[i][j];
}


PetscErrorCode IceModelVec2V::read(const char filename[], const unsigned int time) {
  PetscErrorCode ierr;

  ierr = checkAllocated(); CHKERRQ(ierr);
  if (grid->da2 == PETSC_NULL)
    SETERRQ(1, "IceModelVec::read_from_netcdf: grid.da2 is NULL.");

  Vec tmp;			// a temporary one-component vector,
				// distributed across processors the same way v is
  ierr = DACreateGlobalVector(component_da, &tmp); CHKERRQ(ierr);

  for (int j = 0; j < dof; ++j) {
    // read the first component:
    ierr = vars[j].read(filename, time, tmp); CHKERRQ(ierr);
    ierr = IceModelVec2::set_component(j, tmp); CHKERRQ(ierr);
  }
  
  // The calls above only set the values owned by a processor, so we need to
  // communicate if localp == true:
  if (localp) {
    ierr = beginGhostComm(); CHKERRQ(ierr);
    ierr = endGhostComm(); CHKERRQ(ierr);
  }

  // Clean up:
  ierr = VecDestroy(tmp); CHKERRQ(ierr);
  return 0;  
}

PetscErrorCode IceModelVec2V::write(const char filename[], nc_type nctype) {
  PetscErrorCode ierr;

  ierr = checkAllocated(); CHKERRQ(ierr);
  if (grid->da2 == PETSC_NULL)
    SETERRQ(1, "IceModelVec::read_from_netcdf: grid.da2 is NULL.");

  Vec tmp;			// a temporary one-component vector,
				// distributed across processors the same way v is
  ierr = DACreateGlobalVector(component_da, &tmp); CHKERRQ(ierr);

  for (int j = 0; j < dof; ++j) {
    vars[j].time_independent = time_independent;
    ierr = IceModelVec2::get_component(j, tmp); CHKERRQ(ierr);
    ierr = vars[j].write(filename, nctype,
			 write_in_glaciological_units, tmp);
  }

  // Clean up:
  ierr = VecDestroy(tmp);
  return 0;
}

PetscErrorCode IceModelVec2V::magnitude(IceModelVec2S &result) {
  PetscErrorCode ierr;
  PISMVector2** a;
  PetscScalar **mag;

  ierr = result.get_array(mag); CHKERRQ(ierr);
  ierr = get_array(a);

  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j) {
      mag[i][j] = sqrt(PetscSqr(a[i][j].u) + PetscSqr(a[i][j].v));
    }
  }

  ierr = result.end_access(); CHKERRQ(ierr);
  ierr = end_access();
  return 0;
}

PetscErrorCode IceModelVec2V::regrid(const char filename[], LocalInterpCtx &lic, bool critical) {
  PetscErrorCode ierr;
  MaskInterp *m = NULL;

  Vec tmp;			// a temporary one-component vector,
				// distributed across processors the same way v is
  ierr = DACreateGlobalVector(component_da, &tmp); CHKERRQ(ierr);

  if (use_interpolation_mask) m = &interpolation_mask;

  for (int j = 0; j < dof; ++j) {
    ierr = vars[j].regrid(filename, lic, critical, false, 0.0, m, tmp); CHKERRQ(ierr);
    ierr = IceModelVec2::set_component(j, tmp); CHKERRQ(ierr);
  }

  // The calls above only set the values owned by a processor, so we need to
  // communicate if localp == true:
  if (localp) {
    ierr = beginGhostComm(); CHKERRQ(ierr);
    ierr = endGhostComm(); CHKERRQ(ierr);
  }

  // Clean up:
  ierr = VecDestroy(tmp);
  return 0;
}

PetscErrorCode IceModelVec2V::regrid(const char filename[], LocalInterpCtx &lic, PetscScalar default_value) {
  PetscErrorCode ierr;
  MaskInterp *m = NULL;

  Vec tmp;			// a temporary one-component vector,
				// distributed across processors the same way v is
  ierr = DACreateGlobalVector(component_da, &tmp); CHKERRQ(ierr);

  if (use_interpolation_mask) m = &interpolation_mask;

  for (int j = 0; j < dof; ++j) {
    ierr = vars[j].regrid(filename, lic, false, true, default_value, m, tmp); CHKERRQ(ierr);
    ierr = IceModelVec2::set_component(j, tmp); CHKERRQ(ierr);
  }

  // The calls above only set the values owned by a processor, so we need to
  // communicate if localp == true:
  if (localp) {
    ierr = beginGhostComm(); CHKERRQ(ierr);
    ierr = endGhostComm(); CHKERRQ(ierr);
  }

  // Clean up:
  ierr = VecDestroy(tmp);
  return 0;
}

bool IceModelVec2V::is_valid(PetscScalar U, PetscScalar V) {
  return vars[0].is_valid(U) && vars[1].is_valid(V);
}


PetscErrorCode IceModelVec2V::get_component(int n, IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = IceModelVec2::get_component(n, result.v); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec2V::set_component(int n, IceModelVec2S &source) {
  PetscErrorCode ierr;

  ierr = IceModelVec2::set_component(n, source.v); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec2V::set_name(const char new_name[], int /*component = 0*/) {
  string tmp = new_name;
  reset_attrs(0);
  reset_attrs(1);
  
  name = "vel" + tmp;

  vars[0].short_name = "u" + tmp;
  vars[1].short_name = "v" + tmp;

  return 0;
}

//! \brief View a 2D vector field. Allocates and de-allocates g2, the temporary global
//! vector; performance should not matter here.
PetscErrorCode IceModelVec2V::view(PetscInt viewer_size) {
  PetscErrorCode ierr;
  Vec g2;

  ierr = DACreateGlobalVector(grid->da2, &g2); CHKERRQ(ierr);

  string prefixes[2];
  prefixes[0] = "u";
  prefixes[1] = "v";
  
  for (int j = 0; j < dof; ++j) {
    string c_name = prefixes[j] + name;
    if ((*map_viewers)[c_name] == PETSC_NULL) {
      string title = string_attr("long_name", j) + " (" + string_attr("glaciological_units", j) + ")";

      ierr = create_viewer(viewer_size, title, (*map_viewers)[c_name]); CHKERRQ(ierr);
    }

    ierr = IceModelVec2::get_component(j, g2); CHKERRQ(ierr);

    ierr = vars[j].to_glaciological_units(g2); CHKERRQ(ierr);

    ierr = VecView(g2, (*map_viewers)[c_name]); CHKERRQ(ierr);
  }


  ierr = VecDestroy(g2); CHKERRQ(ierr);

  return 0;
}
