// Copyright (C) 2008--2010 Ed Bueler and Constantine Khroulev
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

#include <cstring>
#include <cstdlib>
#include <petscda.h>
#include <netcdf.h>
#include "nc_util.hh"

#include "iceModelVec.hh"

// this file contains methods for derived classes IceModelVec2 and IceModelVec2Mask

// methods for base class IceModelVec are in "iceModelVec.cc"

IceModelVec2::IceModelVec2() : IceModelVec() {}

IceModelVec2::IceModelVec2(const IceModelVec2 &other) : IceModelVec(other) {}


PetscErrorCode  IceModelVec2::create(IceGrid &my_grid, const char my_name[], bool local) {
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }

  if (v != PETSC_NULL) {
    SETERRQ1(2,"IceModelVec2 with name='%s' already allocated\n", my_name);
  }
  PetscErrorCode ierr = create(my_grid, my_name, local, DA_STENCIL_BOX, 1); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode  IceModelVec2::create(IceGrid &my_grid, const char my_name[], bool local,
                                     DAStencilType my_sten, int stencil_width) {
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }
  if (v != PETSC_NULL) {
    SETERRQ1(2,"IceModelVec2 with name='%s' already allocated\n", my_name);
  }

  grid = &my_grid;
  dims = GRID_2D;
  
  PetscInt       M, N, m, n;
  PetscErrorCode ierr;
  ierr = DAGetInfo(my_grid.da2, PETSC_NULL, &N, &M, PETSC_NULL, &n, &m, PETSC_NULL,
                   PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);
  ierr = DACreate2d(my_grid.com, DA_XYPERIODIC, my_sten, N, M, n, m,
		    1,		// dof
		    stencil_width,
                    PETSC_NULL, PETSC_NULL, &da); CHKERRQ(ierr);

  if (local) {
    ierr = DACreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  localp = local;
  name = my_name;

  var1.init(my_name, my_grid, GRID_2D);

  //  ierr = this->set(GSL_NAN); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec2::get_array(PetscScalar** &a) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  a = (PetscScalar**) array;
  return 0;
}

//! Puts a local IceModelVec2 on processor 0.
/*!
 <ul>
 <li> onp0 and ctx should be created by calling VecScatterCreateToZero or be identical to one,
 <li> g2 is a preallocated temporary global vector,
 <li> g2natural is a preallocated temporary global vector with natural ordering.
 </ul>
*/
PetscErrorCode IceModelVec2::put_on_proc0(Vec onp0, VecScatter ctx, Vec g2, Vec g2natural) {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);

  if (!localp)
    SETERRQ1(1, "Can't put a global IceModelVec '%s' on proc 0.", name.c_str());

  ierr =        DALocalToGlobal(da, v,  INSERT_VALUES, g2);        CHKERRQ(ierr);
  ierr = DAGlobalToNaturalBegin(grid->da2, g2, INSERT_VALUES, g2natural); CHKERRQ(ierr);
  ierr =   DAGlobalToNaturalEnd(grid->da2, g2, INSERT_VALUES, g2natural); CHKERRQ(ierr);

  ierr = VecScatterBegin(ctx, g2natural, onp0, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);
  ierr =   VecScatterEnd(ctx, g2natural, onp0, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(ierr);

  return 0;
}

//! Gets a local IceModelVec2 from processor 0.
/*!
 <ul>
 <li> onp0 and ctx should be created by calling VecScatterCreateToZero or be identical to one,
 <li> g2 is a preallocated temporary global vector,
 <li> g2natural is a preallocated temporary global vector with natural ordering.
 </ul>
*/
PetscErrorCode IceModelVec2::get_from_proc0(Vec onp0, VecScatter ctx, Vec g2, Vec g2natural) {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);

  if (!localp)
    SETERRQ1(1, "Can't get a global IceModelVec '%s' from proc 0.", name.c_str());

  ierr = VecScatterBegin(ctx, onp0, g2natural, INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);
  ierr =   VecScatterEnd(ctx, onp0, g2natural, INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = DANaturalToGlobalBegin(grid->da2, g2natural, INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr =   DANaturalToGlobalEnd(grid->da2, g2natural, INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr =   DAGlobalToLocalBegin(da, g2,               INSERT_VALUES, v);  CHKERRQ(ierr);
  ierr =     DAGlobalToLocalEnd(da, g2,               INSERT_VALUES, v);  CHKERRQ(ierr);

  return 0;
}

//! Sets an IceModelVec2 to the magnitude of a 2D vector field with components \c v_x and \c v_y.
/*! Computes the magnitude \b pointwise, so any of v_x, v_y and the IceModelVec
  this is called on can be the same.

  Does not communicate.
 */
PetscErrorCode IceModelVec2::set_to_magnitude(IceModelVec2 &v_x, IceModelVec2 &v_y) {
  PetscErrorCode ierr;
  PetscScalar **mag = NULL;
  ierr = v_x.begin_access(); CHKERRQ(ierr);
  ierr = v_y.begin_access(); CHKERRQ(ierr);
  ierr = get_array(mag); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j) {
      mag[i][j] = sqrt( PetscSqr(v_x(i,j)) + PetscSqr(v_y(i,j)) );
    }
  }
  ierr = v_x.end_access(); CHKERRQ(ierr);
  ierr = v_y.end_access(); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);
  
  return 0;
}

//! Masks out all the areas where \f$ M \le 0 \f$ by setting them to \c fill. 
PetscErrorCode IceModelVec2::mask_by(IceModelVec2 &M, PetscScalar fill) {
  PetscErrorCode ierr;
  PetscScalar **a = NULL;
  ierr = get_array(a); CHKERRQ(ierr);
  ierr = M.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j) {
      if (M(i,j) <= 0.0)
	a[i][j] = fill;
    }
  }
  ierr = end_access(); CHKERRQ(ierr);
  ierr = M.end_access(); CHKERRQ(ierr);

  return 0;
}


//! \brief View a 2D field. Allocates and de-allocates g2, the temporary global
//! vector; performance should not matter here.
PetscErrorCode IceModelVec2::view(PetscInt viewer_size) {
  PetscErrorCode ierr;
  Vec g2;

  ierr = DACreateGlobalVector(grid->da2, &g2); CHKERRQ(ierr);
  
  if ((*map_viewers)[name] == PETSC_NULL) {
    string title = string_attr("long_name") + " (" + string_attr("glaciological_units") + ")";

    ierr = create_viewer(viewer_size, title, (*map_viewers)[name]); CHKERRQ(ierr);
  }

  if (localp) {
    ierr = copy_to_global(g2); CHKERRQ(ierr);
  } else {
    ierr = VecCopy(v, g2); CHKERRQ(ierr);
  }

  ierr = var1.to_glaciological_units(g2); CHKERRQ(ierr);

  ierr = VecView(g2, (*map_viewers)[name]); CHKERRQ(ierr);

  ierr = VecDestroy(g2); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec2::view_matlab(PetscViewer my_viewer) {
  PetscErrorCode ierr;
  string long_name = var1.get_string("long_name");
  Vec g2;

  ierr = DACreateGlobalVector(grid->da2, &g2); CHKERRQ(ierr);

  if (localp) {
    ierr = copy_to_global(g2); CHKERRQ(ierr);
  } else {
    ierr = VecCopy(v, g2); CHKERRQ(ierr);
  }

  ierr = var1.to_glaciological_units(g2); CHKERRQ(ierr);

  // add Matlab comment before listing, using short title

  ierr = PetscViewerASCIIPrintf(my_viewer, "\n%%%% %s = %s\n",
				name.c_str(), long_name.c_str()); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) g2, name.c_str()); CHKERRQ(ierr);

  ierr = VecView(g2, my_viewer); CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(my_viewer,"\n%s = reshape(%s,%d,%d);\n\n",
				name.c_str(), name.c_str(), grid->My, grid->Mx); CHKERRQ(ierr);

  ierr = VecDestroy(g2); CHKERRQ(ierr);

  return 0;
}



//! Provides access (both read and write) to the internal PetscScalar array.
/*!
  Note that i corresponds to the x direction and j to the y.
 */
inline
PetscScalar& IceModelVec2::operator() (int i, int j) {
  return static_cast<PetscScalar**>(array)[i][j];
}

///// IceModelVec2Mask

//! Returns the mask value; does not check ownership.
PismMask IceModelVec2Mask::value(int i, int j) {
  PetscScalar **a = (PetscScalar**) array;
  const PetscInt ival = static_cast<int>(floor(a[i][j] + 0.5));
  return static_cast<PismMask>(ival);
}


bool IceModelVec2Mask::is_grounded(int i, int j) {
  PismMask m = value(i, j);

  return (m == MASK_SHEET) || (m == MASK_DRAGGING_SHEET) || (m == MASK_ICE_FREE_BEDROCK);
}


bool IceModelVec2Mask::is_floating(int i, int j) {
  PismMask m = value(i, j);

  return (m == MASK_FLOATING) || (m == MASK_ICE_FREE_OCEAN) || (m == MASK_OCEAN_AT_TIME_0);
}

PetscErrorCode IceModelVec2Mask::fill_where_grounded(IceModelVec2 &v, const PetscScalar fillval) {
  PetscErrorCode ierr;

  ierr = begin_access(); CHKERRQ(ierr);
  ierr = v.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j) {
      if (is_grounded(i,j)) {
        v(i,j) = fillval;
      }
    }
  }
  ierr = v.end_access(); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);
  
  return 0;
}


PetscErrorCode IceModelVec2Mask::fill_where_floating(IceModelVec2 &v, const PetscScalar fillval) {
  PetscErrorCode ierr;

  ierr = begin_access(); CHKERRQ(ierr);
  ierr = v.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j) {
      if (is_floating(i,j)) {
        v(i,j) = fillval;
      }
    }
  }
  ierr = v.end_access(); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);
  
  return 0;
}

//! \brief Returns the x-derivative at i,j approximated using centered finite
//! differences.
inline
PetscScalar IceModelVec2::diff_x(int i, int j) {
  return ( (*this)(i + 1,j) - (*this)(i - 1,j) ) / (2 * grid->dx);
}

//! \brief Returns the y-derivative at i,j approximated using centered finite
//! differences.
inline
PetscScalar IceModelVec2::diff_y(int i, int j) {
  return ( (*this)(i,j + 1) - (*this)(i,j - 1) ) / (2 * grid->dy);
}

//! \brief Returns the x-derivative at i,j approximated using centered finite
//! differences. Respects grid periodicity and uses one-sided FD at grid edges
//! if necessary.
PetscScalar IceModelVec2::diff_x_p(int i, int j) {
  if (grid->periodicity & X_PERIODIC)
    return diff_x(i,j);
  
  if (i == 0)
    return ( (*this)(i + 1,j) - (*this)(i,j) ) / (grid->dx);
  else if (i == grid->Mx - 1)
    return ( (*this)(i,j) - (*this)(i - 1,j) ) / (grid->dx);
  else
    return diff_x(i,j);
}

//! \brief Returns the y-derivative at i,j approximated using centered finite
//! differences. Respects grid periodicity and uses one-sided FD at grid edges
//! if necessary.
PetscScalar IceModelVec2::diff_y_p(int i, int j) {
  if (grid->periodicity & Y_PERIODIC)
    return diff_y(i,j);
  
  if (j == 0)
    return ( (*this)(i,j + 1) - (*this)(i,j) ) / (grid->dy);
  else if (j == grid->My - 1)
    return ( (*this)(i,j) - (*this)(i,j - 1) ) / (grid->dy);
  else
    return diff_y(i,j);
}

//! Sums up all the values in an IceModelVec2 object. Ignores ghosts.
/*! Avoids copying to a "global" vector.
 */
PetscErrorCode IceModelVec2::sum(PetscScalar &result) {
  PetscErrorCode ierr;
  PetscScalar my_result = 0;

  // sum up the local part:
  ierr = begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (PetscInt j=grid->ys; j<grid->ys+grid->ym; ++j) {
      my_result += (*this)(i,j);
    }
  }
  ierr = end_access(); CHKERRQ(ierr);

  // find the global sum:
  ierr = PetscGlobalSum(&my_result, &result, grid->com); CHKERRQ(ierr);

  return 0;
}
