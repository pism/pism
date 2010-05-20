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
#include "nc_util.hh"

#include "iceModelVec.hh"

// this file contains methods for derived classes IceModelVec2S and IceModelVec2Mask

// methods for base class IceModelVec are in "iceModelVec.cc"

PetscErrorCode  IceModelVec2S::create(IceGrid &my_grid, const char my_name[], bool local, int width) {
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }

  if (v != PETSC_NULL) {
    SETERRQ1(2,"IceModelVec2S with name='%s' already allocated\n", my_name);
  }
  PetscErrorCode ierr = IceModelVec2::create(my_grid, my_name, local, DA_STENCIL_BOX, width, dof); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IceModelVec2S::get_array(PetscScalar** &a) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  a = static_cast<PetscScalar**>(array);
  return 0;
}

//! Puts a local IceModelVec2S on processor 0.
/*!
 <ul>
 <li> onp0 and ctx should be created by calling VecScatterCreateToZero or be identical to one,
 <li> g2 is a preallocated temporary global vector,
 <li> g2natural is a preallocated temporary global vector with natural ordering.
 </ul>
*/
PetscErrorCode IceModelVec2S::put_on_proc0(Vec onp0, VecScatter ctx, Vec g2, Vec g2natural) {
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
PetscErrorCode IceModelVec2S::get_from_proc0(Vec onp0, VecScatter ctx, Vec g2, Vec g2natural) {
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
PetscErrorCode IceModelVec2S::set_to_magnitude(IceModelVec2S &v_x, IceModelVec2S &v_y) {
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
PetscErrorCode IceModelVec2S::mask_by(IceModelVec2S &M, PetscScalar fill) {
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

  if (dof != 1)
    SETERRQ(1, "This method only supports IceModelVecs with dof == 1.");

  ierr = DACreateGlobalVector(grid->da2, &g2); CHKERRQ(ierr);
  
  if ((*map_viewers)[name] == PETSC_NULL) {
    string title = string_attr("long_name") + " (" + string_attr("glaciological_units") + ")";

    ierr = create_viewer(viewer_size, title, (*map_viewers)[name]); CHKERRQ(ierr);
  }

  if (localp) {
    ierr = copy_to(g2); CHKERRQ(ierr);
  } else {
    ierr = VecCopy(v, g2); CHKERRQ(ierr);
  }

  ierr = vars[0].to_glaciological_units(g2); CHKERRQ(ierr);

  ierr = VecView(g2, (*map_viewers)[name]); CHKERRQ(ierr);

  ierr = VecDestroy(g2); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec2S::view_matlab(PetscViewer my_viewer) {
  PetscErrorCode ierr;
  string long_name = vars[0].get_string("long_name");
  Vec g2;

  ierr = DACreateGlobalVector(grid->da2, &g2); CHKERRQ(ierr);

  if (localp) {
    ierr = copy_to(g2); CHKERRQ(ierr);
  } else {
    ierr = VecCopy(v, g2); CHKERRQ(ierr);
  }

  ierr = vars[0].to_glaciological_units(g2); CHKERRQ(ierr);

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
PetscScalar& IceModelVec2S::operator() (int i, int j) {
  return static_cast<PetscScalar**>(array)[i][j];
}


//! Checks if the current IceModelVec2S has NANs and reports if it does.
/*! Up to a fixed number of messages are printed at stdout.  Returns the full
 count of NANs (which is a nonzero) on this rank. */
PetscErrorCode IceModelVec2S::has_nan() {
  PetscErrorCode ierr;
  const PetscInt max_print_this_rank=10;
  PetscInt retval=0;

  ierr = begin_access(); CHKERRQ(ierr);
  PetscInt i, j;
  for (i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (j=grid->ys; j<grid->ys+grid->ym; ++j) {
      if (gsl_isnan((*this)(i,j))) {
        retval++;
        if (retval <= max_print_this_rank) {
          ierr = PetscSynchronizedPrintf(grid->com, 
             "IceModelVec2S %s: NAN (or uninitialized) at i = %d, j = %d on rank = %d\n",
             name.c_str(), i, j, grid->rank); CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = end_access(); CHKERRQ(ierr);

  if (retval > 0) {
    ierr = PetscSynchronizedPrintf(grid->com, 
       "IceModelVec2S %s: detected %d NANs (or uninitialized) on rank = %d\n",
             name.c_str(), retval, grid->rank); CHKERRQ(ierr);
  }

  ierr = PetscSynchronizedFlush(grid->com); CHKERRQ(ierr);
  return retval;
}


// IceModelVec2Mask

PetscErrorCode  IceModelVec2Mask::create(IceGrid &my_grid, const char my_name[], bool local, int width) {
  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }

  if (v != PETSC_NULL) {
    SETERRQ1(2,"IceModelVec2Mask with name='%s' already allocated\n", my_name);
  }
  PetscErrorCode ierr = IceModelVec2::create(my_grid, my_name, local, DA_STENCIL_BOX, width, dof); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IceModelVec2Mask::get_array(PetscScalar** &a) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  a = (PetscScalar**) array;
  return 0;
}

//! Returns the mask value; does not check ownership.
PismMask IceModelVec2Mask::value(int i, int j) {
  const PetscScalar **a = (const PetscScalar**) array;
  const PetscInt ival = static_cast<int>(floor(a[i][j] + 0.5));
  return static_cast<PismMask>(ival);
}

PetscScalar& IceModelVec2Mask::operator() (int i, int j) {
  return static_cast<PetscScalar**>(array)[i][j];
}

bool IceModelVec2Mask::is_grounded(int i, int j) {
  const PismMask m = value(i, j);

  return (m == MASK_SHEET) || (m == MASK_DRAGGING_SHEET) || (m == MASK_ICE_FREE_BEDROCK);
}


bool IceModelVec2Mask::is_floating(int i, int j) {
  const PismMask m = value(i, j);

  return (m == MASK_FLOATING) || (m == MASK_ICE_FREE_OCEAN) || (m == MASK_OCEAN_AT_TIME_0);
}

PetscErrorCode IceModelVec2Mask::fill_where_grounded(IceModelVec2S &v, const PetscScalar fillval) {
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


PetscErrorCode IceModelVec2Mask::fill_where_floating(IceModelVec2S &v, const PetscScalar fillval) {
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
PetscScalar IceModelVec2S::diff_x(int i, int j) {
  return ( (*this)(i + 1,j) - (*this)(i - 1,j) ) / (2 * grid->dx);
}

//! \brief Returns the y-derivative at i,j approximated using centered finite
//! differences.
PetscScalar IceModelVec2S::diff_y(int i, int j) {
  return ( (*this)(i,j + 1) - (*this)(i,j - 1) ) / (2 * grid->dy);
}


//! \brief Returns the x-derivative at East staggered point i+1/2,j approximated 
//! using centered (obvious) finite differences.
PetscScalar IceModelVec2S::diff_x_stagE(int i, int j) {
  return ( (*this)(i+1,j) - (*this)(i,j) ) / (grid->dx);
}

//! \brief Returns the y-derivative at East staggered point i+1/2,j approximated 
//! using centered [\ref Mahaffy] finite differences.
PetscScalar IceModelVec2S::diff_y_stagE(int i, int j) {
  return (   (*this)(i+1,j+1) + (*this)(i,j+1)
           - (*this)(i+1,j-1) - (*this)(i,j-1) ) / ( 4* grid->dy);
}

//! \brief Returns the x-derivative at North staggered point i,j+1/2 approximated 
//! using centered [\ref Mahaffy] finite differences.
PetscScalar IceModelVec2S::diff_x_stagN(int i, int j) {
  return (   (*this)(i+1,j+1) + (*this)(i+1,j)
           - (*this)(i-1,j+1) - (*this)(i-1,j) ) / ( 4* grid->dx);
}

//! \brief Returns the y-derivative at North staggered point i,j+1/2 approximated 
//! using centered (obvious) finite differences.
PetscScalar IceModelVec2S::diff_y_stagN(int i, int j) {
  return ( (*this)(i,j+1) - (*this)(i,j) ) / (grid->dy);
}


//! \brief Returns the x-derivative at i,j approximated using centered finite
//! differences. Respects grid periodicity and uses one-sided FD at grid edges
//! if necessary.
PetscScalar IceModelVec2S::diff_x_p(int i, int j) {
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
PetscScalar IceModelVec2S::diff_y_p(int i, int j) {
  if (grid->periodicity & Y_PERIODIC)
    return diff_y(i,j);
  
  if (j == 0)
    return ( (*this)(i,j + 1) - (*this)(i,j) ) / (grid->dy);
  else if (j == grid->My - 1)
    return ( (*this)(i,j) - (*this)(i,j - 1) ) / (grid->dy);
  else
    return diff_y(i,j);
}

//! Sums up all the values in an IceModelVec2S object. Ignores ghosts.
/*! Avoids copying to a "global" vector.
 */
PetscErrorCode IceModelVec2S::sum(PetscScalar &result) {
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


// IceModelVec2

PetscErrorCode IceModelVec2::get_component(int N, Vec result) {
  PetscErrorCode ierr;
  void *tmp_res = NULL, *tmp_v;

  if (N < 0 || N >= dof)
    SETERRQ(1, "invalid argument (N)");

  ierr = DAVecGetArray(grid->da2, result, &tmp_res); CHKERRQ(ierr);
  PetscScalar **res = static_cast<PetscScalar**>(tmp_res);

  ierr = DAVecGetArrayDOF(da, v, &tmp_v); CHKERRQ(ierr);
  PetscScalar ***a_dof = static_cast<PetscScalar***>(tmp_v);

  for (PetscInt i = grid->xs; i < grid->xs+grid->xm; ++i)
    for (PetscInt j = grid->ys; j < grid->ys+grid->ym; ++j)
      res[i][j] = a_dof[i][j][N];


  ierr = DAVecRestoreArray(grid->da2, result, &tmp_res); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(da,        v,      &tmp_v); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec2::set_component(int N, Vec source) {
  PetscErrorCode ierr;
  void *tmp_src = NULL, *tmp_v;

  if (N < 0 || N >= dof)
    SETERRQ(1, "invalid argument (N)");

  ierr = DAVecGetArray(grid->da2, source, &tmp_src); CHKERRQ(ierr);
  PetscScalar **src = static_cast<PetscScalar**>(tmp_src);

  ierr = DAVecGetArrayDOF(da, v, &tmp_v); CHKERRQ(ierr);
  PetscScalar ***a_dof = static_cast<PetscScalar***>(tmp_v);

  for (PetscInt i = grid->xs; i < grid->xs+grid->xm; ++i)
    for (PetscInt j = grid->ys; j < grid->ys+grid->ym; ++j)
      a_dof[i][j][N] = src[i][j];

  ierr = DAVecRestoreArray(grid->da2, source, &tmp_src); CHKERRQ(ierr);
  ierr = DAVecRestoreArray(da,        v,      &tmp_v); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode  IceModelVec2::create(IceGrid &my_grid, const char my_name[], bool local,
                                     DAStencilType my_sten, int stencil_width, int my_dof) {
  PetscErrorCode ierr;

  if (!utIsInit()) {
    SETERRQ(1, "PISM ERROR: UDUNITS *was not* initialized.\n");
  }
  if (v != PETSC_NULL) {
    SETERRQ1(2,"IceModelVec2 with name='%s' already allocated\n", my_name);
  }

  dof  = my_dof;
  grid = &my_grid;
  dims = GRID_2D;

  ierr = DACreate2d(my_grid.com, DA_XYPERIODIC, my_sten,
		    grid->My, grid->Mx,
		    grid->Ny, grid->Nx,
		    dof, stencil_width,
                    grid->procs_y, grid->procs_x, &da); CHKERRQ(ierr);

  if (local) {
    ierr = DACreateLocalVector(da, &v); CHKERRQ(ierr);
  } else {
    ierr = DACreateGlobalVector(da, &v); CHKERRQ(ierr);
  }

  localp = local;
  name = my_name;

  vars.resize(dof);
  for (int j = 0; j < dof; ++j)
    vars[j].init(my_name, my_grid, GRID_2D);

  //  ierr = this->set(GSL_NAN); CHKERRQ(ierr);

  return 0;
}

// IceModelVec2Stag

PetscErrorCode  IceModelVec2Stag::create(IceGrid &my_grid, const char my_short_name[], bool local,
					 int stencil_width) {

  PetscErrorCode ierr = IceModelVec2::create(my_grid, my_short_name, local, DA_STENCIL_BOX,
					     stencil_width, dof); CHKERRQ(ierr);
  string s_name = name;
  vars[0].init(s_name + "[0]", my_grid, GRID_2D);
  vars[1].init(s_name + "[1]", my_grid, GRID_2D);

  return 0;
}

PetscErrorCode  IceModelVec2Stag::begin_access() {
  PetscErrorCode ierr;
#ifdef PISM_DEBUG
  ierr = checkAllocated(); CHKERRQ(ierr);

  if (access_counter < 0)
    SETERRQ(1, "IceModelVec::begin_access(): access_counter < 0");
#endif

  if (access_counter == 0) {
    ierr = DAVecGetArrayDOF(da, v, &array); CHKERRQ(ierr);
  }

  access_counter++;

  return 0;
}

//! Checks if an IceModelVec is allocated and calls DAVecRestoreArray.
PetscErrorCode  IceModelVec2Stag::end_access() {
  PetscErrorCode ierr;
  access_counter--;

#ifdef PISM_DEBUG
  ierr = checkAllocated(); CHKERRQ(ierr);

  if (access_counter < 0)
    SETERRQ(1, "IceModelVec::end_access(): access_counter < 0");
#endif

  if (access_counter == 0) {
    ierr = DAVecRestoreArrayDOF(da, v, &array); CHKERRQ(ierr);
    array = NULL;
  }

  return 0;
}

PetscScalar& IceModelVec2Stag::operator() (int i, int j, int k) {
  return static_cast<PetscScalar***>(array)[i][j][k];
}

PetscErrorCode IceModelVec2Stag::get_array(PetscScalar*** &a) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  a = static_cast<PetscScalar***>(array);
  return 0;
}

