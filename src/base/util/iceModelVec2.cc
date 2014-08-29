// Copyright (C) 2008--2014 Ed Bueler and Constantine Khroulev
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

#include <cstring>
#include <cstdlib>
#include <petscdmda.h>
#include <gsl/gsl_math.h>

#include "PIO.hh"
#include "iceModelVec.hh"
#include "IceGrid.hh"
#include "LocalInterpCtx.hh"
#include "iceModelVec_helpers.hh"

#include <assert.h>

namespace pism {

// this file contains methods for derived classes IceModelVec2S and IceModelVec2Int

// methods for base class IceModelVec are in "iceModelVec.cc"

PetscErrorCode  IceModelVec2S::create(IceGrid &my_grid, const std::string &my_name, IceModelVecKind ghostedp, int width) {
  assert(v == NULL);
  PetscErrorCode ierr = IceModelVec2::create(my_grid, my_name, ghostedp, width, m_dof); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IceModelVec2S::get_array(double** &a) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  a = static_cast<double**>(array);
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
  assert(v != NULL);

  PISMDM::Ptr da2;
  ierr = grid->get_dm(1, this->get_stencil_width(), da2); CHKERRQ(ierr);

  if (m_has_ghosts == true) {
    ierr = DMLocalToGlobalBegin(m_da->get(), v, INSERT_VALUES, g2); CHKERRQ(ierr);
    ierr =   DMLocalToGlobalEnd(m_da->get(), v, INSERT_VALUES, g2); CHKERRQ(ierr);
  } else {
    ierr = this->copy_to_vec(g2); CHKERRQ(ierr);
  }

  ierr = DMDAGlobalToNaturalBegin(da2->get(), g2, INSERT_VALUES, g2natural); CHKERRQ(ierr);
  ierr =   DMDAGlobalToNaturalEnd(da2->get(), g2, INSERT_VALUES, g2natural); CHKERRQ(ierr);

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
  assert(v != NULL);

  PISMDM::Ptr da2;
  ierr = grid->get_dm(1, grid->max_stencil_width, da2); CHKERRQ(ierr);

  ierr = VecScatterBegin(ctx, onp0, g2natural, INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);
  ierr =   VecScatterEnd(ctx, onp0, g2natural, INSERT_VALUES, SCATTER_REVERSE); CHKERRQ(ierr);

  ierr = DMDANaturalToGlobalBegin(da2->get(), g2natural, INSERT_VALUES, g2); CHKERRQ(ierr);
  ierr =   DMDANaturalToGlobalEnd(da2->get(), g2natural, INSERT_VALUES, g2); CHKERRQ(ierr);

  if (m_has_ghosts == true) {
    ierr =   DMGlobalToLocalBegin(m_da->get(), g2, INSERT_VALUES, v); CHKERRQ(ierr);
    ierr =     DMGlobalToLocalEnd(m_da->get(), g2, INSERT_VALUES, v); CHKERRQ(ierr);
  } else {
    ierr = this->copy_from_vec(g2); CHKERRQ(ierr);
  }

  return 0;
}

//! Sets an IceModelVec2 to the magnitude of a 2D vector field with components `v_x` and `v_y`.
/*! Computes the magnitude \b pointwise, so any of v_x, v_y and the IceModelVec
  this is called on can be the same.

  Does not communicate.
 */
PetscErrorCode IceModelVec2S::set_to_magnitude(IceModelVec2S &v_x, IceModelVec2S &v_y) {
  PetscErrorCode ierr;
  double **mag = NULL;
  ierr = v_x.begin_access(); CHKERRQ(ierr);
  ierr = v_y.begin_access(); CHKERRQ(ierr);
  ierr = get_array(mag); CHKERRQ(ierr);
  for (int i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (int j=grid->ys; j<grid->ys+grid->ym; ++j) {
      mag[i][j] = sqrt(PetscSqr(v_x(i,j)) + PetscSqr(v_y(i,j)));
    }
  }
  ierr = v_x.end_access(); CHKERRQ(ierr);
  ierr = v_y.end_access(); CHKERRQ(ierr);
  ierr = end_access(); CHKERRQ(ierr);
  
  return 0;
}

//! Masks out all the areas where \f$ M \le 0 \f$ by setting them to `fill`. 
PetscErrorCode IceModelVec2S::mask_by(IceModelVec2S &M, double fill) {
  PetscErrorCode ierr;
  double **a = NULL;
  ierr = get_array(a); CHKERRQ(ierr);
  ierr = M.begin_access(); CHKERRQ(ierr);
  for (int i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (int j=grid->ys; j<grid->ys+grid->ym; ++j) {
      if (M(i,j) <= 0.0)
        a[i][j] = fill;
    }
  }
  ierr = end_access(); CHKERRQ(ierr);
  ierr = M.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec2::write_impl(const PIO &nc, IO_Type nctype) {
  PetscErrorCode ierr;

  assert(v != NULL);

  Vec tmp;                      // a temporary one-component vector,
                                // distributed across processors the same way v is

  // The simplest case:
  if ((m_dof == 1) && (m_has_ghosts == false)) {
    ierr = IceModelVec::write_impl(nc, nctype); CHKERRQ(ierr);
    return 0;
  }

  PISMDM::Ptr da2;
  ierr = grid->get_dm(1, grid->max_stencil_width, da2); CHKERRQ(ierr);

  ierr = DMGetGlobalVector(da2->get(), &tmp); CHKERRQ(ierr);

  if (getVerbosityLevel() > 3) {
    ierr = PetscPrintf(grid->com, "  Writing %s...\n", m_name.c_str()); CHKERRQ(ierr);
  }

  for (unsigned int j = 0; j < m_dof; ++j) {
    ierr = IceModelVec2::get_component(j, tmp); CHKERRQ(ierr);

    ierr = m_metadata[j].write(nc, nctype, write_in_glaciological_units, tmp);
  }

  // Clean up:
  ierr = DMRestoreGlobalVector(da2->get(), &tmp); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IceModelVec2::read_impl(const PIO &nc, const unsigned int time) {
  PetscErrorCode ierr;

  if ((m_dof == 1) && (m_has_ghosts == false)) {
    ierr = IceModelVec::read_impl(nc, time); CHKERRQ(ierr);
    return 0;
  }

  if (getVerbosityLevel() > 3) {
    ierr = PetscPrintf(grid->com, "  Reading %s...\n", m_name.c_str()); CHKERRQ(ierr);
  }

  assert(v != NULL);

  PISMDM::Ptr da2;
  ierr = grid->get_dm(1, grid->max_stencil_width, da2); CHKERRQ(ierr);

  Vec tmp;                      // a temporary one-component vector,
                                // distributed across processors the same way v is
  ierr = DMGetGlobalVector(da2->get(), &tmp); CHKERRQ(ierr);

  for (unsigned int j = 0; j < m_dof; ++j) {
    ierr = m_metadata[j].read(nc, time, tmp); CHKERRQ(ierr);
    ierr = IceModelVec2::set_component(j, tmp); CHKERRQ(ierr);
  }
  
  // The calls above only set the values owned by a processor, so we need to
  // communicate if m_has_ghosts == true:
  if (m_has_ghosts) {
    ierr = update_ghosts(); CHKERRQ(ierr);
  }

  // Clean up:
  ierr = DMRestoreGlobalVector(da2->get(), &tmp); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode IceModelVec2::regrid_impl(const PIO &nc, RegriddingFlag flag,
                                         double default_value) {
  PetscErrorCode ierr;
  LocalInterpCtx *lic = NULL;

  if ((m_dof == 1) && (m_has_ghosts == false)) {
    ierr = IceModelVec::regrid_impl(nc, flag, default_value); CHKERRQ(ierr);
    return 0;
  }

  if (getVerbosityLevel() > 3) {
    ierr = PetscPrintf(grid->com, "  Regridding %s...\n", m_name.c_str()); CHKERRQ(ierr);
  }

  PISMDM::Ptr da2;
  ierr = grid->get_dm(1, grid->max_stencil_width, da2); CHKERRQ(ierr);

  Vec tmp;                      // a temporary one-component vector,
                                // distributed across processors the same way v is
  ierr = DMGetGlobalVector(da2->get(), &tmp); CHKERRQ(ierr);

  for (unsigned int j = 0; j < m_dof; ++j) {
    ierr = m_metadata[j].regrid(nc, flag, m_report_range, default_value, tmp); CHKERRQ(ierr);
    ierr = IceModelVec2::set_component(j, tmp); CHKERRQ(ierr);
  }

  // The calls above only set the values owned by a processor, so we need to
  // communicate if m_has_ghosts == true:
  if (m_has_ghosts == true) {
    ierr = update_ghosts(); CHKERRQ(ierr);
  }

  // Clean up:
  ierr = DMRestoreGlobalVector(da2->get(), &tmp); CHKERRQ(ierr);
  delete lic;
  return 0;
}

//! \brief View a 2D field.
PetscErrorCode IceModelVec2::view(int viewer_size) {
  PetscErrorCode ierr;
  PetscViewer viewers[2] = {NULL, NULL};

  if (m_dof > 2) SETERRQ(grid->com, 1, "dof > 2 is not supported");

  for (unsigned int j = 0; j < m_dof; ++j) {
    std::string c_name = m_metadata[j].get_name(),
      long_name = m_metadata[j].get_string("long_name"),
      units = m_metadata[j].get_string("glaciological_units"),
      title = long_name + " (" + units + ")";

    if (map_viewers[c_name] == NULL) {
      ierr = grid->create_viewer(viewer_size, title, map_viewers[c_name]); CHKERRQ(ierr);
    }

    viewers[j] = map_viewers[c_name];
  }

  ierr = view(viewers[0], viewers[1]); CHKERRQ(ierr); 

  return 0;
}

//! \brief View a 2D vector field using existing PETSc viewers.
//! Allocates and de-allocates g2, the temporary global vector; performance
//! should not matter here.
PetscErrorCode IceModelVec2::view(PetscViewer v1, PetscViewer v2) {
  PetscErrorCode ierr;
  Vec tmp;

  PISMDM::Ptr da2;
  ierr = grid->get_dm(1, grid->max_stencil_width, da2); CHKERRQ(ierr);

  ierr = DMGetGlobalVector(da2->get(), &tmp); CHKERRQ(ierr);

  PetscViewer viewers[2] = {v1, v2};

  for (unsigned int i = 0; i < m_dof; ++i) {
    std::string long_name = m_metadata[i].get_string("long_name"),
      units = m_metadata[i].get_string("glaciological_units"),
      title = long_name + " (" + units + ")";

    PetscDraw draw;
    ierr = PetscViewerDrawGetDraw(viewers[i], 0, &draw); CHKERRQ(ierr);
    ierr = PetscDrawSetTitle(draw, title.c_str()); CHKERRQ(ierr);

    ierr = IceModelVec2::get_component(i, tmp); CHKERRQ(ierr);

    ierr = convert_vec(tmp,
                       m_metadata[i].get_units(),
                       m_metadata[i].get_glaciological_units()); CHKERRQ(ierr);

    ierr = VecView(tmp, viewers[i]); CHKERRQ(ierr);
  }

  ierr = DMRestoreGlobalVector(da2->get(), &tmp); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec2S::view_matlab(PetscViewer my_viewer) {
  PetscErrorCode ierr;
  std::string long_name = metadata().get_string("long_name");
  Vec tmp;

  PISMDM::Ptr da2;
  ierr = grid->get_dm(1, grid->max_stencil_width, da2); CHKERRQ(ierr);

  ierr = DMGetGlobalVector(da2->get(), &tmp); CHKERRQ(ierr);

  if (m_has_ghosts) {
    ierr = copy_to_vec(tmp); CHKERRQ(ierr);
  } else {
    ierr = VecCopy(v, tmp); CHKERRQ(ierr);
  }

  ierr = convert_vec(tmp,
                     metadata().get_units(),
                     metadata().get_glaciological_units()); CHKERRQ(ierr);

  // add Matlab comment before listing, using short title

  ierr = PetscViewerASCIIPrintf(my_viewer, "\n%%%% %s = %s\n",
                                m_name.c_str(), long_name.c_str()); CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) tmp, m_name.c_str()); CHKERRQ(ierr);

  ierr = VecView(tmp, my_viewer); CHKERRQ(ierr);

  ierr = PetscViewerASCIIPrintf(my_viewer,"\n%s = reshape(%s,%d,%d);\n\n",
                                m_name.c_str(), m_name.c_str(), grid->My, grid->Mx); CHKERRQ(ierr);

  ierr = DMRestoreGlobalVector(da2->get(), &tmp); CHKERRQ(ierr);

  return 0;
}

//! Checks if the current IceModelVec2S has NANs and reports if it does.
/*! Up to a fixed number of messages are printed at stdout.  Returns the full
 count of NANs (which is a nonzero) on this rank. */
PetscErrorCode IceModelVec2S::has_nan() {
  PetscErrorCode ierr;
  const int max_print_this_rank=10;
  int retval=0;

  ierr = begin_access(); CHKERRQ(ierr);
  int i, j;
  for (i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (j=grid->ys; j<grid->ys+grid->ym; ++j) {
      if (gsl_isnan((*this)(i,j))) {
        retval++;
        if (retval <= max_print_this_rank) {
          ierr = PetscSynchronizedPrintf(grid->com, 
             "IceModelVec2S %s: NAN (or uninitialized) at i = %d, j = %d on rank = %d\n",
             m_name.c_str(), i, j, grid->rank); CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = end_access(); CHKERRQ(ierr);

  if (retval > 0) {
    ierr = PetscSynchronizedPrintf(grid->com, 
       "IceModelVec2S %s: detected %d NANs (or uninitialized) on rank = %d\n",
             m_name.c_str(), retval, grid->rank); CHKERRQ(ierr);
  }

  ierr = PetscSynchronizedFlush(grid->com); CHKERRQ(ierr);
  return retval;
}

//! \brief Returns the x-derivative at i,j approximated using centered finite
//! differences.
double IceModelVec2S::diff_x(int i, int j) {
  return ((*this)(i + 1,j) - (*this)(i - 1,j)) / (2 * grid->dx);
}

//! \brief Returns the y-derivative at i,j approximated using centered finite
//! differences.
double IceModelVec2S::diff_y(int i, int j) {
  return ((*this)(i,j + 1) - (*this)(i,j - 1)) / (2 * grid->dy);
}


//! \brief Returns the x-derivative at East staggered point i+1/2,j approximated 
//! using centered (obvious) finite differences.
double IceModelVec2S::diff_x_stagE(int i, int j) {
  return ((*this)(i+1,j) - (*this)(i,j)) / (grid->dx);
}

//! \brief Returns the y-derivative at East staggered point i+1/2,j approximated 
//! using centered finite differences.
double IceModelVec2S::diff_y_stagE(int i, int j) {
  return ((*this)(i+1,j+1) + (*this)(i,j+1)
           - (*this)(i+1,j-1) - (*this)(i,j-1)) / (4* grid->dy);
}

//! \brief Returns the x-derivative at North staggered point i,j+1/2 approximated 
//! using centered finite differences.
double IceModelVec2S::diff_x_stagN(int i, int j) {
  return ((*this)(i+1,j+1) + (*this)(i+1,j)
           - (*this)(i-1,j+1) - (*this)(i-1,j)) / (4* grid->dx);
}

//! \brief Returns the y-derivative at North staggered point i,j+1/2 approximated 
//! using centered (obvious) finite differences.
double IceModelVec2S::diff_y_stagN(int i, int j) {
  return ((*this)(i,j+1) - (*this)(i,j)) / (grid->dy);
}


//! \brief Returns the x-derivative at i,j approximated using centered finite
//! differences. Respects grid periodicity and uses one-sided FD at grid edges
//! if necessary.
double IceModelVec2S::diff_x_p(int i, int j) {
  if (grid->periodicity & X_PERIODIC)
    return diff_x(i,j);
  
  if (i == 0)
    return ((*this)(i + 1,j) - (*this)(i,j)) / (grid->dx);
  else if (i == grid->Mx - 1)
    return ((*this)(i,j) - (*this)(i - 1,j)) / (grid->dx);
  else
    return diff_x(i,j);
}

//! \brief Returns the y-derivative at i,j approximated using centered finite
//! differences. Respects grid periodicity and uses one-sided FD at grid edges
//! if necessary.
double IceModelVec2S::diff_y_p(int i, int j) {
  if (grid->periodicity & Y_PERIODIC)
    return diff_y(i,j);
  
  if (j == 0)
    return ((*this)(i,j + 1) - (*this)(i,j)) / (grid->dy);
  else if (j == grid->My - 1)
    return ((*this)(i,j) - (*this)(i,j - 1)) / (grid->dy);
  else
    return diff_y(i,j);
}

//! Sums up all the values in an IceModelVec2S object. Ignores ghosts.
/*! Avoids copying to a "global" vector.
 */
PetscErrorCode IceModelVec2S::sum(double &result) {
  PetscErrorCode ierr;
  double my_result = 0;

  // sum up the local part:
  ierr = begin_access(); CHKERRQ(ierr);
  for (int i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (int j=grid->ys; j<grid->ys+grid->ym; ++j) {
      my_result += (*this)(i,j);
    }
  }
  ierr = end_access(); CHKERRQ(ierr);

  // find the global sum:
  ierr = GlobalSum(&my_result, &result, grid->com); CHKERRQ(ierr);

  return 0;
}


//! Finds maximum over all the values in an IceModelVec2S object.  Ignores ghosts.
PetscErrorCode IceModelVec2S::max(double &result) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  double my_result = (*this)(grid->xs,grid->ys);
  for (int i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (int j=grid->ys; j<grid->ys+grid->ym; ++j) {
      my_result = PetscMax(my_result,(*this)(i,j));
    }
  }
  ierr = end_access(); CHKERRQ(ierr);
  ierr = GlobalMax(&my_result, &result, grid->com); CHKERRQ(ierr);
  return 0;
}


//! Finds maximum over all the absolute values in an IceModelVec2S object.  Ignores ghosts.
PetscErrorCode IceModelVec2S::absmax(double &result) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  double my_result = 0.0;
  for (int i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (int j=grid->ys; j<grid->ys+grid->ym; ++j) {
      my_result = PetscMax(my_result,PetscAbs((*this)(i,j)));
    }
  }
  ierr = end_access(); CHKERRQ(ierr);
  ierr = GlobalMax(&my_result, &result, grid->com); CHKERRQ(ierr);
  return 0;
}


//! Finds minimum over all the values in an IceModelVec2S object.  Ignores ghosts.
PetscErrorCode IceModelVec2S::min(double &result) {
  PetscErrorCode ierr;
  ierr = begin_access(); CHKERRQ(ierr);
  double my_result = (*this)(grid->xs,grid->ys);
  for (int i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (int j=grid->ys; j<grid->ys+grid->ym; ++j) {
      my_result = PetscMin(my_result,(*this)(i,j));
    }
  }
  ierr = end_access(); CHKERRQ(ierr);
  ierr = GlobalMin(&my_result, &result, grid->com); CHKERRQ(ierr);
  return 0;
}


// IceModelVec2

/*!
 * This could be implemented using VecStrideGather, but our code is more
 * flexible: `source` and the current IceModelVec2 need not be both local or
 * global.
 */

PetscErrorCode IceModelVec2::get_component(unsigned int N, Vec result) {
  PetscErrorCode ierr;
  void *tmp_res = NULL, *tmp_v;

  if (N >= m_dof)
    SETERRQ(grid->com, 1, "invalid argument (N)");

  PISMDM::Ptr da2;
  ierr = grid->get_dm(1, grid->max_stencil_width, da2); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(da2->get(), result, &tmp_res); CHKERRQ(ierr);
  double **res = static_cast<double**>(tmp_res);

  ierr = DMDAVecGetArrayDOF(m_da->get(), v, &tmp_v); CHKERRQ(ierr);
  double ***a_dof = static_cast<double***>(tmp_v);

  for (int i = grid->xs; i < grid->xs+grid->xm; ++i)
    for (int j = grid->ys; j < grid->ys+grid->ym; ++j)
      res[i][j] = a_dof[i][j][N];

  ierr = DMDAVecRestoreArray(da2->get(), result, &tmp_res); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(m_da->get(), v, &tmp_v); CHKERRQ(ierr);

  return 0;
}

/*!
 * This could be implemented using VecStrideScatter, but our code is more
 * flexible: `source` and the current IceModelVec2 need not be both local or
 * global.
 */
PetscErrorCode IceModelVec2::set_component(unsigned int N, Vec source) {
  PetscErrorCode ierr;
  void *tmp_src = NULL, *tmp_v;

  if (N >= m_dof)
    SETERRQ(grid->com, 1, "invalid argument (N)");

  PISMDM::Ptr da2;
  ierr = grid->get_dm(1, grid->max_stencil_width, da2); CHKERRQ(ierr);

  ierr = DMDAVecGetArray(da2->get(), source, &tmp_src); CHKERRQ(ierr);
  double **src = static_cast<double**>(tmp_src);

  ierr = DMDAVecGetArrayDOF(m_da->get(), v, &tmp_v); CHKERRQ(ierr);
  double ***a_dof = static_cast<double***>(tmp_v);

  for (int i = grid->xs; i < grid->xs+grid->xm; ++i)
    for (int j = grid->ys; j < grid->ys+grid->ym; ++j)
      a_dof[i][j][N] = src[i][j];

  ierr = DMDAVecRestoreArray(da2->get(), source, &tmp_src); CHKERRQ(ierr);
  ierr = DMDAVecRestoreArrayDOF(m_da->get(), v, &tmp_v); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec2::get_component(unsigned int n, IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = IceModelVec2::get_component(n, result.v); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModelVec2::set_component(unsigned int n, IceModelVec2S &source) {
  PetscErrorCode ierr;

  ierr = IceModelVec2::set_component(n, source.v); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode  IceModelVec2::create(IceGrid &my_grid, const std::string & my_name, IceModelVecKind ghostedp,
                                     unsigned int stencil_width, int my_dof) {
  PetscErrorCode ierr;

  assert(v == NULL);

  m_dof  = my_dof;
  grid = &my_grid;

  if ((m_dof != 1) || (stencil_width > grid->max_stencil_width)) {
    m_da_stencil_width = stencil_width;
  } else {
    m_da_stencil_width = grid->max_stencil_width;
  }

  // initialize the da member:
  ierr = grid->get_dm(this->m_dof, this->m_da_stencil_width, m_da); CHKERRQ(ierr);

  if (ghostedp) {
    ierr = DMCreateLocalVector(m_da->get(), &v); CHKERRQ(ierr);
  } else {
    ierr = DMCreateGlobalVector(m_da->get(), &v); CHKERRQ(ierr);
  }

  m_has_ghosts = (ghostedp == WITH_GHOSTS);
  m_name       = my_name;

  m_metadata.resize(m_dof, NCSpatialVariable(grid->get_unit_system()));

  if (m_dof == 1) {
    metadata().init_2d(my_name, my_grid);
  } else {

    for (unsigned int j = 0; j < m_dof; ++j) {
      char tmp[TEMPORARY_STRING_LENGTH];

      snprintf(tmp, TEMPORARY_STRING_LENGTH, "%s[%d]",
               m_name.c_str(), j);
      metadata(j).init_2d(tmp, my_grid);
    }
  }

  return 0;
}

PetscErrorCode IceModelVec2S::add(double alpha, IceModelVec &x) {
  return add_2d<IceModelVec2S>(this, alpha, &x, this);
}

PetscErrorCode IceModelVec2S::add(double alpha, IceModelVec &x, IceModelVec &result) {
  return add_2d<IceModelVec2S>(this, alpha, &x, &result);
}

PetscErrorCode IceModelVec2S::copy_to(IceModelVec &destination) {
  return copy_2d<IceModelVec2S>(this, &destination);
}

// IceModelVec2Stag
PetscErrorCode IceModelVec2Stag::create(IceGrid &my_grid, const std::string &my_short_name, IceModelVecKind ghostedp,
                                        unsigned int stencil_width) {
  PetscErrorCode ierr;

  ierr = IceModelVec2::create(my_grid, my_short_name, ghostedp, stencil_width, m_dof); CHKERRQ(ierr);

  return 0;
}

//! Averages staggered grid values of a scalar field and puts them on a regular grid.
/*!
 * The current IceModelVec needs to have ghosts.
 */
PetscErrorCode IceModelVec2Stag::staggered_to_regular(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = begin_access(); CHKERRQ(ierr);
  for (int i = grid->xs; i < grid->xs+grid->xm; ++i) {
    for (int j = grid->ys; j < grid->ys+grid->ym; ++j) {
      result(i,j) = 0.25 * ((*this)(i,j,0) + (*this)(i,j,1)
                              + (*this)(i,j-1,1) + (*this)(i-1,j,0));
    } // j
  }   // i
  ierr = end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

//! \brief Averages staggered grid values of a 2D vector field (u on the
//! i-offset, v on the j-offset) and puts them on a regular grid.
/*!
 * The current IceModelVec needs to have ghosts.
 */
PetscErrorCode IceModelVec2Stag::staggered_to_regular(IceModelVec2V &result) {
  PetscErrorCode ierr;

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = begin_access(); CHKERRQ(ierr);
  for (int i = grid->xs; i < grid->xs+grid->xm; ++i) {
    for (int j = grid->ys; j < grid->ys+grid->ym; ++j) {
        result(i,j).u = 0.5 * ((*this)(i-1,j,0) + (*this)(i,j,0));
        result(i,j).v = 0.5 * ((*this)(i,j-1,1) + (*this)(i,j,1));
    } // j
  }   // i
  ierr = end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}


//! For each component, finds the maximum over all the absolute values.  Ignores ghosts.
/*!
Assumes z is allocated.
 */
PetscErrorCode IceModelVec2Stag::absmaxcomponents(double* z) {
  PetscErrorCode ierr;
  double my_z[2] = {0.0, 0.0};
  ierr = begin_access(); CHKERRQ(ierr);
  for (int i=grid->xs; i<grid->xs+grid->xm; ++i) {
    for (int j=grid->ys; j<grid->ys+grid->ym; ++j) {
      my_z[0] = PetscMax(my_z[0],PetscAbs((*this)(i,j,0)));
      my_z[1] = PetscMax(my_z[1],PetscAbs((*this)(i,j,1)));
    }
  }
  ierr = end_access(); CHKERRQ(ierr);
  ierr = GlobalMax(&(my_z[0]), &(z[0]), grid->com); CHKERRQ(ierr);
  ierr = GlobalMax(&(my_z[1]), &(z[1]), grid->com); CHKERRQ(ierr);
  return 0;
}

} // end of namespace pism
