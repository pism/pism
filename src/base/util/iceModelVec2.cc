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

#include <cassert>

#include "PIO.hh"
#include "iceModelVec.hh"
#include "IceGrid.hh"
#include "LocalInterpCtx.hh"
#include "PISMConfig.hh"

#include "error_handling.hh"
#include "iceModelVec_helpers.hh"

namespace pism {

// this file contains methods for derived classes IceModelVec2S and IceModelVec2Int

// methods for base class IceModelVec are in "iceModelVec.cc"

void  IceModelVec2S::create(IceGrid &my_grid, const std::string &my_name, IceModelVecKind ghostedp, int width) {
  assert(m_v == NULL);
  IceModelVec2::create(my_grid, my_name, ghostedp, width, m_dof);
}

void IceModelVec2S::get_array(double** &a) {
  begin_access();
  a = static_cast<double**>(array);
}

/*! Allocate a copy on processor zero and the scatter needed to move data.
 *
 * The caller is responsible for de-allocating result by calling VecDestroy.
 */
void IceModelVec2S::allocate_proc0_copy(Vec &result) const {
  PetscErrorCode ierr;
  Vec v_proc0 = NULL;

  ierr = PetscObjectQuery((PetscObject)m_da->get(), "v_proc0", (PetscObject*)&v_proc0);
  PISM_PETSC_CHK(ierr, "PetscObjectQuery")
                                                                                          ;
  if (v_proc0 == NULL) {

    Vec natural_work = NULL;
    // create a work vector with natural ordering:
    DMDACreateNaturalVector(*m_da, &natural_work);
    // this increments the reference counter of natural_work
    ierr = PetscObjectCompose((PetscObject)m_da->get(), "natural_work",
                              (PetscObject)natural_work);
    PISM_PETSC_CHK(ierr, "PetscObjectCompose");

    // initialize the scatter to processor 0 and create storage on processor 0
    VecScatter scatter_to_zero = NULL;
    ierr = VecScatterCreateToZero(natural_work, &scatter_to_zero, &v_proc0);
    PISM_PETSC_CHK(ierr, "VecScatterCreateToZero");

    // decrement the reference counter; will be destroyed once m_da is destroyed
    ierr = VecDestroy(&natural_work);
    PISM_PETSC_CHK(ierr, "VecDestroy");

    // this increments the reference counter of scatter_to_zero
    ierr = PetscObjectCompose((PetscObject)m_da->get(), "scatter_to_zero",
                              (PetscObject)scatter_to_zero);
    PISM_PETSC_CHK(ierr, "PetscObjectCompose");
    // decrement the reference counter; will be destroyed once m_da is destroyed
    ierr = VecScatterDestroy(&scatter_to_zero);
    PISM_PETSC_CHK(ierr, "VecScatterDestroy");

    // this increments the reference counter of v_proc0
    ierr = PetscObjectCompose((PetscObject)m_da->get(), "v_proc0",
                              (PetscObject)v_proc0);
    PISM_PETSC_CHK(ierr, "PetscObjectCompose");

    // We DO NOT call VecDestroy(v_proc0): the caller is expected to call VecDestroy
    result = v_proc0;
  } else {
    ierr = VecDuplicate(v_proc0, &result);
    PISM_PETSC_CHK(ierr, "VecDuplicate");
  }
}

//! Puts a local IceModelVec2S on processor 0.
void IceModelVec2S::put_on_proc0(Vec onp0) const {
  PetscErrorCode ierr;
  assert(m_v != NULL);

  VecScatter scatter_to_zero = NULL;
  Vec natural_work = NULL;
  ierr = PetscObjectQuery((PetscObject)m_da->get(), "scatter_to_zero",
                          (PetscObject*)&scatter_to_zero);
  PISM_PETSC_CHK(ierr, "PetscObjectQuery");
  ierr = PetscObjectQuery((PetscObject)m_da->get(), "natural_work",
                          (PetscObject*)&natural_work);
  PISM_PETSC_CHK(ierr, "PetscObjectQuery");

  if (natural_work == NULL || scatter_to_zero == NULL) {
    throw RuntimeError("call allocate_proc0_copy() before calling put_on_proc0");
  }

  Vec global = NULL;

  if (m_has_ghosts) {
    DMGetGlobalVector(*m_da, &global);
    this->copy_to_vec(m_da, global);
  } else {
    global = m_v;
  }

  DMDAGlobalToNaturalBegin(*m_da, global, INSERT_VALUES, natural_work);
  DMDAGlobalToNaturalEnd(*m_da, global, INSERT_VALUES, natural_work);

  if (m_has_ghosts) {
    DMRestoreGlobalVector(*m_da, &global);
  }

  ierr = VecScatterBegin(scatter_to_zero, natural_work, onp0,
                         INSERT_VALUES, SCATTER_FORWARD);
  PISM_PETSC_CHK(ierr, "VecScatterBegin");
  VecScatterEnd(scatter_to_zero, natural_work, onp0,
                INSERT_VALUES, SCATTER_FORWARD);
}

//! Gets a local IceModelVec2 from processor 0.
void IceModelVec2S::get_from_proc0(Vec onp0) {
  PetscErrorCode ierr;
  assert(m_v != NULL);

  VecScatter scatter_to_zero = NULL;
  Vec natural_work = NULL;
  ierr = PetscObjectQuery((PetscObject)m_da->get(), "scatter_to_zero",
                          (PetscObject*)&scatter_to_zero);
  PISM_PETSC_CHK(ierr, "PetscObjectQuery");
  ierr = PetscObjectQuery((PetscObject)m_da->get(), "natural_work",
                          (PetscObject*)&natural_work);
  PISM_PETSC_CHK(ierr, "PetscObjectQuery");

  if (natural_work == NULL || scatter_to_zero == NULL) {
    throw RuntimeError("call allocate_proc0_copy() before calling get_from_proc0");
  }

  ierr = VecScatterBegin(scatter_to_zero, onp0, natural_work,
                         INSERT_VALUES, SCATTER_REVERSE);
  PISM_PETSC_CHK(ierr, "VecScatterBegin");
  VecScatterEnd(scatter_to_zero, onp0, natural_work,
                INSERT_VALUES, SCATTER_REVERSE);

  Vec global = NULL;

  if (m_has_ghosts) {
    DMGetGlobalVector(*m_da, &global);
  } else {
    global = m_v;
  }

  DMDANaturalToGlobalBegin(*m_da, natural_work, INSERT_VALUES, global);
  DMDANaturalToGlobalEnd(*m_da, natural_work, INSERT_VALUES, global);

  if (m_has_ghosts) {
    this->copy_from_vec(global);
    DMRestoreGlobalVector(*m_da, &global);
  }

  inc_state_counter();          // mark as modified
}

//! Sets an IceModelVec2 to the magnitude of a 2D vector field with components `v_x` and `v_y`.
/*! Computes the magnitude \b pointwise, so any of v_x, v_y and the IceModelVec
  this is called on can be the same.

  Does not communicate.
 */
void IceModelVec2S::set_to_magnitude(IceModelVec2S &v_x, IceModelVec2S &v_y) {
  IceModelVec::AccessList list(*this);
  list.add(v_x);
  list.add(v_y);

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*this)(i,j) = sqrt(PetscSqr(v_x(i,j)) + PetscSqr(v_y(i,j)));
  }

  inc_state_counter();          // mark as modified
  
}

//! Masks out all the areas where \f$ M \le 0 \f$ by setting them to `fill`. 
void IceModelVec2S::mask_by(IceModelVec2S &M, double fill) {
  IceModelVec::AccessList list(*this);
  list.add(M);

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (M(i,j) <= 0.0) {
      (*this)(i,j) = fill;
    }
  }

  inc_state_counter();          // mark as modified
}

void IceModelVec2::write_impl(const PIO &nc, IO_Type nctype) const {
  PetscErrorCode ierr;

  assert(m_v != NULL);

  Vec tmp;                      // a temporary one-component vector,
                                // distributed across processors the same way v is

  // The simplest case:
  if ((m_dof == 1) && (m_has_ghosts == false)) {
    IceModelVec::write_impl(nc, nctype);
    return;
  }

  // Get the dof=1, stencil_width=0 DMDA (components are always scalar
  // and we just need a global Vec):
  PISMDM::Ptr da2 = grid->get_dm(1, 0);

  DMGetGlobalVector(*da2, &tmp);

  if (getVerbosityLevel() > 3) {
    ierr = PetscPrintf(grid->com, "  Writing %s...\n", m_name.c_str());
    PISM_PETSC_CHK(ierr, "PetscPrintf");
  }

  for (unsigned int j = 0; j < m_dof; ++j) {
    IceModelVec2::get_dof(da2, tmp, j);

    ierr = m_metadata[j].write(nc, nctype, write_in_glaciological_units, tmp);
  }

  // Clean up:
  DMRestoreGlobalVector(*da2, &tmp);
}

void IceModelVec2::read_impl(const PIO &nc, const unsigned int time) {
  PetscErrorCode ierr;

  if ((m_dof == 1) && (m_has_ghosts == false)) {
    IceModelVec::read_impl(nc, time);
    return;
  }

  if (getVerbosityLevel() > 3) {
    ierr = PetscPrintf(grid->com, "  Reading %s...\n", m_name.c_str());
    PISM_PETSC_CHK(ierr, "PetscPrintf");
  }

  assert(m_v != NULL);

  // Get the dof=1, stencil_width=0 DMDA (components are always scalar
  // and we just need a global Vec):
  PISMDM::Ptr da2 = grid->get_dm(1, 0);

  Vec tmp;                      // a temporary one-component vector,
                                // distributed across processors the same way v is
  DMGetGlobalVector(*da2, &tmp);

  for (unsigned int j = 0; j < m_dof; ++j) {
    m_metadata[j].read(nc, time, tmp);
    IceModelVec2::set_dof(da2, tmp, j);
  }
  
  // The calls above only set the values owned by a processor, so we need to
  // communicate if m_has_ghosts == true:
  if (m_has_ghosts) {
    update_ghosts();
  }

  // Clean up:
  DMRestoreGlobalVector(*da2, &tmp);
}

void IceModelVec2::regrid_impl(const PIO &nc, RegriddingFlag flag,
                                         double default_value) {
  PetscErrorCode ierr;

  if ((m_dof == 1) && (m_has_ghosts == false)) {
    IceModelVec::regrid_impl(nc, flag, default_value);
    return;
  }

  if (getVerbosityLevel() > 3) {
    ierr = PetscPrintf(grid->com, "  Regridding %s...\n", m_name.c_str());
    PISM_PETSC_CHK(ierr, "PetscPrintf");
  }

  // Get the dof=1, stencil_width=0 DMDA (components are always scalar
  // and we just need a global Vec):
  PISMDM::Ptr da2 = grid->get_dm(1, 0);

  Vec tmp;                      // a temporary one-component vector,
                                // distributed across processors the same way v is
  DMGetGlobalVector(*da2, &tmp);

  for (unsigned int j = 0; j < m_dof; ++j) {
    m_metadata[j].regrid(nc, flag, m_report_range, default_value, tmp);
    IceModelVec2::set_dof(da2, tmp, j);
  }

  // The calls above only set the values owned by a processor, so we need to
  // communicate if m_has_ghosts == true:
  if (m_has_ghosts == true) {
    update_ghosts();
  }

  // Clean up:
  DMRestoreGlobalVector(*da2, &tmp);
}

//! \brief View a 2D field.
void IceModelVec2::view(int viewer_size) const {
  PetscViewer viewers[2] = {NULL, NULL};

  if (m_dof > 2) {
    throw RuntimeError("dof > 2 is not supported");
  }

  for (unsigned int j = 0; j < m_dof; ++j) {
    std::string c_name = m_metadata[j].get_name(),
      long_name = m_metadata[j].get_string("long_name"),
      units = m_metadata[j].get_string("glaciological_units"),
      title = long_name + " (" + units + ")";

    if (map_viewers[c_name] == NULL) {
      grid->create_viewer(viewer_size, title, map_viewers[c_name]);
    }

    viewers[j] = map_viewers[c_name];
  }

  view(viewers[0], viewers[1]); 
}

//! \brief View a 2D vector field using existing PETSc viewers.
//! Allocates and de-allocates g2, the temporary global vector; performance
//! should not matter here.
void IceModelVec2::view(PetscViewer v1, PetscViewer v2) const {
  PetscErrorCode ierr;
  Vec tmp;

  // Get the dof=1, stencil_width=0 DMDA (components are always scalar
  // and we just need a global Vec):
  PISMDM::Ptr da2 = grid->get_dm(1, 0);

  DMGetGlobalVector(*da2, &tmp);

  PetscViewer viewers[2] = {v1, v2};

  for (unsigned int i = 0; i < m_dof; ++i) {
    std::string long_name = m_metadata[i].get_string("long_name"),
      units = m_metadata[i].get_string("glaciological_units"),
      title = long_name + " (" + units + ")";

    PetscDraw draw;
    ierr = PetscViewerDrawGetDraw(viewers[i], 0, &draw);
    PISM_PETSC_CHK(ierr, "PetscViewerDrawGetDraw");
    ierr = PetscDrawSetTitle(draw, title.c_str());
    PISM_PETSC_CHK(ierr, "PetscDrawSetTitle");

    IceModelVec2::get_dof(da2, tmp, i);

    convert_vec(tmp,
                m_metadata[i].get_units(),
                m_metadata[i].get_glaciological_units());

    ierr = VecView(tmp, viewers[i]);
    PISM_PETSC_CHK(ierr, "VecView");
  }

  DMRestoreGlobalVector(*da2, &tmp);
}

//! \brief Returns the x-derivative at i,j approximated using centered finite
//! differences.
double IceModelVec2S::diff_x(int i, int j) const {
  return ((*this)(i + 1,j) - (*this)(i - 1,j)) / (2 * grid->dx);
}

//! \brief Returns the y-derivative at i,j approximated using centered finite
//! differences.
double IceModelVec2S::diff_y(int i, int j) const {
  return ((*this)(i,j + 1) - (*this)(i,j - 1)) / (2 * grid->dy);
}


//! \brief Returns the x-derivative at East staggered point i+1/2,j approximated 
//! using centered (obvious) finite differences.
double IceModelVec2S::diff_x_stagE(int i, int j) const {
  return ((*this)(i+1,j) - (*this)(i,j)) / (grid->dx);
}

//! \brief Returns the y-derivative at East staggered point i+1/2,j approximated 
//! using centered finite differences.
double IceModelVec2S::diff_y_stagE(int i, int j) const {
  return ((*this)(i+1,j+1) + (*this)(i,j+1)
           - (*this)(i+1,j-1) - (*this)(i,j-1)) / (4* grid->dy);
}

//! \brief Returns the x-derivative at North staggered point i,j+1/2 approximated 
//! using centered finite differences.
double IceModelVec2S::diff_x_stagN(int i, int j) const {
  return ((*this)(i+1,j+1) + (*this)(i+1,j)
           - (*this)(i-1,j+1) - (*this)(i-1,j)) / (4* grid->dx);
}

//! \brief Returns the y-derivative at North staggered point i,j+1/2 approximated 
//! using centered (obvious) finite differences.
double IceModelVec2S::diff_y_stagN(int i, int j) const {
  return ((*this)(i,j+1) - (*this)(i,j)) / (grid->dy);
}


//! \brief Returns the x-derivative at i,j approximated using centered finite
//! differences. Respects grid periodicity and uses one-sided FD at grid edges
//! if necessary.
double IceModelVec2S::diff_x_p(int i, int j) const {
  if (grid->periodicity & X_PERIODIC) {
    return diff_x(i,j);
  }
  
  if (i == 0) {
    return ((*this)(i + 1,j) - (*this)(i,j)) / (grid->dx);
  } else if (i == (int)grid->Mx() - 1) {
    return ((*this)(i,j) - (*this)(i - 1,j)) / (grid->dx);
  } else {
    return diff_x(i,j);
 }
}

//! \brief Returns the y-derivative at i,j approximated using centered finite
//! differences. Respects grid periodicity and uses one-sided FD at grid edges
//! if necessary.
double IceModelVec2S::diff_y_p(int i, int j) const {
  if (grid->periodicity & Y_PERIODIC) {
    return diff_y(i,j);
  }
  
  if (j == 0) {
    return ((*this)(i,j + 1) - (*this)(i,j)) / (grid->dy);
  } else if (j == (int)grid->My - 1) {
    return ((*this)(i,j) - (*this)(i,j - 1)) / (grid->dy);
  } else {
    return diff_y(i,j);
  }
}

//! Sums up all the values in an IceModelVec2S object. Ignores ghosts.
/*! Avoids copying to a "global" vector.
 */
void IceModelVec2S::sum(double &result) {
  double my_result = 0;

  IceModelVec::AccessList list(*this);
  
  // sum up the local part:
  for (Points p(*grid); p; p.next()) {
    my_result += (*this)(p.i(), p.j());
  }

  // find the global sum:
  GlobalSum(grid->com, &my_result,  &result);
}


//! Finds maximum over all the values in an IceModelVec2S object.  Ignores ghosts.
void IceModelVec2S::max(double &result) const {
  IceModelVec::AccessList list(*this);

  double my_result = (*this)(grid->xs(),grid->ys());
  for (Points p(*grid); p; p.next()) {
    my_result = PetscMax(my_result,(*this)(p.i(), p.j()));
  }

  GlobalMax(grid->com, &my_result, &result);
}


//! Finds maximum over all the absolute values in an IceModelVec2S object.  Ignores ghosts.
void IceModelVec2S::absmax(double &result) const {

  IceModelVec::AccessList list(*this);
  double my_result = 0.0;
  for (Points p(*grid); p; p.next()) {
    my_result = PetscMax(my_result,PetscAbs((*this)(p.i(), p.j())));
  }

  GlobalMax(grid->com, &my_result,  &result);
}


//! Finds minimum over all the values in an IceModelVec2S object.  Ignores ghosts.
void IceModelVec2S::min(double &result) const {
  IceModelVec::AccessList list(*this);

  double my_result = (*this)(grid->xs(),grid->ys());
  for (Points p(*grid); p; p.next()) {
    my_result = PetscMin(my_result,(*this)(p.i(), p.j()));
  }

  GlobalMin(grid->com, &my_result,  &result);
}


// IceModelVec2

void IceModelVec2::get_component(unsigned int n, IceModelVec2S &result) const {

  IceModelVec2::get_dof(result.get_dm(), result.m_v, n);
}

void IceModelVec2::set_component(unsigned int n, IceModelVec2S &source) {

  IceModelVec2::set_dof(source.get_dm(), source.m_v, n);
}

void  IceModelVec2::create(IceGrid &my_grid, const std::string & my_name, IceModelVecKind ghostedp,
                                     unsigned int stencil_width, int my_dof) {

  assert(m_v == NULL);

  m_dof  = my_dof;
  grid = &my_grid;

  if ((m_dof != 1) || (stencil_width > grid->config.get("grid_max_stencil_width"))) {
    m_da_stencil_width = stencil_width;
  } else {
    m_da_stencil_width = grid->config.get("grid_max_stencil_width");
  }

  // initialize the da member:
  m_da = grid->get_dm(this->m_dof, this->m_da_stencil_width);

  if (ghostedp) {
    DMCreateLocalVector(*m_da, &m_v);
  } else {
    DMCreateGlobalVector(*m_da, &m_v);
  }

  m_has_ghosts = (ghostedp == WITH_GHOSTS);
  m_name       = my_name;

  if (m_dof == 1) {
    m_metadata.push_back(NCSpatialVariable(grid->get_unit_system(),
                                           my_name, *grid));
  } else {

    for (unsigned int j = 0; j < m_dof; ++j) {
      char tmp[TEMPORARY_STRING_LENGTH];

      snprintf(tmp, TEMPORARY_STRING_LENGTH, "%s[%d]",
               m_name.c_str(), j);
      m_metadata.push_back(NCSpatialVariable(grid->get_unit_system(),
                                             tmp, *grid));
    }
  }
}

void IceModelVec2S::add(double alpha, IceModelVec &x) {
  add_2d<IceModelVec2S>(this, alpha, &x, this);
}

void IceModelVec2S::add(double alpha, const IceModelVec &x, IceModelVec &result) const {
  add_2d<IceModelVec2S>(this, alpha, &x, &result);
}

void IceModelVec2S::copy_to(IceModelVec &destination) const {
  copy_2d<IceModelVec2S>(this, &destination);
}

// IceModelVec2Stag
void IceModelVec2Stag::create(IceGrid &my_grid, const std::string &my_short_name, IceModelVecKind ghostedp,
                                        unsigned int stencil_width) {

  IceModelVec2::create(my_grid, my_short_name, ghostedp, stencil_width, m_dof);
}

//! Averages staggered grid values of a scalar field and puts them on a regular grid.
/*!
 * The current IceModelVec needs to have ghosts.
 */
void IceModelVec2Stag::staggered_to_regular(IceModelVec2S &result) const {
  IceModelVec::AccessList list(*this);
  list.add(result);

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i,j) = 0.25 * ((*this)(i,j,0) + (*this)(i,j,1)
                          + (*this)(i,j-1,1) + (*this)(i-1,j,0));
  }
}

//! \brief Averages staggered grid values of a 2D vector field (u on the
//! i-offset, v on the j-offset) and puts them on a regular grid.
/*!
 * The current IceModelVec needs to have ghosts.
 */
void IceModelVec2Stag::staggered_to_regular(IceModelVec2V &result) const {
  IceModelVec::AccessList list(*this);
  list.add(result);

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i,j).u = 0.5 * ((*this)(i-1,j,0) + (*this)(i,j,0));
    result(i,j).v = 0.5 * ((*this)(i,j-1,1) + (*this)(i,j,1));
  }
}


//! For each component, finds the maximum over all the absolute values.  Ignores ghosts.
/*!
Assumes z is allocated.
 */
void IceModelVec2Stag::absmaxcomponents(double* z) const {
  double my_z[2] = {0.0, 0.0};

  IceModelVec::AccessList list(*this);
  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    my_z[0] = PetscMax(my_z[0],PetscAbs((*this)(i,j,0)));
    my_z[1] = PetscMax(my_z[1],PetscAbs((*this)(i,j,1)));
  }

  GlobalMax(grid->com, &(my_z[0]), &(z[0]));
  GlobalMax(grid->com, &(my_z[1]), &(z[1]));
}

} // end of namespace pism
