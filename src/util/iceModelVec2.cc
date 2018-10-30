// Copyright (C) 2008--2018 Ed Bueler and Constantine Khroulev
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

#include <cassert>

#include <memory>
using std::dynamic_pointer_cast;

#include <petscdraw.h>
#include <petscdmda.h>

#include "pism/util/io/PIO.hh"
#include "iceModelVec.hh"
#include "IceGrid.hh"
#include "ConfigInterface.hh"

#include "error_handling.hh"
#include "iceModelVec_helpers.hh"

#include "pism_utilities.hh"
#include "io/io_helpers.hh"
#include "pism/util/Logger.hh"

namespace pism {

// this file contains methods for derived classes IceModelVec2S and IceModelVec2Int

// methods for base class IceModelVec are in "iceModelVec.cc"

IceModelVec2::IceModelVec2()
  : IceModelVec() {
  // empty
}

IceModelVec2::Ptr IceModelVec2::To2D(IceModelVec::Ptr input) {
  IceModelVec2::Ptr result = dynamic_pointer_cast<IceModelVec2,IceModelVec>(input);
  if (not (bool)result) {
    throw RuntimeError(PISM_ERROR_LOCATION, "dynamic cast failure");
  }
  return result;
}

IceModelVec2V::~IceModelVec2V() {
  // empty
}

IceModelVec2S::Ptr IceModelVec2S::To2DScalar(IceModelVec::Ptr input) {
  IceModelVec2S::Ptr result = dynamic_pointer_cast<IceModelVec2S,IceModelVec>(input);
  if (not (bool)result) {
    throw RuntimeError(PISM_ERROR_LOCATION, "dynamic cast failure");
  }
  return result;
}

IceModelVec2S::IceModelVec2S() {
  m_begin_end_access_use_dof = false;
}

IceModelVec2S::IceModelVec2S(IceGrid::ConstPtr grid, const std::string &name,
                             IceModelVecKind ghostedp, int width) {
  m_begin_end_access_use_dof = false;

  create(grid, name, ghostedp, width);
}


void  IceModelVec2S::create(IceGrid::ConstPtr grid, const std::string &name, IceModelVecKind ghostedp, int width) {
  assert(m_v == NULL);
  IceModelVec2::create(grid, name, ghostedp, width, m_dof);
}


double** IceModelVec2S::get_array() {
  begin_access();
  return static_cast<double**>(m_array);
}

//! Sets an IceModelVec2 to the magnitude of a 2D vector field with components `v_x` and `v_y`.
/*! Computes the magnitude \b pointwise, so any of v_x, v_y and the IceModelVec
  this is called on can be the same.

  Does not communicate.
 */
void IceModelVec2S::set_to_magnitude(const IceModelVec2S &v_x,
                                     const IceModelVec2S &v_y) {
  IceModelVec::AccessList list{this, &v_x, &v_y};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*this)(i,j) = sqrt(PetscSqr(v_x(i,j)) + PetscSqr(v_y(i,j)));
  }

  inc_state_counter();          // mark as modified
  
}

void IceModelVec2S::set_to_magnitude(const IceModelVec2V &input) {
  IceModelVec::AccessList list{this, &input};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    (*this)(i, j) = input(i, j).magnitude();
  }
}

//! Masks out all the areas where \f$ M \le 0 \f$ by setting them to `fill`. 
void IceModelVec2S::mask_by(const IceModelVec2S &M, double fill) {
  IceModelVec::AccessList list{this, &M};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (M(i,j) <= 0.0) {
      (*this)(i,j) = fill;
    }
  }

  inc_state_counter();          // mark as modified
}

void IceModelVec2::write_impl(const PIO &file) const {

  assert(m_v != NULL);

  // The simplest case:
  if ((m_dof == 1) and (not m_has_ghosts)) {
    IceModelVec::write_impl(file);
    return;
  }

  // Get the dof=1, stencil_width=0 DMDA (components are always scalar
  // and we just need a global Vec):
  petsc::DM::Ptr da2 = m_grid->get_dm(1, 0);

  // a temporary one-component vector, distributed across processors
  // the same way v is
  petsc::TemporaryGlobalVec tmp(da2);

  Logger::ConstPtr log = m_grid->ctx()->log();
  log->message(4, "  Writing %s...\n", m_name.c_str());

  for (unsigned int j = 0; j < m_dof; ++j) {
    IceModelVec2::get_dof(da2, tmp, j);

    petsc::VecArray tmp_array(tmp);
    io::write_spatial_variable(m_metadata[j], *m_grid, file,
                           tmp_array.get());
  }
}

void IceModelVec2::read_impl(const PIO &nc, const unsigned int time) {

  if ((m_dof == 1) and (not m_has_ghosts)) {
    IceModelVec::read_impl(nc, time);
    return;
  }

  Logger::ConstPtr log = m_grid->ctx()->log();
  log->message(4, "  Reading %s...\n", m_name.c_str());

  assert(m_v != NULL);

  // Get the dof=1, stencil_width=0 DMDA (components are always scalar
  // and we just need a global Vec):
  petsc::DM::Ptr da2 = m_grid->get_dm(1, 0);

  // a temporary one-component vector, distributed across processors
  // the same way v is
  petsc::TemporaryGlobalVec tmp(da2);

  for (unsigned int j = 0; j < m_dof; ++j) {

    {
      petsc::VecArray tmp_array(tmp);
      io::read_spatial_variable(m_metadata[j], *m_grid, nc, time, tmp_array.get());
    }

    IceModelVec2::set_dof(da2, tmp, j);
  }
  
  // The calls above only set the values owned by a processor, so we need to
  // communicate if m_has_ghosts == true:
  if (m_has_ghosts) {
    update_ghosts();
  }
}

void IceModelVec2::regrid_impl(const PIO &file, RegriddingFlag flag,
                               double default_value) {
  if ((m_dof == 1) and (not m_has_ghosts)) {
    IceModelVec::regrid_impl(file, flag, default_value);
    return;
  }

  // Get the dof=1, stencil_width=0 DMDA (components are always scalar
  // and we just need a global Vec):
  petsc::DM::Ptr da2 = m_grid->get_dm(1, 0);

  // a temporary one-component vector, distributed across processors
  // the same way v is
  petsc::TemporaryGlobalVec tmp(da2);

  const bool allow_extrapolation = m_grid->ctx()->config()->get_boolean("grid.allow_extrapolation");

  for (unsigned int j = 0; j < m_dof; ++j) {
    {
      petsc::VecArray tmp_array(tmp);
      io::regrid_spatial_variable(m_metadata[j], *m_grid,  file, flag,
                                  m_report_range, allow_extrapolation,
                                  default_value, m_interpolation_type, tmp_array.get());
    }

    IceModelVec2::set_dof(da2, tmp, j);
  }

  // The calls above only set the values owned by a processor, so we need to
  // communicate if m_has_ghosts == true:
  if (m_has_ghosts == true) {
    update_ghosts();
  }
}

//! \brief View a 2D field.
void IceModelVec2::view(int viewer_size) const {
  petsc::Viewer::Ptr viewers[2];

  if (m_dof > 2) {
    throw RuntimeError(PISM_ERROR_LOCATION, "dof > 2 is not supported");
  }

  for (unsigned int j = 0; j < std::min(m_dof, 2U); ++j) {
    std::string
      c_name              = m_metadata[j].get_name(),
      long_name           = m_metadata[j].get_string("long_name"),
      glaciological_units = m_metadata[j].get_string("glaciological_units"),
      title               = long_name + " (" + glaciological_units + ")";

    if (not m_map_viewers[c_name]) {
      m_map_viewers[c_name].reset(new petsc::Viewer(m_grid->com, title, viewer_size,
                                                    m_grid->Lx(), m_grid->Ly()));
    }

    viewers[j] = m_map_viewers[c_name];
  }

  view(viewers[0], viewers[1]);
}

//! \brief View a 2D vector field using existing PETSc viewers.
void IceModelVec2::view(petsc::Viewer::Ptr v1, petsc::Viewer::Ptr v2) const {
  PetscErrorCode ierr;

  petsc::Viewer::Ptr viewers[2] = {v1, v2};

  // Get the dof=1, stencil_width=0 DMDA (components are always scalar
  // and we just need a global Vec):
  petsc::DM::Ptr da2 = m_grid->get_dm(1, 0);

  petsc::TemporaryGlobalVec tmp(da2);

  for (unsigned int i = 0; i < std::min(m_dof, 2U); ++i) {
    std::string
      long_name           = m_metadata[i].get_string("long_name"),
      units               = m_metadata[i].get_string("units"),
      glaciological_units = m_metadata[i].get_string("glaciological_units"),
      title               = long_name + " (" + glaciological_units + ")";

    if (not (bool)viewers[i]) {
      continue;
    }

    PetscViewer v = *viewers[i].get();

    PetscDraw draw = NULL;
    ierr = PetscViewerDrawGetDraw(v, 0, &draw);
    PISM_CHK(ierr, "PetscViewerDrawGetDraw");

    ierr = PetscDrawSetTitle(draw, title.c_str());
    PISM_CHK(ierr, "PetscDrawSetTitle");

    IceModelVec2::get_dof(da2, tmp, i);

    convert_vec(tmp, m_metadata[i].unit_system(),
                units, glaciological_units);

    double bounds[2] = {0.0, 0.0};
    ierr = VecMin(tmp, NULL, &bounds[0]); PISM_CHK(ierr, "VecMin");
    ierr = VecMax(tmp, NULL, &bounds[1]); PISM_CHK(ierr, "VecMax");

    if (bounds[0] == bounds[1]) {
      bounds[0] = -1.0;
      bounds[1] =  1.0;
    }

    ierr = PetscViewerDrawSetBounds(v, 1, bounds);
    PISM_CHK(ierr, "PetscViewerDrawSetBounds");

    ierr = VecView(tmp, v);
    PISM_CHK(ierr, "VecView");
  }
}

//! \brief Returns the x-derivative at i,j approximated using centered finite
//! differences.
double IceModelVec2S::diff_x(int i, int j) const {
  return ((*this)(i + 1,j) - (*this)(i - 1,j)) / (2 * m_grid->dx());
}

//! \brief Returns the y-derivative at i,j approximated using centered finite
//! differences.
double IceModelVec2S::diff_y(int i, int j) const {
  return ((*this)(i,j + 1) - (*this)(i,j - 1)) / (2 * m_grid->dy());
}


//! \brief Returns the x-derivative at East staggered point i+1/2,j approximated 
//! using centered (obvious) finite differences.
double IceModelVec2S::diff_x_stagE(int i, int j) const {
  return ((*this)(i+1,j) - (*this)(i,j)) / (m_grid->dx());
}

//! \brief Returns the y-derivative at East staggered point i+1/2,j approximated 
//! using centered finite differences.
double IceModelVec2S::diff_y_stagE(int i, int j) const {
  return ((*this)(i+1,j+1) + (*this)(i,j+1)
           - (*this)(i+1,j-1) - (*this)(i,j-1)) / (4* m_grid->dy());
}

//! \brief Returns the x-derivative at North staggered point i,j+1/2 approximated 
//! using centered finite differences.
double IceModelVec2S::diff_x_stagN(int i, int j) const {
  return ((*this)(i+1,j+1) + (*this)(i+1,j)
           - (*this)(i-1,j+1) - (*this)(i-1,j)) / (4* m_grid->dx());
}

//! \brief Returns the y-derivative at North staggered point i,j+1/2 approximated 
//! using centered (obvious) finite differences.
double IceModelVec2S::diff_y_stagN(int i, int j) const {
  return ((*this)(i,j+1) - (*this)(i,j)) / (m_grid->dy());
}


//! \brief Returns the x-derivative at i,j approximated using centered finite
//! differences. Respects grid periodicity and uses one-sided FD at grid edges
//! if necessary.
double IceModelVec2S::diff_x_p(int i, int j) const {
  if (m_grid->periodicity() & X_PERIODIC) {
    return diff_x(i,j);
  }
  
  if (i == 0) {
    return ((*this)(i + 1,j) - (*this)(i,j)) / (m_grid->dx());
  } else if (i == (int)m_grid->Mx() - 1) {
    return ((*this)(i,j) - (*this)(i - 1,j)) / (m_grid->dx());
  } else {
    return diff_x(i,j);
 }
}

//! \brief Returns the y-derivative at i,j approximated using centered finite
//! differences. Respects grid periodicity and uses one-sided FD at grid edges
//! if necessary.
double IceModelVec2S::diff_y_p(int i, int j) const {
  if (m_grid->periodicity() & Y_PERIODIC) {
    return diff_y(i,j);
  }
  
  if (j == 0) {
    return ((*this)(i,j + 1) - (*this)(i,j)) / (m_grid->dy());
  } else if (j == (int)m_grid->My() - 1) {
    return ((*this)(i,j) - (*this)(i,j - 1)) / (m_grid->dy());
  } else {
    return diff_y(i,j);
  }
}

//! Sums up all the values in an IceModelVec2S object. Ignores ghosts.
/*! Avoids copying to a "global" vector.
 */
double IceModelVec2S::sum() const {
  double result = 0;

  IceModelVec::AccessList list(*this);
  
  // sum up the local part:
  for (Points p(*m_grid); p; p.next()) {
    result += (*this)(p.i(), p.j());
  }

  // find the global sum:
  return GlobalSum(m_grid->com, result);
}


//! Finds maximum over all the values in an IceModelVec2S object.  Ignores ghosts.
double IceModelVec2S::max() const {
  IceModelVec::AccessList list(*this);

  double result = (*this)(m_grid->xs(),m_grid->ys());
  for (Points p(*m_grid); p; p.next()) {
    result = std::max(result,(*this)(p.i(), p.j()));
  }

  return GlobalMax(m_grid->com, result);
}


//! Finds maximum over all the absolute values in an IceModelVec2S object.  Ignores ghosts.
double IceModelVec2S::absmax() const {

  double result = 0.0;

  IceModelVec::AccessList list(*this);
  for (Points p(*m_grid); p; p.next()) {
    result = std::max(result, fabs((*this)(p.i(), p.j())));
  }

  return GlobalMax(m_grid->com, result);
}


//! Finds minimum over all the values in an IceModelVec2S object.  Ignores ghosts.
double IceModelVec2S::min() const {
  IceModelVec::AccessList list(*this);

  double result = (*this)(m_grid->xs(), m_grid->ys());
  for (Points p(*m_grid); p; p.next()) {
    result = std::min(result,(*this)(p.i(), p.j()));
  }

  return GlobalMin(m_grid->com, result);
}


// IceModelVec2

void IceModelVec2::get_component(unsigned int n, IceModelVec2S &result) const {

  IceModelVec2::get_dof(result.dm(), result.m_v, n);
}

void IceModelVec2::set_component(unsigned int n, const IceModelVec2S &source) {

  IceModelVec2::set_dof(source.dm(), source.m_v, n);
}

void IceModelVec2::create(IceGrid::ConstPtr grid, const std::string & name,
                           IceModelVecKind ghostedp,
                           unsigned int stencil_width, int dof) {
  PetscErrorCode ierr;
  assert(m_v == NULL);

  m_dof  = dof;
  m_grid = grid;

  if ((m_dof != 1) || (stencil_width > m_grid->ctx()->config()->get_double("grid.max_stencil_width"))) {
    m_da_stencil_width = stencil_width;
  } else {
    m_da_stencil_width = m_grid->ctx()->config()->get_double("grid.max_stencil_width");
  }

  // initialize the da member:
  m_da = m_grid->get_dm(this->m_dof, this->m_da_stencil_width);

  if (ghostedp) {
    ierr = DMCreateLocalVector(*m_da, m_v.rawptr());
    PISM_CHK(ierr, "DMCreateLocalVector");
  } else {
    ierr = DMCreateGlobalVector(*m_da, m_v.rawptr());
    PISM_CHK(ierr, "DMCreateGlobalVector");
  }

  m_has_ghosts = (ghostedp == WITH_GHOSTS);
  m_name       = name;

  if (m_dof == 1) {
    m_metadata.push_back(SpatialVariableMetadata(m_grid->ctx()->unit_system(),
                                                 name));
  } else {

    for (unsigned int j = 0; j < m_dof; ++j) {
      m_metadata.push_back(SpatialVariableMetadata(m_grid->ctx()->unit_system(),
                                                   pism::printf("%s[%d]", m_name.c_str(), j)));
    }
  }
}

void IceModelVec2S::add(double alpha, const IceModelVec &x) {
  add_2d<IceModelVec2S>(this, alpha, &x, this);
}

void IceModelVec2S::add(double alpha, const IceModelVec &x, IceModelVec &result) const {
  add_2d<IceModelVec2S>(this, alpha, &x, &result);
}

void IceModelVec2S::copy_from(const IceModelVec &source) {
  copy_2d<IceModelVec2S>(&source, this);
}

// IceModelVec2Stag

IceModelVec2Stag::IceModelVec2Stag()
  : IceModelVec2() {
  m_dof = 2;
  m_begin_end_access_use_dof = true;
}

IceModelVec2Stag::Ptr IceModelVec2Stag::ToStaggered(IceModelVec::Ptr input) {
  IceModelVec2Stag::Ptr result = dynamic_pointer_cast<IceModelVec2Stag,IceModelVec>(input);
  if (not (bool)result) {
    throw RuntimeError(PISM_ERROR_LOCATION, "dynamic cast failure");
  }
  return result;
}

IceModelVec2Stag::IceModelVec2Stag(IceGrid::ConstPtr grid, const std::string &name,
                                   IceModelVecKind ghostedp,
                                   unsigned int stencil_width)
  : IceModelVec2() {
  m_dof = 2;
  m_begin_end_access_use_dof = true;

  create(grid, name, ghostedp, stencil_width);
}

void IceModelVec2Stag::create(IceGrid::ConstPtr grid, const std::string &name,
                              IceModelVecKind ghostedp,
                              unsigned int stencil_width) {

  IceModelVec2::create(grid, name, ghostedp, stencil_width, m_dof);
}

//! Averages staggered grid values of a scalar field and puts them on a regular grid.
/*!
 * The current IceModelVec needs to have ghosts.
 */
void IceModelVec2Stag::staggered_to_regular(IceModelVec2S &result) const {
  IceModelVec::AccessList list{this, &result};

  for (Points p(*m_grid); p; p.next()) {
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
  IceModelVec::AccessList list{this, &result};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    result(i,j).u = 0.5 * ((*this)(i-1,j,0) + (*this)(i,j,0));
    result(i,j).v = 0.5 * ((*this)(i,j-1,1) + (*this)(i,j,1));
  }
}


//! For each component, finds the maximum over all the absolute values.  Ignores ghosts.
std::vector<double> IceModelVec2Stag::absmaxcomponents() const {
  std::vector<double> z(2, 0.0);

  IceModelVec::AccessList list(*this);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    z[0] = std::max(z[0], fabs((*this)(i, j, 0)));
    z[1] = std::max(z[1], fabs((*this)(i, j, 1)));
  }

  z[0] = GlobalMax(m_grid->com, z[0]);
  z[1] = GlobalMax(m_grid->com, z[1]);

  return z;
}

IceModelVec2Int::IceModelVec2Int()
  : IceModelVec2S() {
  m_interpolation_type = NEAREST;
}


IceModelVec2Int::IceModelVec2Int(IceGrid::ConstPtr grid, const std::string &name,
                                 IceModelVecKind ghostedp, int width)
  : IceModelVec2S(grid, name, ghostedp, width) {
  m_interpolation_type = NEAREST;
}


} // end of namespace pism
