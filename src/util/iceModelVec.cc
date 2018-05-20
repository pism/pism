// Copyright (C) 2008--2018 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#include <cassert>

#include "pism_utilities.hh"
#include "iceModelVec.hh"
#include "pism/util/io/PIO.hh"
#include "Time.hh"
#include "IceGrid.hh"
#include "ConfigInterface.hh"

#include "error_handling.hh"
#include "iceModelVec_helpers.hh"
#include "io/io_helpers.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/petscwrappers/VecScatter.hh"

namespace pism {

IceModelVec::IceModelVec() {
  m_access_counter = 0;
  m_array = NULL;

  m_da.reset();
  m_da_stencil_width = 1;
  m_dof = 1;                    // default
  m_begin_end_access_use_dof = true;

  m_has_ghosts = true;

  m_name = "unintialized variable";

  // would resize "vars", but "grid" is not initialized, and so we
  // cannot get the unit system:
  // vars.resize(dof);
  reset_attrs(0);

  m_state_counter = 0;
  m_interpolation_type = BILINEAR;

  m_zlevels.resize(1);
  m_zlevels[0] = 0.0;
}

//! \brief Get the object state counter.
/*!
 * This method returns the "revision number" of an IceModelVec.
 *
 * It can be used to determine it a field was updated and if a certain
 * computation needs to be re-done. One example is computing the smoothed bed
 * for the SIA computation, which is only necessary if the bed deformation code
 * fired.
 *
 * See also inc_state_counter().
 */
int IceModelVec::state_counter() const {
  return m_state_counter;
}

IceGrid::ConstPtr IceModelVec::grid() const {
  return m_grid;
}

unsigned int IceModelVec::ndof() const {
  return m_dof;
}

std::vector<double> IceModelVec::levels() const {
  return m_zlevels;
}


//! \brief Increment the object state counter.
/*!
 * See the documentation of get_state_counter(). This method is the
 * *only* way to manually increment the state counter. It is also
 * automatically updated by IceModelVec methods that are known to
 * change stored values.
 */
void IceModelVec::inc_state_counter() {
  m_state_counter++;
}

IceModelVec::~IceModelVec() {
  assert(m_access_counter == 0);
}

//! Returns true if create() was called and false otherwise.
bool IceModelVec::was_created() const {
  return (m_v != NULL);
}

//! Returns the grid type of an IceModelVec. (This is the way to figure out if an IceModelVec is 2D or 3D).
unsigned int IceModelVec::ndims() const {
  if (m_zlevels.size() > 1) {
    return 3;
  }

  return 2;
}

std::vector<int> IceModelVec::shape() const {

  if (ndims() == 3) {
    return {(int)m_grid->My(), (int)m_grid->Mx(), (int)levels().size()};
  } else {
    if (ndof() == 1) {
      return {(int)m_grid->My(), (int)m_grid->Mx()};
    } else {
      return {(int)m_grid->My(), (int)m_grid->Mx(), (int)ndof()};
    }
  }
}

//! Set the time independent flag for all variables corresponding to this IceModelVec instance.
/** A "time independent" IceModelVec will be saved to a NetCDF
    variable which does not depend on the "time" dimension.
 */
void IceModelVec::set_time_independent(bool flag) {
  for (unsigned int j = 0; j < m_dof; ++j) {
    m_metadata[j].set_time_independent(flag);
  }
}

//! Result: min <- min(v[j]), max <- max(v[j]).
/*!
PETSc manual correctly says "VecMin and VecMax are collective on Vec" but
GlobalMax,GlobalMin \e are needed, when m_has_ghosts==true, to get correct
values because Vecs created with DACreateLocalVector() are of type
VECSEQ and not VECMPI.  See src/trypetsc/localVecMax.c.
 */
Range IceModelVec::range() const {
  Range result;
  PetscErrorCode ierr;
  assert(m_v != NULL);

  ierr = VecMin(m_v, NULL, &result.min);
  PISM_CHK(ierr, "VecMin");

  ierr = VecMax(m_v, NULL, &result.max);
  PISM_CHK(ierr, "VecMax");

  if (m_has_ghosts) {
    // needs a reduce operation; use GlobalMin and GlobalMax;
    result.min = GlobalMin(m_grid->com, result.min);
    result.max = GlobalMax(m_grid->com, result.max);
  }
  return result;
}

/** Convert from `int` to PETSc's `NormType`.
 * 
 *
 * @param[in] input norm type as an integer
 *
 * @return norm type as PETSc's `NormType`.
 */
NormType IceModelVec::int_to_normtype(int input) const {
  assert(input == NORM_1 || input == NORM_2 || input == NORM_INFINITY);

  switch (input) {
  case NORM_1:
    return NORM_1;
  case NORM_2:
    return NORM_2;
  default:
  case NORM_INFINITY:
    return NORM_INFINITY;
  }
}

//! Computes the norm of an dof==1 IceModelVec by calling PETSc VecNorm.
/*!
See comment for range(); because local Vecs are VECSEQ, needs a reduce operation.
See src/trypetsc/localVecMax.c.

@note This method works for all IceModelVecs, including ones with
dof > 1. You might want to use norm_all() for IceModelVec2Stag,
though.
 */
double IceModelVec::norm(int n) const {
  return this->norm_all(n)[0];
}

//! Result: v <- sqrt(v), elementwise.  Calls VecSqrt(v).
/*!
Name avoids clash with sqrt() in math.h.
 */
void IceModelVec::squareroot() {
  assert(m_v != NULL);

  PetscErrorCode ierr = VecSqrtAbs(m_v);
  PISM_CHK(ierr, "VecSqrtAbs");
}


//! Result: v <- v + alpha * x. Calls VecAXPY.
void IceModelVec::add(double alpha, const IceModelVec &x) {
  assert(m_v != NULL && x.m_v != NULL);

  checkCompatibility("add", x);

  PetscErrorCode ierr = VecAXPY(m_v, alpha, x.m_v);
  PISM_CHK(ierr, "VecAXPY");

  inc_state_counter();          // mark as modified
}

//! Result: v[j] <- v[j] + alpha for all j. Calls VecShift.
void IceModelVec::shift(double alpha) {
  assert(m_v != NULL);

  PetscErrorCode ierr = VecShift(m_v, alpha);
  PISM_CHK(ierr, "VecShift");

  inc_state_counter();          // mark as modified
}

//! Result: v <- v * alpha. Calls VecScale.
void IceModelVec::scale(double alpha) {
  assert(m_v != NULL);

  PetscErrorCode ierr = VecScale(m_v, alpha);
  PISM_CHK(ierr, "VecScale");

  inc_state_counter();          // mark as modified
}

//! Copies v to a global vector 'destination'. Ghost points are discarded.
/*! This is potentially dangerous: make sure that `destination` has the same
    dimensions as the current IceModelVec.

    DMLocalToGlobalBegin/End is broken in PETSc 3.5, so we roll our
    own.
 */
void  IceModelVec::copy_to_vec(petsc::DM::Ptr destination_da, Vec destination) const {
  assert(m_v != NULL);

  // m_dof > 1 for vector, staggered grid 2D fields, etc. In this case
  // zlevels.size() == 1. For 3D fields, m_dof == 1 (all 3D fields are
  // scalar) and zlevels.size() corresponds to dof of the underlying PETSc
  // DM object. So we want the bigger of the two numbers here.
  unsigned int N = std::max((size_t)m_dof, m_zlevels.size());

  this->get_dof(destination_da, destination, 0, N);
}

//! \brief Copies data from a Vec `source` to this IceModelVec. Updates ghost
//! points if necessary.
/*!
  Unlike DMLocalToGlobalBegin/End, DMGlobalToLocalBegin/End is *not*
  broken in PETSc 3.5 (and ealier), so we can use it here.
 */
void IceModelVec::copy_from_vec(Vec source) {
  PetscErrorCode ierr;
  assert(m_v != NULL);

  if (m_has_ghosts) {
    global_to_local(m_da, source, m_v);
  } else {
    ierr = VecCopy(source, m_v);
    PISM_CHK(ierr, "VecCopy");
  }

  inc_state_counter();          // mark as modified
}


void IceModelVec::get_dof(petsc::DM::Ptr da_result, Vec result,
                          unsigned int start, unsigned int count) const {
  if (start >= m_dof) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid argument (start); got %d", start);
  }

  petsc::DMDAVecArrayDOF tmp_res(da_result, result), tmp_v(m_da, m_v);

  double
    ***result_a = static_cast<double***>(tmp_res.get()),
    ***source_a = static_cast<double***>(tmp_v.get());

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      PetscErrorCode ierr = PetscMemcpy(result_a[j][i], &source_a[j][i][start],
                                        count*sizeof(PetscScalar));
      PISM_CHK(ierr, "PetscMemcpy");
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();
}

void IceModelVec::set_dof(petsc::DM::Ptr da_source, Vec source,
                          unsigned int start, unsigned int count) {
  if (start >= m_dof) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid argument (start); got %d", start);
  }

  petsc::DMDAVecArrayDOF tmp_src(da_source, source), tmp_v(m_da, m_v);
  
  double
    ***source_a = static_cast<double***>(tmp_src.get()),
    ***result_a = static_cast<double***>(tmp_v.get());

  ParallelSection loop(m_grid->com);
  try {
    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();
      PetscErrorCode ierr = PetscMemcpy(&result_a[j][i][start], source_a[j][i],
                                        count*sizeof(PetscScalar));
      PISM_CHK(ierr, "PetscMemcpy");
    }
  } catch (...) {
    loop.failed();
  }
  loop.check();

  inc_state_counter();          // mark as modified
}

//! Result: v <- source.  Leaves metadata alone but copies values in Vec.  Uses VecCopy.
void  IceModelVec::copy_from(const IceModelVec &source) {
  PetscErrorCode ierr;
  assert(m_v != NULL && source.m_v != NULL);

  checkCompatibility("copy_from", source);

  ierr = VecCopy(source.m_v, m_v);
  PISM_CHK(ierr, "VecCopy");

  this->inc_state_counter();          // mark as modified
}

//! @brief Get the stencil width of the current IceModelVec. Returns 0
//! if ghosts are not available.
unsigned int IceModelVec::stencil_width() const {
  if (m_has_ghosts) {
    return m_da_stencil_width;
  } else {
    return 0;
  }
}

Vec IceModelVec::vec() {
  return m_v;
}

petsc::DM::Ptr IceModelVec::dm() const {
  return m_da;
}

//! Sets the variable name to `name`.
/**
 * This is the "overall" name of a field. This is **not** the same as the
 * NetCDF variable name. Use `metadata(...).set_name(...)` to set that.
 */
void IceModelVec::set_name(const std::string &name) {
  m_name = name;
}

//! @brief Get the name of an IceModelVec object.
/**
 * This is the name used to refer to this object in PISM (e.g. via the
 * Vars class), **not** the name of the corresponding NetCDF variable.
 * (The problem is that one IceModelVec instance may correspond to
 * several NetCDF variables. Use metadata(...).get_name() to get the
 * name of NetCDF variables an IceModelVec is saved to.)
 */
const std::string& IceModelVec::get_name() const {
  return m_name;
}

//! Resets most IceModelVec attributes.
void IceModelVec::reset_attrs(unsigned int N) {

  m_report_range = true;

  if (N > 0 && N < m_metadata.size()) {
    metadata(N).clear_all_strings();
    metadata(N).clear_all_doubles();
  }
}

//! Sets NetCDF attributes of an IceModelVec object.
/*! Call set_attrs("new pism_intent", "new long name", "new units", "") if a
  variable does not have a standard name. Similarly, by putting "" in an
  appropriate spot, it is possible tp leave long_name, units or pism_intent
  unmodified.

  If units != "", this also resets glaciological_units, so that they match
  internal units.
 */
void IceModelVec::set_attrs(const std::string &pism_intent,
                            const std::string &long_name,
                            const std::string &units,
                            const std::string &standard_name,
                            unsigned int N) {

  metadata(N).set_string("long_name", long_name);

  metadata(N).set_string("units", units);

  metadata(N).set_string("pism_intent", pism_intent);

  metadata(N).set_string("standard_name", standard_name);
}

//! Gets an IceModelVec from a file `nc`, interpolating onto the current grid.
/*! Stops if the variable was not found and `critical` == true.
 */
void IceModelVec::regrid_impl(const PIO &file, RegriddingFlag flag,
                              double default_value) {
  if (m_dof != 1) {
    throw RuntimeError(PISM_ERROR_LOCATION, "This method (IceModelVec::regrid_impl)"
                       " only supports IceModelVecs with dof == 1.");
  }

  bool allow_extrapolation = m_grid->ctx()->config()->get_boolean("grid.allow_extrapolation");

  if (m_has_ghosts) {
    petsc::TemporaryGlobalVec tmp(m_da);
    petsc::VecArray tmp_array(tmp);

    io::regrid_spatial_variable(metadata(0), *m_grid, file, flag,
                                m_report_range, allow_extrapolation,
                                default_value, m_interpolation_type, tmp_array.get());

    global_to_local(m_da, tmp, m_v);
  } else {
    petsc::VecArray v_array(m_v);
    io::regrid_spatial_variable(metadata(0), *m_grid,  file, flag,
                                m_report_range, allow_extrapolation,
                                default_value, m_interpolation_type, v_array.get());
  }
}

//! Reads appropriate NetCDF variable(s) into an IceModelVec.
void IceModelVec::read_impl(const PIO &nc, const unsigned int time) {

  m_grid->ctx()->log()->message(3, "  Reading %s...\n", m_name.c_str());

  if (m_dof != 1) {
    throw RuntimeError(PISM_ERROR_LOCATION, "This method (IceModelVec::read_impl) only supports"
                       " IceModelVecs with dof == 1.");
  }

  if (m_has_ghosts) {
    petsc::TemporaryGlobalVec tmp(m_da);
    petsc::VecArray tmp_array(tmp);

    io::read_spatial_variable(metadata(0), *m_grid, nc, time, tmp_array.get());

    global_to_local(m_da, tmp, m_v);
  } else {
    petsc::VecArray v_array(m_v);
    io::read_spatial_variable(metadata(0), *m_grid, nc, time, v_array.get());
  }
}

//! \brief Define variables corresponding to an IceModelVec in a file opened using `nc`.
void IceModelVec::define(const PIO &nc, IO_Type default_type) const {
  std::string order = m_grid->ctx()->config()->get_string("output.variable_order");
  for (unsigned int j = 0; j < m_dof; ++j) {
    IO_Type type = metadata(j).get_output_type();
    type = type == PISM_NAT ? default_type : type;
    io::define_spatial_variable(metadata(j), *m_grid, nc, type, order);
  }
}

//! \brief Read attributes from the corresponding variable in `nc`.
/*! Note that unlike read() and regrid(), this method does not use the standard
  name to find the variable to read attributes from.
 */
void IceModelVec::read_attributes(const std::string &filename, int N) {
  PIO nc(m_grid->com, "netcdf3", filename, PISM_READONLY); // OK to use netcdf3
  io::read_attributes(nc, metadata(N).get_name(), metadata(N));
}


//! @brief Returns a reference to the SpatialVariableMetadata object
//! containing metadata for the compoment N.
SpatialVariableMetadata& IceModelVec::metadata(unsigned int N) {
  assert(N < m_dof);
  return m_metadata[N];
}

const SpatialVariableMetadata& IceModelVec::metadata(unsigned int N) const {
  assert(N < m_dof);
  return m_metadata[N];
}

//! Writes an IceModelVec to a NetCDF file.
void IceModelVec::write_impl(const PIO &file) const {

  if (m_dof != 1) {
    throw RuntimeError(PISM_ERROR_LOCATION, "This method (IceModelVec::write_impl) only supports"
                       " IceModelVecs with dof == 1");
  }

  if (m_has_ghosts) {
    petsc::TemporaryGlobalVec tmp(m_da);

    this->copy_to_vec(m_da, tmp);

    petsc::VecArray tmp_array(tmp);

    io::write_spatial_variable(metadata(0), *m_grid,  file,
                               tmp_array.get());
  } else {
    petsc::VecArray v_array(m_v);
    io::write_spatial_variable(metadata(0), *m_grid, file,
                               v_array.get());
  }
}

//! Dumps a variable to a file, overwriting this file's contents (for debugging).
void IceModelVec::dump(const char filename[]) const {
  PIO file(m_grid->com, m_grid->ctx()->config()->get_string("output.format"),
           filename, PISM_READWRITE_CLOBBER);

  io::define_time(file, *m_grid->ctx());
  io::append_time(file, *m_grid->ctx()->config(), m_grid->ctx()->time()->current());

  define(file, PISM_DOUBLE);
  write(file);
}

//! Checks if two IceModelVecs have compatible sizes, dimensions and numbers of degrees of freedom.
void IceModelVec::checkCompatibility(const char* func, const IceModelVec &other) const {
  PetscErrorCode ierr;
  // We have to use PetscInt because of VecGetSizes below.
  PetscInt X_size, Y_size;

  if (m_dof != other.m_dof) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "IceModelVec::%s(...): operands have different numbers of degrees of freedom",
                                  func);
  }

  ierr = VecGetSize(m_v, &X_size);
  PISM_CHK(ierr, "VecGetSize");

  ierr = VecGetSize(other.m_v, &Y_size);
  PISM_CHK(ierr, "VecGetSize");

  if (X_size != Y_size) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "IceModelVec::%s(...): incompatible Vec sizes (called as %s.%s(%s))",
                                  func, m_name.c_str(), func, other.m_name.c_str());
  }
}

//! Checks if an IceModelVec is allocated and calls DAVecGetArray.
void  IceModelVec::begin_access() const {
  assert(m_v != NULL);

  if (m_access_counter < 0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "IceModelVec::begin_access(): m_access_counter < 0");
  }

  if (m_access_counter == 0) {
    PetscErrorCode ierr;
    if (m_begin_end_access_use_dof == true) {
      ierr = DMDAVecGetArrayDOF(*m_da, m_v, &m_array);
      PISM_CHK(ierr, "DMDAVecGetArrayDOF");
    } else {
      ierr = DMDAVecGetArray(*m_da, m_v, &m_array);
      PISM_CHK(ierr, "DMDAVecGetArray");
    }
  }

  m_access_counter++;
}

//! Checks if an IceModelVec is allocated and calls DAVecRestoreArray.
void  IceModelVec::end_access() const {
  PetscErrorCode ierr;
  assert(m_v != NULL);

  if (m_array == NULL) {
    throw RuntimeError(PISM_ERROR_LOCATION, "IceModelVec::end_access(): a == NULL (looks like begin_acces() was not called)");
  }

  if (m_access_counter < 0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "IceModelVec::end_access(): m_access_counter < 0");
  }


  m_access_counter--;
  if (m_access_counter == 0) {
    if (m_begin_end_access_use_dof == true) {
      ierr = DMDAVecRestoreArrayDOF(*m_da, m_v, &m_array);
      PISM_CHK(ierr, "DMDAVecRestoreArrayDOF");
    } else {
      ierr = DMDAVecRestoreArray(*m_da, m_v, &m_array);
      PISM_CHK(ierr, "DMDAVecRestoreArray");
    }
    m_array = NULL;
  }
}

//! Updates ghost points.
void  IceModelVec::update_ghosts() {
  PetscErrorCode ierr;
  if (not m_has_ghosts) {
    return;
  }

  assert(m_v != NULL);

  ierr = DMLocalToLocalBegin(*m_da, m_v, INSERT_VALUES, m_v);
  PISM_CHK(ierr, "DMLocalToLocalBegin");
  
  ierr = DMLocalToLocalEnd(*m_da, m_v, INSERT_VALUES, m_v);
  PISM_CHK(ierr, "DMLocalToLocalEnd");
}

void IceModelVec::global_to_local(petsc::DM::Ptr dm, Vec source, Vec destination) const {
  PetscErrorCode ierr;

  ierr = DMGlobalToLocalBegin(*dm, source, INSERT_VALUES, destination);
  PISM_CHK(ierr, "DMGlobalToLocalBegin");

  ierr = DMGlobalToLocalEnd(*dm, source, INSERT_VALUES, destination);
  PISM_CHK(ierr, "DMGlobalToLocalEnd");
}


//! Scatters ghost points to IceModelVec destination.
void  IceModelVec::update_ghosts(IceModelVec &destination) const {
  PetscErrorCode ierr;

  // Make sure it is allocated:
  assert(m_v != NULL);
  // Make sure "destination" has ghosts to update.
  assert(destination.m_has_ghosts);

  if (m_has_ghosts and destination.m_has_ghosts) {
    ierr = DMLocalToLocalBegin(*m_da, m_v, INSERT_VALUES, destination.vec());
    PISM_CHK(ierr, "DMLocalToLocalBegin");

    ierr = DMLocalToLocalEnd(*m_da, m_v, INSERT_VALUES, destination.vec());
    PISM_CHK(ierr, "DMLocalToLocalEnd");

    return;
  }

  if (not m_has_ghosts and destination.m_has_ghosts) {
    global_to_local(destination.dm(), m_v, destination.vec());

    return;
  }

  destination.inc_state_counter();          // mark as modified
}

//! Result: v[j] <- c for all j.
void  IceModelVec::set(const double c) {
  assert(m_v != NULL);

  PetscErrorCode ierr = VecSet(m_v,c);
  PISM_CHK(ierr, "VecSet");

  inc_state_counter();          // mark as modified
}

void IceModelVec::check_array_indices(int i, int j, unsigned int k) const {
  double ghost_width = 0;
  if (m_has_ghosts) {
    ghost_width = m_da_stencil_width;
  }
  // m_dof > 1 for vector, staggered grid 2D fields, etc. In this case
  // zlevels.size() == 1. For 3D fields, m_dof == 1 (all 3D fields are
  // scalar) and zlevels.size() corresponds to dof of the underlying PETSc
  // DM object. So we want the bigger of the two numbers here.
  unsigned int N = std::max((size_t)m_dof, m_zlevels.size());

  bool out_of_range = (i < m_grid->xs() - ghost_width) ||
    (i > m_grid->xs() + m_grid->xm() + ghost_width) ||
    (j < m_grid->ys() - ghost_width) ||
    (j > m_grid->ys() + m_grid->ym() + ghost_width) ||
    (k >= N);

  if (out_of_range) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "%s(%d, %d, %d) is out of bounds",
                                  m_name.c_str(), i, j, k);
  }

  if (m_array == NULL) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "%s: begin_access() was not called", m_name.c_str());
  }
}

//! \brief Compute parameters for 2D loop computations involving 3
//! IceModelVecs.
/*!
 * Here we assume that z is updated using a local (point-wise) computation
 * involving x and y.
 *
 * "ghosts" is the width of the stencil that can be updated locally.
 * "scatter" is false if all ghosts can be updated locally.
 */
void compute_params(const IceModelVec* const x, const IceModelVec* const y,
		    const IceModelVec* const z, int &ghosts, bool &scatter) {

  // We have 2^3=8 cases here (x,y,z having or not having ghosts).
  if (z->stencil_width() == 0) {
    // z has no ghosts; we can update everything locally
    // (This covers 4 cases.)
    ghosts = 0;
    scatter = false;
  } else if (x->stencil_width() == 0 ||
             y->stencil_width() == 0) {
    // z has ghosts, but at least one of x and y does not. we have to scatter
    // ghosts.
    // (This covers 3 cases.)
    ghosts = 0;
    scatter = true;
  } else {
    // all of x, y, z have ghosts
    // (The remaining 8-th case.)
    if (z->stencil_width() <= x->stencil_width() &&
        z->stencil_width() <= y->stencil_width()) {
      // x and y have enough ghosts to update ghosts of z locally
      ghosts = z->stencil_width();
      scatter = false;
    } else {
      // z has ghosts, but at least one of x and y doesn't have a wide enough
      // stencil
      ghosts = 0;
      scatter = true;
    }
  }
}

//! \brief Computes the norm of all components.
std::vector<double> IceModelVec::norm_all(int n) const {
  assert(m_v != NULL);

  std::vector<double> result(m_dof);

  NormType type = this->int_to_normtype(n);

  if (m_dof > 1) {
    PetscErrorCode ierr = VecStrideNormAll(m_v, type, &result[0]);
    PISM_CHK(ierr, "VecStrideNormAll");
  } else {
    PetscErrorCode ierr = VecNorm(m_v, type, &result[0]);
    PISM_CHK(ierr, "VecNorm");
  }

  if (m_has_ghosts) {
    // needs a reduce operation; use GlobalMax() if NORM_INFINITY,
    // otherwise GlobalSum; carefully in NORM_2 case
    switch (type) {
    case NORM_1_AND_2: {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "IceModelVec::norm_all(...): NORM_1_AND_2"
                                    " not implemented (called as %s.norm_all(...))",
                                    m_name.c_str());

    }
    case NORM_1: {
      for (unsigned int k = 0; k < m_dof; ++k) {
        result[k] = GlobalSum(m_grid->com, result[k]);
      }
      return result;
    }
    case NORM_2: {
      for (unsigned int k = 0; k < m_dof; ++k) {
        // undo sqrt in VecNorm before sum; sum up; take sqrt
        result[k] = sqrt(GlobalSum(m_grid->com, result[k]*result[k]));
      }
      return result;
    }
    case NORM_INFINITY: {
      for (unsigned int k = 0; k < m_dof; ++k) {
        result[k] = GlobalMax(m_grid->com, result[k]);
      }
      return result;
    }
    default: {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "IceModelVec::norm_all(...): unknown norm type"
                                    " (called as %s.norm_all(...))",
                                    m_name.c_str());
    }
    } // switch
  } else {
    return result;
  }
}

void IceModelVec::write(const std::string &filename) const {
  // We expect the file to be present and ready to write into.
  PIO nc(m_grid->com, m_grid->ctx()->config()->get_string("output.format"),
         filename, PISM_READWRITE);

  this->write(nc);
}

void IceModelVec::read(const std::string &filename, unsigned int time) {
  PIO nc(m_grid->com, "guess_mode", filename, PISM_READONLY);
  this->read(nc, time);
}

void IceModelVec::regrid(const std::string &filename, RegriddingFlag flag,
                                   double default_value) {
  PIO nc(m_grid->com, "guess_mode", filename, PISM_READONLY);

  try {
    this->regrid(nc, flag, default_value);
  } catch (RuntimeError &e) {
    e.add_context("regridding '%s' from '%s'",
                  this->get_name().c_str(), filename.c_str());
    throw;
  }
}

/** Read a field from a file, interpolating onto the current grid.
 *
 * When `flag` is set to `CRITICAL`, stop if could not find the variable
 * in the provided input file; `default_value` is ignored.
 *
 * When `flag` is set to `OPTIONAL`, fill this IceModelVec with
 * `default_value` if could not find the variable in the provided
 * input file.
 *
 * When `flag` is set to `CRITICAL_FILL_MISSING`, replace missing
 * values matching the `_FillValue` attribute with `default_value`,
 * stop if could not find the variable.
 *
 * When `flag` is set to `OPTIONAL_FILL_MISSING`, replace missing
 * values matching the `_FillValue` attribute with `default_value`;
 * fill the whole IceModelVec with `default_value` if could not find
 * the variable.
 *
 * @param nc input file
 * @param flag regridding mode, see above
 * @param default_value default value, meaning depends on the
 *        regridding mode flag
 *
 * @return 0 on success
 */
void IceModelVec::regrid(const PIO &nc, RegriddingFlag flag,
                         double default_value) {
  m_grid->ctx()->log()->message(3, "  [%s] Regridding %s...\n",
                                timestamp(m_grid->com).c_str(), m_name.c_str());
  double start_time = get_time();
  m_grid->ctx()->profiling().begin("io.regridding");
  {
    this->regrid_impl(nc, flag, default_value);
    inc_state_counter();          // mark as modified
  }
  m_grid->ctx()->profiling().end("io.regridding");
  double
    end_time   = get_time(),
    time_spent = end_time - start_time;

  if (time_spent > 1.0) {
    m_grid->ctx()->log()->message(3, "  done in %f seconds.\n", time_spent);
  } else {
    m_grid->ctx()->log()->message(3, "  done.\n");
  }
}

void IceModelVec::read(const PIO &nc, const unsigned int time) {
  this->read_impl(nc, time);
  inc_state_counter();          // mark as modified
}

void IceModelVec::write(const PIO &nc) const {
  define(nc);

  m_grid->ctx()->log()->message(3, "  [%s] Writing %s...",
                               timestamp(m_grid->com).c_str(),
                               m_name.c_str());

  double start_time = get_time();
  write_impl(nc);
  double end_time = get_time();

  const double
    time_spent = end_time - start_time,
    megabyte = pow(2, 20),
    mb_double = sizeof(double) * size() / megabyte,
    mb_float =  sizeof(float) * size() / megabyte;

  std::string timestamp = pism::timestamp(m_grid->com);
  std::string spacer(timestamp.size(), ' ');
  if (time_spent > 1) {
    m_grid->ctx()->log()->message(3,
                                  "\n"
                                  "  [%s] Done writing %s (%f Mb double, %f Mb float)\n"
                                  "   %s  in %f seconds (%f minutes).\n"
                                  "   %s  Effective throughput: double: %f Mb/s, float: %f Mb/s.\n",
                                  timestamp.c_str(), m_name.c_str(), mb_double, mb_float,
                                  spacer.c_str(), time_spent, time_spent / 60.0,
                                  spacer.c_str(),
                                  mb_double / time_spent, mb_float / time_spent);
  } else {
    m_grid->ctx()->log()->message(3, " done.\n");
  }
}

AccessList::AccessList() {
  // empty
}

AccessList::~AccessList() {
  while (not m_vecs.empty()) {
    try {
      m_vecs.back()->end_access();
      m_vecs.pop_back();
    } catch (...) {
      handle_fatal_errors(MPI_COMM_SELF);
    }
  }
}

AccessList::AccessList(std::initializer_list<const PetscAccessible *> vecs) {
  for (auto j : vecs) {
    assert(j != nullptr);
    add(*j);
  }
}

AccessList::AccessList(const PetscAccessible &vec) {
  add(vec);
}

void AccessList::add(const PetscAccessible &vec) {
  vec.begin_access();
  m_vecs.push_back(&vec);
}

void AccessList::add(const std::vector<const PetscAccessible*> vecs) {
  for (auto v : vecs) {
    assert(v != nullptr);
    add(*v);
  }
}

//! Return the total number of elements in the *owned* part of an array.
size_t IceModelVec::size() const {
  // m_dof > 1 for vector, staggered grid 2D fields, etc. In this case
  // zlevels.size() == 1. For 3D fields, m_dof == 1 (all 3D fields are
  // scalar) and zlevels.size() corresponds to dof of the underlying PETSc
  // DM object.

  size_t
    Mx = m_grid->Mx(),
    My = m_grid->My(),
    Mz = m_zlevels.size(),
    dof = m_dof;

  return Mx * My * Mz * dof;
}

struct VecAndScatter {
  VecScatter scatter;
  Vec v;
};

/*!
 * Allocate the scatter from a part of a parallel Vec to a target rank.
 *
 * The caller is responsible for de-allocating both the scatter and the target Vec.
 */
VecAndScatter scatter_part(Vec v_in, int start, int length, int target_rank) {
  PetscErrorCode ierr;
  int rank;
  VecAndScatter result;
  IS is;

  MPI_Comm_rank(PetscObjectComm((PetscObject)v_in), &rank);

  if (rank != target_rank) {
    length = 0;
  }

  ierr = VecCreateSeq(PETSC_COMM_SELF, length, &result.v);
  PISM_CHK(ierr, "VecCreateSeq");

  ierr = ISCreateStride(PETSC_COMM_SELF, length, start, 1, &is);
  PISM_CHK(ierr, "ISCreateStride");

  ierr = VecScatterCreate(v_in, is, result.v, NULL, &result.scatter);
  PISM_CHK(ierr, "VecScatterCreate");

  ierr = ISDestroy(&is);
  PISM_CHK(ierr, "ISDestroy");

  return result;
}

/*!
 * Allocate a natural Vec for a given DM.
 *
 * The caller is responsible for de-allocating the Vec returned by this function.
 */
Vec get_natural_work(DM dm) {
  PetscErrorCode ierr;
  Vec result;

  ierr = PetscObjectQuery((PetscObject)dm, "natural_work", (PetscObject*)&result);
  PISM_CHK(ierr, "PetscObjectQuery");

  if (result == NULL) {
    Vec v = NULL;
    ierr = DMDACreateNaturalVector(dm, &v);
    PISM_CHK(ierr, "DMDACreateNaturalVector");

    ierr = PetscObjectCompose((PetscObject)dm, "natural_work", (PetscObject)(v));
    PISM_CHK(ierr, "PetscObjectCompose");

    result = v;

    ierr = VecDestroy(&v);
    PISM_CHK(ierr, "VecDestroy");
  }

  return result;
}

/*!
 * Given a DM, allocate a rank 0 target Vec that can be used to gather a part of a
 * "global" Vec on rank 0. Arguments "start" and "length" define the part in question.
 *
 * The caller is responsible for de-allocating the Vec returned by this function.
 */
Vec proc0_copy(DM dm, int start, int length) {
  Vec v_proc0 = NULL;
  PetscErrorCode ierr = 0;

  ierr = PetscObjectQuery((PetscObject)dm, "v_proc0", (PetscObject*)&v_proc0);
  PISM_CHK(ierr, "PetscObjectQuery");
                                                                                          ;
  if (v_proc0 == NULL) {

    // natural_work will be destroyed at the end of scope, but it will
    // only decrement the reference counter incremented by
    // PetscObjectCompose below.
    auto natural_work = get_natural_work(dm);

    // scatter_to_zero will be destroyed at the end of scope, but it
    // will only decrement the reference counter incremented by
    // PetscObjectCompose below.
    auto vs = scatter_part(natural_work, start, length, 0);

    // this increments the reference counter of scatter_to_zero
    ierr = PetscObjectCompose((PetscObject)dm, "scatter_to_zero",
                              (PetscObject)(vs.scatter));
    PISM_CHK(ierr, "PetscObjectCompose");

    // this increments the reference counter of v_proc0
    ierr = PetscObjectCompose((PetscObject)dm, "v_proc0",
                              (PetscObject)vs.v);
    PISM_CHK(ierr, "PetscObjectCompose");

    VecScatterDestroy(&vs.scatter);

    // We DO NOT call VecDestroy(v_proc0): the petsc::Vec wrapper will
    // take care of this.
    return vs.v;
  }
  return v_proc0;
}

/*! Allocate a copy on processor zero and the scatter needed to move data.
 */
petsc::Vec::Ptr IceModelVec::allocate_proc0_copy() const {
  PetscErrorCode ierr;
  Vec v_proc0 = NULL;
  Vec result = NULL;

  ierr = PetscObjectQuery((PetscObject)m_da->get(), "v_proc0", (PetscObject*)&v_proc0);
  PISM_CHK(ierr, "PetscObjectQuery")
                                                                                          ;
  if (v_proc0 == NULL) {

    // natural_work will be destroyed at the end of scope, but it will
    // only decrement the reference counter incremented by
    // PetscObjectCompose below.
    petsc::Vec natural_work;
    // create a work vector with natural ordering:
    ierr = DMDACreateNaturalVector(*m_da, natural_work.rawptr());
    PISM_CHK(ierr, "DMDACreateNaturalVector");

    // this increments the reference counter of natural_work
    ierr = PetscObjectCompose((PetscObject)m_da->get(), "natural_work",
                              (PetscObject)((::Vec)natural_work));
    PISM_CHK(ierr, "PetscObjectCompose");

    // scatter_to_zero will be destroyed at the end of scope, but it
    // will only decrement the reference counter incremented by
    // PetscObjectCompose below.
    petsc::VecScatter scatter_to_zero;

    // initialize the scatter to processor 0 and create storage on processor 0
    ierr = VecScatterCreateToZero(natural_work, scatter_to_zero.rawptr(),
                                  &v_proc0);
    PISM_CHK(ierr, "VecScatterCreateToZero");

    // this increments the reference counter of scatter_to_zero
    ierr = PetscObjectCompose((PetscObject)m_da->get(), "scatter_to_zero",
                              (PetscObject)((::VecScatter)scatter_to_zero));
    PISM_CHK(ierr, "PetscObjectCompose");

    // this increments the reference counter of v_proc0
    ierr = PetscObjectCompose((PetscObject)m_da->get(), "v_proc0",
                              (PetscObject)v_proc0);
    PISM_CHK(ierr, "PetscObjectCompose");

    // We DO NOT call VecDestroy(v_proc0): the petsc::Vec wrapper will
    // take care of this.
    result = v_proc0;
  } else {
    // We DO NOT call VecDestroy(result): the petsc::Vec wrapper will
    // take care of this.
    ierr = VecDuplicate(v_proc0, &result);
    PISM_CHK(ierr, "VecDuplicate");
  }
  return petsc::Vec::Ptr(new petsc::Vec(result));
}

void IceModelVec::put_on_proc0(Vec parallel, Vec onp0) const {
  PetscErrorCode ierr = 0;
  VecScatter scatter_to_zero = NULL;
  Vec natural_work = NULL;

  ierr = PetscObjectQuery((PetscObject)m_da->get(), "scatter_to_zero",
                          (PetscObject*)&scatter_to_zero);
  PISM_CHK(ierr, "PetscObjectQuery");

  ierr = PetscObjectQuery((PetscObject)m_da->get(), "natural_work",
                          (PetscObject*)&natural_work);
  PISM_CHK(ierr, "PetscObjectQuery");

  if (natural_work == NULL || scatter_to_zero == NULL) {
    throw RuntimeError(PISM_ERROR_LOCATION, "call allocate_proc0_copy() before calling put_on_proc0");
  }

  ierr = DMDAGlobalToNaturalBegin(*m_da, parallel, INSERT_VALUES, natural_work);
  PISM_CHK(ierr, "DMDAGlobalToNaturalBegin");

  ierr = DMDAGlobalToNaturalEnd(*m_da, parallel, INSERT_VALUES, natural_work);
  PISM_CHK(ierr, "DMDAGlobalToNaturalEnd");

  ierr = VecScatterBegin(scatter_to_zero, natural_work, onp0,
                         INSERT_VALUES, SCATTER_FORWARD);
  PISM_CHK(ierr, "VecScatterBegin");

  ierr = VecScatterEnd(scatter_to_zero, natural_work, onp0,
                       INSERT_VALUES, SCATTER_FORWARD);
  PISM_CHK(ierr, "VecScatterEnd");
}


//! Puts a local IceModelVec2S on processor 0.
void IceModelVec::put_on_proc0(Vec onp0) const {
  if (m_has_ghosts) {
    petsc::TemporaryGlobalVec tmp(m_da);
    this->copy_to_vec(m_da, tmp);
    put_on_proc0(tmp, onp0);
  } else {
    put_on_proc0(m_v, onp0);
  }
}

void IceModelVec::get_from_proc0(Vec onp0, Vec parallel) {
  PetscErrorCode ierr;

  VecScatter scatter_to_zero = NULL;
  Vec natural_work = NULL;
  ierr = PetscObjectQuery((PetscObject)m_da->get(), "scatter_to_zero",
                          (PetscObject*)&scatter_to_zero);
  PISM_CHK(ierr, "PetscObjectQuery");

  ierr = PetscObjectQuery((PetscObject)m_da->get(), "natural_work",
                          (PetscObject*)&natural_work);
  PISM_CHK(ierr, "PetscObjectQuery");

  if (natural_work == NULL || scatter_to_zero == NULL) {
    throw RuntimeError(PISM_ERROR_LOCATION, "call allocate_proc0_copy() before calling get_from_proc0");
  }

  ierr = VecScatterBegin(scatter_to_zero, onp0, natural_work,
                         INSERT_VALUES, SCATTER_REVERSE);
  PISM_CHK(ierr, "VecScatterBegin");

  ierr = VecScatterEnd(scatter_to_zero, onp0, natural_work,
                       INSERT_VALUES, SCATTER_REVERSE);
  PISM_CHK(ierr, "VecScatterEnd");

  ierr = DMDANaturalToGlobalBegin(*m_da, natural_work, INSERT_VALUES, parallel);
  PISM_CHK(ierr, "DMDANaturalToGlobalBegin");

  ierr = DMDANaturalToGlobalEnd(*m_da, natural_work, INSERT_VALUES, parallel);
  PISM_CHK(ierr, "DMDANaturalToGlobalEnd");
}

//! Gets a local IceModelVec2 from processor 0.
void IceModelVec::get_from_proc0(Vec onp0) {
  if (m_has_ghosts) {
    petsc::TemporaryGlobalVec tmp(m_da);
    get_from_proc0(onp0, tmp);
    this->copy_from_vec(tmp);
  } else {
    get_from_proc0(onp0, m_v);
  }
  inc_state_counter();          // mark as modified
}

void convert_vec(Vec v, units::System::Ptr system,
                 const std::string &spec1, const std::string &spec2) {
  units::Converter c(system, spec1, spec2);

  // has to be a PetscInt because of the VecGetLocalSize() call
  PetscInt data_size = 0;
  PetscErrorCode ierr = VecGetLocalSize(v, &data_size);
  PISM_CHK(ierr, "VecGetLocalSize");

  petsc::VecArray data(v);
  c.convert_doubles(data.get(), data_size);
}

bool set_contains(const std::set<std::string> &S, const IceModelVec &field) {
  // Note that this uses IceModelVec::get_name() and not IceModelVec::metadata() and
  // VariableMetadata::get_name(): this is used to check if a possibly multi-variable field was
  // requested.
  return member(field.get_name(), S);
}

} // end of namespace pism
