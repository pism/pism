// Copyright (C) 2008--2015 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#include <gsl/gsl_math.h>       // gsl_isnan()
#include <cassert>

#include "pism_const.hh"
#include "iceModelVec.hh"
#include "PIO.hh"
#include "PISMTime.hh"
#include "IceGrid.hh"
#include "LocalInterpCtx.hh"
#include "PISMConfig.hh"


#include "error_handling.hh"
#include "iceModelVec_helpers.hh"

namespace pism {

IceModelVec::IceModelVec() {
  m_access_counter = 0;
  array = NULL;

  m_da.reset();
  m_da_stencil_width = 1;
  m_dof = 1;                    // default
  begin_end_access_use_dof = true;

  m_grid = NULL;

  m_has_ghosts = true;

  m_n_levels = 1;
  m_name = "unintialized variable";

  // would resize "vars", but "grid" is not initialized, and so we
  // cannot get the unit system:
  // vars.resize(dof);
  reset_attrs(0);

  m_state_counter = 0;

  m_v = NULL;

  zlevels.resize(1);
  zlevels[0] = 0.0;
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
int IceModelVec::get_state_counter() const {
  return m_state_counter;
}

const IceGrid* IceModelVec::get_grid() const {
  return m_grid;
}

unsigned int IceModelVec::get_ndof() const {
  return m_dof;
}

int IceModelVec::nlevels() const {
  return m_n_levels;
}

std::vector<double> IceModelVec::get_levels() const {
  return zlevels;
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
  destroy();
}

//! Returns true if create() was called and false otherwise.
bool IceModelVec::was_created() const {
  return (m_v != NULL);
}

//! Returns the grid type of an IceModelVec. (This is the way to figure out if an IceModelVec is 2D or 3D).
unsigned int IceModelVec::get_ndims() const {
  if (zlevels.size() > 1) {
    return 3;
  }

  return 2;
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

//! \brief De-allocates an IceModelVec object.
void  IceModelVec::destroy() {
  PetscErrorCode ierr;

  if (m_v != NULL) {
    ierr = VecDestroy(&m_v);
    PISM_PETSC_CHK(ierr, "VecDestroy");
    m_v = NULL;
  }

  assert(m_access_counter == 0);
}

//! Result: min <- min(v[j]), max <- max(v[j]).
/*!
PETSc manual correctly says "VecMin and VecMax are collective on Vec" but
GlobalMax,GlobalMin \e are needed, when m_has_ghosts==true, to get correct
values because Vecs created with DACreateLocalVector() are of type
VECSEQ and not VECMPI.  See src/trypetsc/localVecMax.c.
 */
void IceModelVec::range(double &min, double &max) const {
  double my_min = 0.0, my_max = 0.0, gmin = 0.0, gmax = 0.0;
  PetscErrorCode ierr;
  assert(m_v != NULL);

  ierr = VecMin(m_v, NULL, &my_min);
  PISM_PETSC_CHK(ierr, "VecMin");
  ierr = VecMax(m_v, NULL, &my_max);
  PISM_PETSC_CHK(ierr, "VecMax");

  if (m_has_ghosts) {
    // needs a reduce operation; use GlobalMax;
    gmin = GlobalMin(m_grid->com, my_min);
    gmax = GlobalMax(m_grid->com, my_max);
    min = gmin;
    max = gmax;
  } else {
    min = my_min;
    max = my_max;
  }
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

//! Computes the norm of an IceModelVec by calling PETSc VecNorm.
/*!
See comment for range(); because local Vecs are VECSEQ, needs a reduce operation.
See src/trypetsc/localVecMax.c.

@note This method works for all IceModelVecs, including ones with
dof > 1. You might want to use norm_all() for IceModelVec2Stag,
though.
 */
void IceModelVec::norm(int n, double &out) const {
  PetscErrorCode ierr;
  double my_norm, gnorm;
  assert(m_v != NULL);

  NormType type = int_to_normtype(n);

  ierr = VecNorm(m_v, type, &my_norm);
  PISM_PETSC_CHK(ierr, "VecNorm");

  if (m_has_ghosts == true) {
    // needs a reduce operation; use GlobalMax if NORM_INFINITY,
    //   otherwise GlobalSum; carefully in NORM_2 case
    if (n == NORM_1_AND_2) {
      throw RuntimeError::formatted("IceModelVec::norm(...): NORM_1_AND_2 not implemented (called as %s.norm(...))",
         m_name.c_str());
    } else if (n == NORM_1) {
      gnorm = GlobalSum(m_grid->com, my_norm);
    } else if (n == NORM_2) {
      my_norm = PetscSqr(my_norm);  // undo sqrt in VecNorm before sum
      gnorm = GlobalSum(m_grid->com, my_norm);
      gnorm = sqrt(gnorm);
    } else if (n == NORM_INFINITY) {
      gnorm = GlobalMax(m_grid->com, my_norm);
    } else {
      throw RuntimeError::formatted("IceModelVec::norm(...): unknown norm type (called as %s.norm(...))",
                                    m_name.c_str());
    }
    out = gnorm;
  } else {
    out = my_norm;
  }
}

//! Result: v <- sqrt(v), elementwise.  Calls VecSqrt(v).
/*!
Name avoids clash with sqrt() in math.h.
 */
void IceModelVec::squareroot() {
  PetscErrorCode ierr;
  assert(m_v != NULL);

  ierr = VecSqrtAbs(m_v);
  PISM_PETSC_CHK(ierr, "VecSqrtAbs");
}


//! Result: v <- v + alpha * x. Calls VecAXPY.
void IceModelVec::add(double alpha, const IceModelVec &x) {
  assert(m_v != NULL && x.m_v != NULL);

  checkCompatibility("add", x);

  PetscErrorCode ierr = VecAXPY(m_v, alpha, x.m_v);
  PISM_PETSC_CHK(ierr, "VecAXPY");

  inc_state_counter();          // mark as modified
}

//! Result: v[j] <- v[j] + alpha for all j. Calls VecShift.
void IceModelVec::shift(double alpha) {
  assert(m_v != NULL);

  PetscErrorCode ierr = VecShift(m_v, alpha);
  PISM_PETSC_CHK(ierr, "VecShift");

  inc_state_counter();          // mark as modified
}

//! Result: v <- v * alpha. Calls VecScale.
void IceModelVec::scale(double alpha) {
  assert(m_v != NULL);

  PetscErrorCode ierr = VecScale(m_v, alpha);
  PISM_PETSC_CHK(ierr, "VecScale");

  inc_state_counter();          // mark as modified
}

//! Copies v to a global vector 'destination'. Ghost points are discarded.
/*! This is potentially dangerous: make sure that `destination` has the same
    dimensions as the current IceModelVec.

    DMLocalToGlobalBegin/End is broken in PETSc 3.5, so we roll our
    own.
 */
void  IceModelVec::copy_to_vec(PISMDM::Ptr destination_da, Vec destination) const {
  assert(m_v != NULL);

  // m_dof > 1 for vector, staggered grid 2D fields, etc. In this case
  // m_n_levels == 1. For 3D fields, m_dof == 1 (all 3D fields are
  // scalar) and m_n_levels corresponds to dof of the underlying PETSc
  // DM object. So we want the bigger of the two numbers here.
  unsigned int N = std::max(m_dof, m_n_levels);

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
    DMGlobalToLocalBegin(*m_da, source, INSERT_VALUES, m_v);
    DMGlobalToLocalEnd(*m_da, source, INSERT_VALUES, m_v);
  } else {
    ierr = VecCopy(source, m_v);
    PISM_PETSC_CHK(ierr, "VecCopy");
  }

  inc_state_counter();          // mark as modified
}


void IceModelVec::get_dof(PISMDM::Ptr da_result, Vec result,
                                    unsigned int start, unsigned int count) const {
  PetscErrorCode ierr;
  void *tmp_res = NULL, *tmp_v = NULL;

  if (start >= m_dof) {
    throw RuntimeError::formatted("invalid argument (start); got %d", start);
  }

  DMDAVecGetArrayDOF(*da_result, result, &tmp_res);
  double ***result_a = static_cast<double***>(tmp_res);

  DMDAVecGetArrayDOF(*m_da, m_v, &tmp_v);
  double ***source_a = static_cast<double***>(tmp_v);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    ierr = PetscMemcpy(result_a[i][j], &source_a[i][j][start],
                       count*sizeof(PetscScalar));
    PISM_PETSC_CHK(ierr, "PetscMemcpy");
  }

  DMDAVecRestoreArray(*da_result, result, &tmp_res);
  DMDAVecRestoreArrayDOF(*m_da, m_v, &tmp_v);
}

void IceModelVec::set_dof(PISMDM::Ptr da_source, Vec source,
                                    unsigned int start, unsigned int count) {
  PetscErrorCode ierr;
  void *tmp_src = NULL, *tmp_v = NULL;

  if (start >= m_dof) {
    throw RuntimeError::formatted("invalid argument (start); got %d", start);
  }

  DMDAVecGetArrayDOF(*da_source, source, &tmp_src);
  double ***source_a = static_cast<double***>(tmp_src);

  DMDAVecGetArrayDOF(*m_da, m_v, &tmp_v);
  double ***result_a = static_cast<double***>(tmp_v);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();
    ierr = PetscMemcpy(&result_a[i][j][start], source_a[i][j],
                       count*sizeof(PetscScalar));
    PISM_PETSC_CHK(ierr, "PetscMemcpy");
  }

  DMDAVecRestoreArray(*da_source, source, &tmp_src);
  DMDAVecRestoreArrayDOF(*m_da, m_v, &tmp_v);

  inc_state_counter();          // mark as modified
}

//! Result: destination <- v.  Leaves metadata alone but copies values in Vec.  Uses VecCopy.
void  IceModelVec::copy_to(IceModelVec &destination) const {
  PetscErrorCode ierr;
  assert(m_v != NULL && destination.m_v != NULL);

  checkCompatibility("copy_to", destination);

  ierr = VecCopy(m_v, destination.m_v);
  PISM_PETSC_CHK(ierr, "VecCopy");
  destination.inc_state_counter();          // mark as modified
}

//! Result: v <- source.  Leaves metadata alone but copies values in Vec.  Uses VecCopy.
void  IceModelVec::copy_from(const IceModelVec &source) {
  source.copy_to(*this);
}

//! @brief Get the stencil width of the current IceModelVec. Returns 0
//! if ghosts are not available.
unsigned int IceModelVec::get_stencil_width() const {
  if (m_has_ghosts) {
    return m_da_stencil_width;
  } else {
    return 0;
  }
}

Vec IceModelVec::get_vec() {
  return m_v;
}

PISMDM::Ptr IceModelVec::get_dm() const {
  return m_da;
}

//! Sets the variable name to `name` and resets metadata.
void  IceModelVec::set_name(const std::string &new_name, int N) {
  reset_attrs(N);

  if (N == 0) {
    m_name = new_name;
  }

  metadata(N).set_name(new_name);
}

std::string IceModelVec::name() const {
  return m_name;
}

//! Sets the variable's various names without changing any other metadata
void IceModelVec::rename(const std::string &short_name, const std::string &long_name,
                                   const std::string &standard_name, int N) {

  if (short_name.empty() == false) {
    if (N == 0) {
      m_name = short_name;
    }
    metadata(N).set_name(short_name);
  }

  if (long_name.empty() == false) {
    metadata(N).set_string("long_name", long_name);
  }

  if (!standard_name.empty()) {
    metadata(N).set_string("standard_name", standard_name);
  }
}

//! Sets the glaciological units of an IceModelVec.
/*!
This affects NCVariable::report_range() and IceModelVec::write().  In write(),
if IceModelVec::write_in_glaciological_units == true, then that variable is written
with a conversion to the glaciological units set here.
 */
void  IceModelVec::set_glaciological_units(const std::string &my_units) {

  for (unsigned int j = 0; j < m_dof; ++j) {
    metadata(j).set_glaciological_units(my_units);
  }
}

//! Resets most IceModelVec attributes.
void IceModelVec::reset_attrs(unsigned int N) {

  write_in_glaciological_units = false;
  m_report_range                 = true;

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

  If my_units != "", this also resets glaciological_units, so that they match
  internal units.
 */
void IceModelVec::set_attrs(const std::string &my_pism_intent,
                                      const std::string &my_long_name,
                                      const std::string &my_units,
                                      const std::string &my_standard_name,
                                      int N) {

  metadata(N).set_string("long_name", my_long_name);

  metadata(N).set_units(my_units);

  metadata(N).set_string("pism_intent", my_pism_intent);

  metadata(N).set_string("standard_name", my_standard_name);
}


//! Gets an IceModelVec from a file `nc`, interpolating onto the current grid.
/*! Stops if the variable was not found and `critical` == true.
 */
void IceModelVec::regrid_impl(const PIO &nc, RegriddingFlag flag,
                                        double default_value) {
  PetscErrorCode ierr;
  Vec tmp;

  if (getVerbosityLevel() > 3) {
    ierr = PetscPrintf(m_grid->com, "  Regridding %s...\n", m_name.c_str());
    PISM_PETSC_CHK(ierr, "PetscPrintf");
  }

  if (m_dof != 1) {
    throw RuntimeError("This method (IceModelVec::regrid_impl) only supports IceModelVecs with dof == 1.");
  }

  if (m_has_ghosts) {
    DMGetGlobalVector(*m_da, &tmp);

    metadata(0).regrid(nc, flag, m_report_range, default_value, tmp);

    DMGlobalToLocalBegin(*m_da, tmp, INSERT_VALUES, m_v);
    DMGlobalToLocalEnd(*m_da, tmp, INSERT_VALUES, m_v);

    DMRestoreGlobalVector(*m_da, &tmp);
  } else {
    metadata(0).regrid(nc, flag, m_report_range, default_value, m_v);
  }
}

//! Reads appropriate NetCDF variable(s) into an IceModelVec.
void IceModelVec::read_impl(const PIO &nc, const unsigned int time) {
  PetscErrorCode ierr;
  Vec tmp;

  if (getVerbosityLevel() > 3) {
    ierr = PetscPrintf(m_grid->com, "  Reading %s...\n", m_name.c_str());
    PISM_PETSC_CHK(ierr, "PetscPrintf");
  }

  if (m_dof != 1) {
    throw RuntimeError("This method (IceModelVec::read_impl) only supports IceModelVecs with dof == 1.");
  }

  if (m_has_ghosts) {
    DMGetGlobalVector(*m_da, &tmp);

    metadata(0).read(nc, time, tmp);

    DMGlobalToLocalBegin(*m_da, tmp, INSERT_VALUES, m_v);
    DMGlobalToLocalEnd(*m_da, tmp, INSERT_VALUES, m_v);

    DMRestoreGlobalVector(*m_da, &tmp);
  } else {
    metadata(0).read(nc, time, m_v);
  }
}

//! \brief Define variables corresponding to an IceModelVec in a file opened using `nc`.
void IceModelVec::define(const PIO &nc, IO_Type output_datatype) const {

  for (unsigned int j = 0; j < m_dof; ++j) {
    metadata(j).define(nc, output_datatype, write_in_glaciological_units);
  }
}

//! \brief Read attributes from the corresponding variable in `nc`.
/*! Note that unlike read() and regrid(), this method does not use the standard
  name to find the variable to read attributes from.
 */
void IceModelVec::read_attributes(const std::string &filename, int N) {
  PIO nc(*m_grid, "netcdf3");     // OK to use netcdf3

  nc.open(filename, PISM_READONLY);
  nc.read_attributes(metadata(N).get_name(), metadata(N));
  nc.close();
}


//! @brief Returns a reference to the NCSpatialVariable object
//! containing metadata for the compoment N.
NCSpatialVariable& IceModelVec::metadata(unsigned int N) {
  assert(N < m_dof);
  return m_metadata[N];
}

const NCSpatialVariable& IceModelVec::metadata(unsigned int N) const {
  assert(N < m_dof);
  return m_metadata[N];
}

//! Writes an IceModelVec to a NetCDF file.
void IceModelVec::write_impl(const PIO &nc, IO_Type nctype) const {
  PetscErrorCode ierr;
  Vec tmp;

  if (getVerbosityLevel() > 3) {
    ierr = PetscPrintf(m_grid->com, "  Writing %s...\n", m_name.c_str());
    PISM_PETSC_CHK(ierr, "PetscPrintf");
  }

  if (m_dof != 1) {
    throw RuntimeError("This method (IceModelVec::write_impl) only supports IceModelVecs with dof == 1");
  }

  if (m_has_ghosts) {
    DMGetGlobalVector(*m_da, &tmp);

    this->copy_to_vec(m_da, tmp);

    metadata(0).write(nc, nctype, write_in_glaciological_units, tmp);

    DMRestoreGlobalVector(*m_da, &tmp);
  } else {
    metadata(0).write(nc, nctype, write_in_glaciological_units, m_v);
  }
}

//! Dumps a variable to a file, overwriting this file's contents (for debugging).
void IceModelVec::dump(const char filename[]) const {
  PIO nc(*m_grid, m_grid->config.get_string("output_format"));

  nc.open(filename, PISM_READWRITE_CLOBBER);
  nc.def_time(m_grid->config.get_string("time_dimension_name"),
              m_grid->time->calendar(),
              m_grid->time->units_string());
  nc.append_time(m_grid->config.get_string("time_dimension_name"),
                 m_grid->time->current());

  write(nc, PISM_DOUBLE);

  nc.close();
}

//! Checks if two IceModelVecs have compatible sizes, dimensions and numbers of degrees of freedom.
void IceModelVec::checkCompatibility(const char* func, const IceModelVec &other) const {
  PetscErrorCode ierr;
  PetscInt X_size, Y_size;

  if (m_dof != other.m_dof) {
    throw RuntimeError::formatted("IceModelVec::%s(...): operands have different numbers of degrees of freedom",
                                  func);
  }

  ierr = VecGetSize(m_v, &X_size);
  PISM_PETSC_CHK(ierr, "VecGetSize");
  ierr = VecGetSize(other.m_v, &Y_size);
  PISM_PETSC_CHK(ierr, "VecGetSize");
  if (X_size != Y_size) {
    throw RuntimeError::formatted("IceModelVec::%s(...): incompatible Vec sizes (called as %s.%s(%s))",
                                  func, m_name.c_str(), func, other.m_name.c_str());
  }
}

//! Checks if an IceModelVec is allocated and calls DAVecGetArray.
void  IceModelVec::begin_access() const {
  assert(m_v != NULL);

  if (m_access_counter < 0) {
    throw RuntimeError("IceModelVec::begin_access(): m_access_counter < 0");
  }

  if (m_access_counter == 0) {

    if (begin_end_access_use_dof == true) {
      DMDAVecGetArrayDOF(*m_da, m_v, &array);
    } else {
      DMDAVecGetArray(*m_da, m_v, &array);
    }
  }

  m_access_counter++;
}

//! Checks if an IceModelVec is allocated and calls DAVecRestoreArray.
void  IceModelVec::end_access() const {
  assert(m_v != NULL);

  if (array == NULL) {
    throw RuntimeError("IceModelVec::end_access(): a == NULL (looks like begin_acces() was not called)");
  }

  if (m_access_counter < 0) {
    throw RuntimeError("IceModelVec::end_access(): m_access_counter < 0");
  }


  m_access_counter--;
  if (m_access_counter == 0) {
    if (begin_end_access_use_dof == true) {
      DMDAVecRestoreArrayDOF(*m_da, m_v, &array);
    } else {
      DMDAVecRestoreArray(*m_da, m_v, &array);
    }
    array = NULL;
  }
}

//! Updates ghost points.
void  IceModelVec::update_ghosts() {
  if (m_has_ghosts == false) {
    return;
  }

  assert(m_v != NULL);
#if PETSC_VERSION_LT(3,5,0)
  DMDALocalToLocalBegin(*m_da, m_v, INSERT_VALUES, m_v);
  DMDALocalToLocalEnd(*m_da, m_v, INSERT_VALUES, m_v);
#else
  DMLocalToLocalBegin(*m_da, m_v, INSERT_VALUES, m_v);
  DMLocalToLocalEnd(*m_da, m_v, INSERT_VALUES, m_v);
#endif
}

//! Scatters ghost points to IceModelVec destination.
void  IceModelVec::update_ghosts(IceModelVec &destination) const {

  // Make sure it is allocated:
  assert(m_v != NULL);
  // Make sure "destination" has ghosts to update.
  assert(destination.m_has_ghosts == true);

  if (m_has_ghosts == true && destination.m_has_ghosts == true) {
#if PETSC_VERSION_LT(3,5,0)
    DMDALocalToLocalBegin(*m_da, m_v, INSERT_VALUES, destination.m_v);
    DMDALocalToLocalEnd(*m_da, m_v, INSERT_VALUES, destination.m_v);
#else
    DMLocalToLocalBegin(*m_da, m_v, INSERT_VALUES, destination.m_v);
    DMLocalToLocalEnd(*m_da, m_v, INSERT_VALUES, destination.m_v);
#endif
    return;
  }

  if (m_has_ghosts == false && destination.m_has_ghosts == true) {
    DMGlobalToLocalBegin(*destination.m_da, m_v, INSERT_VALUES, destination.m_v);
    DMGlobalToLocalEnd(*destination.m_da, m_v, INSERT_VALUES, destination.m_v);
    return;
  }

  destination.inc_state_counter();          // mark as modified
}

//! Result: v[j] <- c for all j.
void  IceModelVec::set(const double c) {
  PetscErrorCode ierr;
  assert(m_v != NULL);
  ierr = VecSet(m_v,c);
  PISM_PETSC_CHK(ierr, "VecSet");
  inc_state_counter();          // mark as modified
}

void IceModelVec::check_array_indices(int i, int j, unsigned int k) const {
  double ghost_width = 0;
  if (m_has_ghosts) {
    ghost_width = m_da_stencil_width;
  }
  // m_dof > 1 for vector, staggered grid 2D fields, etc. In this case
  // m_n_levels == 1. For 3D fields, m_dof == 1 (all 3D fields are
  // scalar) and m_n_levels corresponds to dof of the underlying PETSc
  // DM object. So we want the bigger of the two numbers here.
  unsigned int N = std::max(m_dof, m_n_levels);

  bool out_of_range = (i < m_grid->xs() - ghost_width) ||
    (i > m_grid->xs() + m_grid->xm() + ghost_width) ||
    (j < m_grid->ys() - ghost_width) ||
    (j > m_grid->ys() + m_grid->ym() + ghost_width) ||
    (k >= N);

  assert(out_of_range == false);

  if (array == NULL) {
    PetscPrintf(m_grid->com,
                "PISM ERROR: IceModelVec::begin_access() was not called (name = '%s')\n",
                m_name.c_str());
  }
  assert(array != NULL);
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
  if (z->get_stencil_width() == 0) {
    // z has no ghosts; we can update everything locally
    // (This covers 4 cases.)
    ghosts = 0;
    scatter = false;
  } else if (x->get_stencil_width() == 0 ||
             y->get_stencil_width() == 0) {
    // z has ghosts, but at least one of x and y does not. we have to scatter
    // ghosts.
    // (This covers 3 cases.)
    ghosts = 0;
    scatter = true;
  } else {
    // all of x, y, z have ghosts
    // (The remaining 8-th case.)
    if (z->get_stencil_width() <= x->get_stencil_width() &&
        z->get_stencil_width() <= y->get_stencil_width()) {
      // x and y have enough ghosts to update ghosts of z locally
      ghosts = z->get_stencil_width();
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
void IceModelVec::norm_all(int n, std::vector<double> &result) const {
  PetscErrorCode ierr;

  assert(n == NORM_1 || n == NORM_2 || n == NORM_INFINITY);

  std::vector<double> norm_result(m_dof);
  result.resize(m_dof);

  NormType type = this->int_to_normtype(n);

  ierr = VecStrideNormAll(m_v, type, &norm_result[0]);
  PISM_PETSC_CHK(ierr, "VecStrideNormAll");

  if (m_has_ghosts) {
    // needs a reduce operation; use GlobalMax() if NORM_INFINITY,
    //   otherwise GlobalSum; carefully in NORM_2 case
    if (n == NORM_1_AND_2) {
      throw RuntimeError::formatted("IceModelVec::norm_all(...): NORM_1_AND_2 not implemented (called as %s.norm_all(...))",
                                    m_name.c_str());
    } else if (n == NORM_1) {

      for (unsigned int k = 0; k < m_dof; ++k) {
        result[k] = GlobalSum(m_grid->com, norm_result[k]);
      }

    } else if (n == NORM_2) {

      for (unsigned int k = 0; k < m_dof; ++k) {
        norm_result[k] = PetscSqr(norm_result[k]);  // undo sqrt in VecNorm before sum
        result[k] = GlobalSum(m_grid->com, norm_result[k]);
        result[k] = sqrt(result[k]);
      }

    } else if (n == NORM_INFINITY) {
      for (unsigned int k = 0; k < m_dof; ++k) {
        result[k] = GlobalMax(m_grid->com, norm_result[k]);
      }
    } else {
      throw RuntimeError::formatted("IceModelVec::norm_all(...): unknown norm type (called as %s.norm_all(...))",
                                    m_name.c_str());
    }
  } else {

    for (unsigned int k = 0; k < m_dof; ++k) {
      result[k] = norm_result[k];
    }

  }
}

void IceModelVec::write(const std::string &filename, IO_Type nctype) const {

  PIO nc(*m_grid, m_grid->config.get_string("output_format"));

  // We expect the file to be present and ready to write into.
  nc.open(filename, PISM_READWRITE);

  this->write(nc, nctype);

  nc.close();
}

void IceModelVec::read(const std::string &filename, unsigned int time) {

  PIO nc(*m_grid, "guess_mode");

  nc.open(filename, PISM_READONLY);

  this->read(nc, time);

  nc.close();
}

void IceModelVec::regrid(const std::string &filename, RegriddingFlag flag,
                                   double default_value) {

  PIO nc(*m_grid, "guess_mode");

  nc.open(filename, PISM_READONLY);

  this->regrid(nc, flag, default_value);

  nc.close();
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
  this->regrid_impl(nc, flag, default_value);
  inc_state_counter();          // mark as modified
}

void IceModelVec::read(const PIO &nc, const unsigned int time) {
  this->read_impl(nc, time);
  inc_state_counter();          // mark as modified
}

void IceModelVec::write(const PIO &nc, IO_Type nctype) const {
  write_impl(nc, nctype);
}

IceModelVec::AccessList::AccessList() {
  // empty
}

IceModelVec::AccessList::~AccessList() {
  while (m_vecs.empty() == false) {
    m_vecs.back()->end_access();
    m_vecs.pop_back();
  }
}

IceModelVec::AccessList::AccessList(const IceModelVec &vec) {
  add(vec);
}

void IceModelVec::AccessList::add(const IceModelVec &vec) {
  vec.begin_access();
  m_vecs.push_back(&vec);
}

void convert_vec(Vec v, Unit from, Unit to) {
  PetscErrorCode ierr;

  UnitConverter c(from, to);

  PetscInt data_size = 0;
  ierr = VecGetLocalSize(v, &data_size);
  PISM_PETSC_CHK(ierr, "VecGetLocalSize");

  double *data = NULL;
  ierr = VecGetArray(v, &data);
  PISM_PETSC_CHK(ierr, "VecGetArray");

  c.convert_doubles(data, data_size);

  ierr = VecRestoreArray(v, &data);
  PISM_PETSC_CHK(ierr, "VecRestoreArray");
}

} // end of namespace pism
