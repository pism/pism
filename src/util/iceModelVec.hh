// Copyright (C) 2008--2022 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#ifndef __IceModelVec_hh
#define __IceModelVec_hh

#include <initializer_list>
#include <memory>               // shared_ptr, dynamic_pointer_cast
#include <cstdint>              // uint64_t
#include <array>

#include "stencils.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/error_handling.hh" // RuntimeError

namespace pism {

enum InterpolationType : int;
class IceGrid;
class File;
class SpatialVariableMetadata;

namespace petsc {
class DM;
class Vec;
class Viewer;
} // end of namespace petsc

//! What "kind" of a vector to create: with or without ghosts.
enum IceModelVecKind {WITHOUT_GHOSTS=0, WITH_GHOSTS=1};

class PetscAccessible {
public:
  virtual ~PetscAccessible() = default;
  virtual void begin_access() const = 0;
  virtual void end_access() const = 0;
};

//! Makes sure that we call begin_access() and end_access() for all accessed IceModelVecs.
class AccessList {
public:
  AccessList();
  AccessList(std::initializer_list<const PetscAccessible *> vecs);
  AccessList(const PetscAccessible &v);
  ~AccessList();
  void add(const PetscAccessible &v);
  void add(const std::vector<const PetscAccessible*> &vecs);
private:
  std::vector<const PetscAccessible*> m_vecs;
};

/*!
 * Interpolation helper. Does not check if points needed for interpolation are within the current
 * processor's sub-domain.
 */
template<class F, typename T>
T interpolate(const F &field, double x, double y) {
  auto grid = field.grid();

  int i_left = 0, i_right = 0, j_bottom = 0, j_top = 0;
  grid->compute_point_neighbors(x, y, i_left, i_right, j_bottom, j_top);

  auto w = grid->compute_interp_weights(x, y);

  return (w[0] * field(i_left,  j_bottom) +
          w[1] * field(i_right, j_bottom) +
          w[2] * field(i_right, j_top) +
          w[3] * field(i_left,  j_top));
}

//! \brief Abstract class for reading, writing, allocating, and accessing a
//! DA-based PETSc Vec (2D and 3D fields) from within IceModel.
/*!
  @anchor icemodelvec_use

  This class represents 2D and 3D fields in PISM. Its methods common to all
  the derived classes can be split (roughly) into six kinds:

  - memory allocation (create)
  - point-wise access (begin_access(), end_access())
  - arithmetic (range(), norm(), add(), shift(), scale(), set(), ...)
  - setting or reading metadata (set_attrs(), metadata())
  - file input/output (read, write, regrid)
  - tracking whether a field was updated (get_state_counter(), inc_state_counter())

  ## Memory allocation

  \code
  array::Scalar var(grid, "var_name", WITH_GHOSTS);
  // var is ready to use
  \endcode

  ("WITH_GHOSTS" means "can be used in computations using map-plane neighbors
  of grid points.)

  It is usually a good idea to set variable metadata right after creating it.
  The method set_attrs() is used throughout PISM to set commonly used
  attributes.

  ## Point-wise access

  PETSc performs some pointer arithmetic magic to allow convenient indexing of
  grid point values. Because of this one needs to surround the code using row,
  column or level indexes with begin_access() and end_access() calls:

  \code
  double foo;
  int i = 0, j = 0;
  array::Scalar var;
  // assume that var was allocated
  ierr = var.begin_access(); CHKERRQ(ierr);
  foo = var(i,j) * 2;
  ierr = var.end_access(); CHKERRQ(ierr);
  \endcode

  Please see [this page](@ref computational_grid) for a discussion of the
  organization of PISM's computational grid and examples of for-loops you will
  probably put between begin_access() and end_access().

  To ensure that ghost values are up to date add the following call
  before the code using ghosts:

  \code
  ierr = var.update_ghosts(); CHKERRQ(ierr);
  \endcode

  ## Reading and writing variables

  PISM can read variables either from files with data on a grid matching the
  current grid (read()) or, using bilinear interpolation, from files
  containing data on a different (but compatible) grid (regrid()).

  To write a field to a "prepared" NetCDF file, use write(). (A file is prepared
  if it contains all the necessary dimensions, coordinate variables and global
  metadata.)

  If you need to "prepare" a file, do:
  \code
  File file(grid.com, PISM_NETCDF3);
  io::prepare_for_output(file, *grid.ctx());
  \endcode

  A note about NetCDF write performance: due to limitations of the NetCDF
  (classic, version 3) format, it is significantly faster to
  \code
  for (all variables)
  var.define(...);

  for (all variables)
  var.write(...);
  \endcode

  as opposed to

  \code
  for (all variables) {
  var.define(...);
  var.write(...);
  }
  \endcode

  IceModelVec::define() is here so that we can use the first approach.

  ## Tracking if a field changed

  It is possible to track if a certain field changed with the help of
  get_state_counter() and inc_state_counter() methods.

  For example, PISM's SIA code re-computes the smoothed bed only if the bed
  deformation code updated it:

  \code
  if (bed->get_state_counter() > bed_state_counter) {
  ierr = bed_smoother->preprocess_bed(...); CHKERRQ(ierr);
  bed_state_counter = bed->get_state_counter();
  }
  \endcode

  The state counter is **not** updated automatically. For the code snippet above
  to work, a bed deformation model has to call inc_state_counter() after an
  update.
*/
class IceModelVec : public PetscAccessible {
public:
  virtual ~IceModelVec();

  typedef std::shared_ptr<IceModelVec> Ptr;
  typedef std::shared_ptr<const IceModelVec> ConstPtr;

  //! `dynamic_pointer_cast` wrapper that checks if the cast succeeded.
  template<class T>
  static typename T::Ptr cast(IceModelVec::Ptr input) {
    auto result = std::dynamic_pointer_cast<T,IceModelVec>(input);

    if (not result) {
      throw RuntimeError(PISM_ERROR_LOCATION, "dynamic cast failure");
    }

    return result;
  }

  IceGrid::ConstPtr grid() const;
  unsigned int ndims() const;
  std::vector<int> shape() const;
  //! @brief Returns the number of degrees of freedom per grid point.
  unsigned int ndof() const;
  unsigned int stencil_width() const;
  std::vector<double> levels() const;

  std::array<double,2> range() const;
  std::vector<double> norm(int n) const;

  void add(double alpha, const IceModelVec &x);
  void shift(double alpha);
  void scale(double alpha);

  petsc::Vec& vec() const;
  std::shared_ptr<petsc::DM> dm() const;

  void set_name(const std::string &name);
  const std::string& get_name() const;

  void set_attrs(const std::string &pism_intent,
                 const std::string &long_name,
                 const std::string &units,
                 const std::string &glaciological_units,
                 const std::string &standard_name,
                 unsigned int component);

  void define(const File &file, IO_Type default_type = PISM_DOUBLE) const;

  void read(const std::string &filename, unsigned int time);
  void read(const File &file, unsigned int time);

  void write(const std::string &filename) const;
  void write(const File &file) const;

  void regrid(const std::string &filename, RegriddingFlag flag,
              double default_value = 0.0);
  void regrid(const File &file, RegriddingFlag flag,
              double default_value = 0.0);

  virtual void begin_access() const;
  virtual void end_access() const;
  void update_ghosts();

  std::shared_ptr<petsc::Vec> allocate_proc0_copy() const;
  void put_on_proc0(petsc::Vec &onp0) const;
  void get_from_proc0(petsc::Vec &onp0);

  void set(double c);

  SpatialVariableMetadata& metadata(unsigned int N = 0);

  const SpatialVariableMetadata& metadata(unsigned int N = 0) const;

  int state_counter() const;
  void inc_state_counter();
  void set_time_independent(bool flag);
  void set_interpolation_type(InterpolationType type);

  void view(std::vector<std::shared_ptr<petsc::Viewer> > viewers) const;

protected:
  IceModelVec(IceGrid::ConstPtr grid,
              const std::string &name,
              IceModelVecKind ghostedp,
              size_t dof,
              size_t stencil_width,
              const std::vector<double> &zlevels);
  struct Impl;
  Impl *m_impl;

  // will be cast to double**, double***, or Vector2** in derived classes
  // This is not hidden in m_impl to make it possible to inline operator()
  mutable void *m_array;

  void set_begin_access_use_dof(bool flag);

  void read_impl(const File &file, unsigned int time);
  void regrid_impl(const File &file, RegriddingFlag flag,
                   double default_value = 0.0);
  void write_impl(const File &file) const;

  void checkCompatibility(const char *function, const IceModelVec &other) const;

  //! @brief Check array indices and warn if they are out of range.
  void check_array_indices(int i, int j, unsigned int k) const;

  void copy_to_vec(std::shared_ptr<petsc::DM> destination_da, petsc::Vec &destination) const;
  void get_dof(std::shared_ptr<petsc::DM> da_result, petsc::Vec &result, unsigned int start,
               unsigned int count=1) const;
  void set_dof(std::shared_ptr<petsc::DM> da_source, petsc::Vec &source, unsigned int start,
               unsigned int count=1);
private:
  size_t size() const;
  // disable copy constructor and the assignment operator:
  IceModelVec(const IceModelVec &other);
  IceModelVec& operator=(const IceModelVec&);
public:
  //! Dump an IceModelVec to a file. *This is for debugging only.*
  //! Uses const char[] to make it easier to call it from gdb.
  void dump(const char filename[]) const;

  uint64_t fletcher64_serial() const;
  uint64_t fletcher64() const;
  std::string checksum(bool serial) const;
  void print_checksum(const char *prefix = "", bool serial=false) const;

  typedef pism::AccessList AccessList;
protected:
  void put_on_proc0(petsc::Vec &parallel, petsc::Vec &onp0) const;
  void get_from_proc0(petsc::Vec &onp0, petsc::Vec &parallel) const;
};

/**
 * Convert a PETSc Vec from the units in `from` into units in `to` (in place).
 *
 * @param v data to convert
 * @param system unit system
 * @param spec1 source unit specification string
 * @param spec2 destination unit specification string
 */
void convert_vec(petsc::Vec &v, std::shared_ptr<units::System> system,
                 const std::string &spec1, const std::string &spec2);

} // end of namespace pism

#endif /* __IceModelVec_hh */
