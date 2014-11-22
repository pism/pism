// Copyright (C) 2008--2014 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#include <cstring>
#include <cstdlib>
#include <petscdmda.h>

#include "NCVariable.hh"
#include "pism_const.hh"

#include "IceGrid.hh"

namespace pism {

class PIO;
class LocalInterpCtx;

//! What "kind" of a vector to create: with or without ghosts.
enum IceModelVecKind {WITHOUT_GHOSTS=0, WITH_GHOSTS=1};

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

  Creating an IceModelVec... object does not allocate memory for storing it
  (some IceModelVecs serve as "references" and don't have their own storage).
  To complete IceModelVec... creation, use the "create()" method:

  \code
  IceModelVec2S var;
  ierr = var.create(grid, "var_name", WITH_GHOSTS); CHKERRQ(ierr);
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
  IceModelVec2S var;
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
  PIO nc(grid.com, grid.config.get_string("output_format"));

  std::string time_name = config.get_string("time_dimension_name");
  ierr = nc.open(filename, PISM_READWRITE); CHKERRQ(ierr); // append == false
  ierr = nc.def_time(time_name, grid.time->calendar(),
  grid.time->CF_units_string()); CHKERRQ(ierr);
  ierr = nc.append_time(time_name, grid.time->current()); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);
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
class IceModelVec {
public:
  IceModelVec();
  virtual ~IceModelVec();

  virtual bool was_created() const;
  IceGrid* get_grid() const {
    return grid;
  }
  unsigned int get_ndims() const;
  //! \brief Returns the number of degrees of freedom per grid point.
  unsigned int get_ndof() const {
    return m_dof;
  }
  unsigned int get_stencil_width() const;
  int nlevels() const {
    return m_n_levels;
  }
  std::vector<double>  get_levels() const {
    return zlevels;
  }

  virtual void  range(double &min, double &max) const;
  virtual void  norm(int n, double &out) const;
  virtual void  norm_all(int n, std::vector<double> &result) const;
  virtual void  add(double alpha, IceModelVec &x);
  virtual void  squareroot();
  virtual void  shift(double alpha);
  virtual void  scale(double alpha);
  void copy_to_vec(PISMDM::Ptr destination_da, Vec destination) const;
  void copy_from_vec(Vec source);
  virtual void copy_to(IceModelVec &destination) const;
  void copy_from(const IceModelVec &source);
  Vec get_vec();
  PISMDM::Ptr get_dm() const;
  virtual void  set_name(const std::string &name, int component = 0);
  virtual std::string name() const;
  virtual void  set_glaciological_units(const std::string &units);
  virtual void  set_attrs(const std::string &my_pism_intent, const std::string &my_long_name,
                                    const std::string &my_units, const std::string &my_standard_name, int component = 0);
  virtual void  rename(const std::string &short_name, const std::string &long_name,
                                 const std::string &standard_name, int component = 0);
  virtual void  read_attributes(const std::string &filename, int component = 0);
  virtual void  define(const PIO &nc, IO_Type output_datatype) const;

  void read(const std::string &filename, unsigned int time);
  void read(const PIO &nc, unsigned int time);

  void  write(const std::string &filename, IO_Type nctype = PISM_DOUBLE) const;
  void  write(const PIO &nc, IO_Type nctype = PISM_DOUBLE) const;

  void  regrid(const std::string &filename, RegriddingFlag flag,
                         double default_value = 0.0);
  void  regrid(const PIO &nc, RegriddingFlag flag,
                         double default_value = 0.0);

  virtual void  begin_access() const;
  virtual void  end_access() const;
  virtual void  update_ghosts();
  virtual void  update_ghosts(IceModelVec &destination) const;

  void  set(double c);

  NCSpatialVariable& metadata(unsigned int N = 0);

  const NCSpatialVariable& metadata(unsigned int N = 0) const;

  int get_state_counter() const;
  void inc_state_counter();
  void set_time_independent(bool flag);

  bool   m_report_range;                 //!< If true, report range when regridding.
  bool   write_in_glaciological_units;
  //!< \brief If true, data is written to a file in "human-friendly" units.

protected:
  virtual void read_impl(const PIO &nc, unsigned int time);
  virtual void regrid_impl(const PIO &nc, RegriddingFlag flag,
                                     double default_value = 0.0);
  virtual void write_impl(const PIO &nc, IO_Type nctype = PISM_DOUBLE) const;
  std::vector<double> zlevels;
  unsigned int m_n_levels;                 //!< number of vertical levels

  Vec  m_v;                       //!< Internal storage
  std::string m_name;

  //! stores metadata (NetCDF variable attributes)
  std::vector<NCSpatialVariable> m_metadata;

  IceGrid *grid;

  unsigned int m_dof;                     //!< number of "degrees of freedom" per grid point
  unsigned int m_da_stencil_width;      //!< stencil width supported by the DA
  bool m_has_ghosts;            //!< m_has_ghosts == true means "has ghosts"
  PISMDM::Ptr m_da;          //!< distributed mesh manager (DM)

  bool begin_end_access_use_dof;

  //! It is a map, because a temporary IceModelVec can be used to view
  //! different quantities
  mutable std::map<std::string,PetscViewer> map_viewers;

  mutable void *array;  // will be cast to double** or double*** in derived classes

  mutable int m_access_counter;           // used in begin_access() and end_access()
  int m_state_counter;            //!< Internal IceModelVec "revision number"

  virtual void destroy();
  virtual void checkCompatibility(const char *function, const IceModelVec &other) const;

  //! \brief Check the array indices and warn if they are out of range.
  void check_array_indices(int i, int j, unsigned int k) const;
  virtual void reset_attrs(unsigned int N);
  NormType int_to_normtype(int input) const;

  void get_dof(PISMDM::Ptr da_result, Vec result, unsigned int n,
                         unsigned int count=1) const;
  void set_dof(PISMDM::Ptr da_source, Vec source, unsigned int n,
                         unsigned int count=1);
private:
  // disable copy constructor and the assignment operator:
  IceModelVec(const IceModelVec &other);
  IceModelVec& operator=(const IceModelVec&);
public:
  //! Dump an IceModelVec to a file. *This is for debugging only.*
  //! Uses const char[] to make it easier to call it from gdb.
  void dump(const char filename[]) const;

public:

  //! Makes sure that we call begin_access() and end_access() for all accessed IceModelVecs.
  class AccessList {
  public:
    AccessList();
    AccessList(const IceModelVec &v);
    ~AccessList();
    void add(const IceModelVec &v);
  private:
    std::vector<const IceModelVec*> m_vecs;
  };
};

enum Direction {North = 0, East, South, West};

//! \brief Star stencil points (in the map-plane).
template <typename T>
struct planeStar {
  T ij, e, w, n, s;
  void set(T input) {
    ij = e = w = n = s = input;
  }

  //! Get the element corresponding to a given direction.
  //! Use foo.ij to get the value at i,j (center of the star).
  inline T& operator[](Direction direction) {
    switch (direction) {
    default:                    // just to silence the warning
    case North:
      return n;
    case East:
      return e;
    case South:
      return s;
    case West:
      return w;
    }
  }
};

class IceModelVec2S;

/** Class for a 2d DA-based Vec.

    As for the difference between IceModelVec2 and IceModelVec2S, the
    former can store fields with more than 1 "degree of freedom" per grid
    point (such as 2D fields on the "staggered" grid, with the first
    degree of freedom corresponding to the i-offset and second to
    j-offset). */
class IceModelVec2 : public IceModelVec {
public:
  IceModelVec2() : IceModelVec() {}
  virtual void view(int viewer_size) const;
  virtual void view(PetscViewer v1, PetscViewer v2) const;
  // component-wise access:
  virtual void get_component(unsigned int n, IceModelVec2S &result) const;
  virtual void set_component(unsigned int n, IceModelVec2S &source);
  inline double& operator() (int i, int j, int k) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, k);
#endif
    return static_cast<double***>(array)[i][j][k];
  }

  inline const double& operator() (int i, int j, int k) const {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, k);
#endif
    return static_cast<double***>(array)[i][j][k];
  }
  virtual void create(IceGrid &my_grid, const std::string &my_short_name,
                                IceModelVecKind ghostedp, unsigned int stencil_width, int dof);
protected:
  virtual void read_impl(const PIO &nc, const unsigned int time);
  virtual void regrid_impl(const PIO &nc, RegriddingFlag flag,
                                     double default_value = 0.0);
  virtual void write_impl(const PIO &nc, IO_Type nctype = PISM_DOUBLE) const;
};

/** A class for storing and accessing scalar 2D fields.
    IceModelVec2S is just IceModelVec2 with "dof == 1" */
class IceModelVec2S : public IceModelVec2 {
  friend class IceModelVec2V;
  friend class IceModelVec2Stag;
public:
  IceModelVec2S() { begin_end_access_use_dof = false; }
  // does not need a copy constructor, because it does not add any new data members
  using IceModelVec2::create;
  virtual void  create(IceGrid &my_grid, const std::string &my_name,
                                 IceModelVecKind ghostedp, int width = 1);
  void allocate_proc0_copy(Vec &result) const;
  void put_on_proc0(Vec onp0) const;
  void get_from_proc0(Vec onp0);
  virtual void  copy_to(IceModelVec &destination) const;
  void  get_array(double** &a);
  virtual void set_to_magnitude(IceModelVec2S &v_x, IceModelVec2S &v_y);
  virtual void mask_by(IceModelVec2S &M, double fill = 0.0);
  virtual void add(double alpha, IceModelVec &x);
  virtual void add(double alpha, const IceModelVec &x, IceModelVec &result) const;
  virtual void sum(double &result);
  virtual void min(double &result) const;
  virtual void max(double &result) const;
  virtual void absmax(double &result) const;
  virtual double diff_x(int i, int j) const;
  virtual double diff_y(int i, int j) const;
  virtual double diff_x_stagE(int i, int j) const;
  virtual double diff_y_stagE(int i, int j) const;
  virtual double diff_x_stagN(int i, int j) const;
  virtual double diff_y_stagN(int i, int j) const;
  virtual double diff_x_p(int i, int j) const;
  virtual double diff_y_p(int i, int j) const;

  //! Provides access (both read and write) to the internal double array.
  /*!
    Note that i corresponds to the x direction and j to the y.
  */
  inline double& operator() (int i, int j) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, 0);
#endif
    return static_cast<double**>(array)[i][j];
  }

  inline const double& operator()(int i, int j) const {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, 0);
#endif
    return static_cast<double**>(array)[i][j];
  }

  inline planeStar<double> star(int i, int j) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, 0);
    check_array_indices(i+1, j, 0);
    check_array_indices(i-1, j, 0);
    check_array_indices(i, j+1, 0);
    check_array_indices(i, j-1, 0);
#endif
    planeStar<double> result;

    result.ij = operator()(i,j);
    result.e =  operator()(i+1,j);
    result.w =  operator()(i-1,j);
    result.n =  operator()(i,j+1);
    result.s =  operator()(i,j-1);

    return result;
  }
};


//! \brief A simple class "hiding" the fact that the mask is stored as
//! floating-point scalars (instead of integers).
class IceModelVec2Int : public IceModelVec2S {
public:
  inline int as_int(int i, int j) const {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, 0);
#endif
    const double **a = (const double**) array;
    return static_cast<int>(floor(a[i][j] + 0.5));
  }

  inline planeStar<int> int_star(int i, int j) const {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, 0);
    check_array_indices(i+1, j, 0);
    check_array_indices(i-1, j, 0);
    check_array_indices(i, j+1, 0);
    check_array_indices(i, j-1, 0);
#endif

    planeStar<int> result;
    result.ij = as_int(i,j);
    result.e =  as_int(i+1,j);
    result.w =  as_int(i-1,j);
    result.n =  as_int(i,j+1);
    result.s =  as_int(i,j-1);

    return result;
  }
};

//! \brief A class representing a horizontal velocity at a certain grid point.
class Vector2 {
public:
  Vector2() : u(0), v(0) {}
  Vector2(double a, double b) : u(a), v(b) {}

  //! Magnitude squared.
  inline double magnitude_squared() const {
    return u*u + v*v;
  }
  //! Magnitude.
  inline double magnitude() const {
    return sqrt(magnitude_squared());
  }

  inline Vector2& operator=(const Vector2 &other) {
    // NOTE: we don't check for self-assignment because there is no memory
    // (de-)allocation here.
    u = other.u;
    v = other.v;
    return *this;
  }

  inline Vector2& operator+=(const Vector2 &other) {
    u += other.u;
    v += other.v;
    return *this;
  }

  inline Vector2& operator-=(const Vector2 &other) {
    u -= other.u;
    v -= other.v;
    return *this;
  }

  inline Vector2& operator*=(const double &a) {
    u *= a;
    v *= a;
    return *this;
  }

  inline Vector2& operator/=(const double &a) {
    u /= a;
    v /= a;
    return *this;
  }

  //! \brief Adds two vectors.
  inline Vector2 operator+(const Vector2 &other) const {
    return Vector2(u + other.u, v + other.v);
  }

  //! \brief Substracts two vectors.
  inline Vector2 operator-(const Vector2 &other) const {
    return Vector2(u - other.u, v - other.v);
  }

  //! \brief Scales a vector.
  inline Vector2 operator*(const double &a) const {
    return Vector2(u * a, v * a);
  }

  //! \brief Scales a vector.
  inline Vector2 operator/(const double &a) const {
    return Vector2(u / a, v / a);
  }

  double u, v;
};

/** Class for storing and accessing 2D vector fields used in IceModel.
    IceModelVec2V is IceModelVec2 with "dof == 2". (Plus some extra methods, of course.)
*/
class IceModelVec2V : public IceModelVec2 {
public:
  IceModelVec2V();
  ~IceModelVec2V() {}

  using IceModelVec2::create;
  virtual void create(IceGrid &my_grid, const std::string &my_short_name,
                                IceModelVecKind ghostedp, unsigned int stencil_width = 1);
  virtual void copy_to(IceModelVec &destination) const;
  virtual void add(double alpha, IceModelVec &x);
  virtual void add(double alpha, IceModelVec &x, IceModelVec &result) const;

  // I/O:
  virtual void get_array(Vector2 ** &a);
  virtual void magnitude(IceModelVec2S &result) const;
  inline Vector2& operator()(int i, int j) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, 0);
#endif
    return static_cast<Vector2**>(array)[i][j];
  }

  inline const Vector2& operator()(int i, int j) const {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, 0);
#endif
    return static_cast<Vector2**>(array)[i][j];
  }

  inline planeStar<Vector2> star(int i, int j) const {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, 0);
    check_array_indices(i+1, j, 0);
    check_array_indices(i-1, j, 0);
    check_array_indices(i, j+1, 0);
    check_array_indices(i, j-1, 0);
#endif
    planeStar<Vector2> result;

    result.ij = operator()(i,j);
    result.e =  operator()(i+1,j);
    result.w =  operator()(i-1,j);
    result.n =  operator()(i,j+1);
    result.s =  operator()(i,j-1);

    return result;
  }

  // Metadata, etc:
  virtual void set_name(const std::string &name, int component = 0);
  virtual void rename(const std::string &short_name, const std::string &long_name,
                                const std::string &standard_name, int component = 0);
  virtual void rename(const std::string &short_name,
                                const std::vector<std::string> &long_names,
                                const std::string &standard_name);
};

//! \brief A class for storing and accessing internal staggered-grid 2D fields.
//! Uses dof=2 storage. This class is identical to IceModelVec2V, except that
//! components are not called `u` and `v` (to avoid confusion).
class IceModelVec2Stag : public IceModelVec2 {
public:
  IceModelVec2Stag() : IceModelVec2() {
    m_dof = 2;
    begin_end_access_use_dof = true;
  }
  using IceModelVec2::create;
  virtual void create(IceGrid &my_grid, const std::string &my_short_name, IceModelVecKind ghostedp,
                                unsigned int stencil_width = 1);
  virtual void staggered_to_regular(IceModelVec2S &result) const;
  virtual void staggered_to_regular(IceModelVec2V &result) const;
  virtual void absmaxcomponents(double* z) const;

  //! Returns the values at interfaces of the cell i,j using the staggered grid.
  /*! The ij member of the return value is set to 0, since it has no meaning in
    this context.
  */
  inline planeStar<double> star(int i, int j) const {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, 0);
    check_array_indices(i+1, j, 0);
    check_array_indices(i-1, j, 0);
    check_array_indices(i, j+1, 0);
    check_array_indices(i, j-1, 0);
#endif
    planeStar<double> result;

    result.ij = 0.0;             // has no meaning in this context
    result.e =  operator()(i, j, 0);
    result.w =  operator()(i-1, j, 0);
    result.n =  operator()(i, j, 1);
    result.s =  operator()(i, j-1, 1);

    return result;
  }
};

//! \brief A virtual class collecting methods common to ice and bedrock 3D
//! fields.
class IceModelVec3D : public IceModelVec
{
public:
  IceModelVec3D();
  virtual ~IceModelVec3D();
public:

  void  setColumn(int i, int j, double c);
  void  setInternalColumn(int i, int j, double *valsIN);
  void  getInternalColumn(int i, int j, double **valsOUT);
  void  getInternalColumn(int i, int j, const double **valsOUT) const;

  virtual double    getValZ(int i, int j, double z) const;
  virtual bool isLegalLevel(double z) const;

  inline double& operator() (int i, int j, int k) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, k);
#endif
    return static_cast<double***>(array)[i][j][k];
  }

  inline const double& operator() (int i, int j, int k) const {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, k);
#endif
    return static_cast<double***>(array)[i][j][k];
  }

protected:
  virtual void allocate(IceGrid &mygrid, const std::string &my_short_name,
                                  IceModelVecKind ghostedp, const std::vector<double> &levels,
                                  unsigned int stencil_width = 1);
};


//! Class for a 3d DA-based Vec for ice scalar quantities.
class IceModelVec3 : public IceModelVec3D {
public:
  IceModelVec3() {}
  virtual ~IceModelVec3() {}

  virtual void create(IceGrid &mygrid, const std::string &my_short_name,
                                IceModelVecKind ghostedp,
                                unsigned int stencil_width = 1);

  // need to call begin_access() before set...(i,j,...) or get...(i,j,...) *and* need call
  // end_access() afterward
  void  getValColumn(int i, int j, unsigned int ks, double *valsOUT) const;
  void  getValColumnQUAD(int i, int j, unsigned int ks, double *valsOUT) const;
  void  getValColumnPL(int i, int j, unsigned int ks, double *valsOUT) const;

  void  setValColumnPL(int i, int j, std::vector<double> &values_fine);

  void  getPlaneStarZ(int i, int j, double z,
                                planeStar<double> *star) const;
  void  getPlaneStar_fine(int i, int j, unsigned int k,
                                    planeStar<double> *star) const;
  void  getPlaneStar(int i, int j, unsigned int k,
                               planeStar<double> *star) const;

  void  getHorSlice(Vec &gslice, double z) const; // used in iMmatlab.cc
  void  getHorSlice(IceModelVec2S &gslice, double z) const;
  void  getSurfaceValues(IceModelVec2S &gsurf, const IceModelVec2S &myH) const;
};

} // end of namespace pism

#endif /* __IceModelVec_hh */

