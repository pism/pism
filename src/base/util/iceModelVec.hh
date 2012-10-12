// Copyright (C) 2008--2012 Ed Bueler, Constantine Khroulev, and David Maxwell
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

#ifndef __IceModelVec_hh
#define __IceModelVec_hh

#include <cstring>
#include <cstdlib>
#include <petscdmda.h>

#include "NCSpatialVariable.hh"
#include "pism_const.hh"

class PIO;
class LocalInterpCtx;

//! \brief Abstract class for reading, writing, allocating, and accessing a
//! DA-based PETSc Vec from within IceModel.
/*!
 This class represents 2D and 3D fields in PISM. Its methods common to all
 the derived classes can be split (roughly) into six kinds:

 \li memory allocation (create)
 \li point-wise access (begin_access, get_array, end_access)
 \li arithmetic (range, norm, add, shift, scale, multiply_by, set)
 \li setting or reading metadata (set_attrs, set_units, etc)
 \li file input/output (read, write, regrid)
 \li tracking whether a field was updated (get_state_counter, inc_state_counter)

 \section imv_allocation Memory allocation

 Creating an IceModelVec... object does not allocate memory for storing it
 (some IceModelVecs serve as "references" and don't have their own storage).
 To complete IceModelVec... creation, use the "create()" method:

 \code
 IceModelVec2S var;
 bool has_ghosts = true;
 ierr = var.create(grid, "var_name", has_ghosts); CHKERRQ(ierr);
 // var is ready to use
 \endcode

 ("Has ghosts" means "can be used in computations using map-plane neighbors
 of grid points.)

 It is usually a good idea to set variable metadata right after creating it.
 The method set_attrs() is used throughout PISM to set commonly used
 attributes.

 \section imv_pointwise Pointwise access

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

 Please see \ref computational_grid "this page" for a discussion of the
 organization of PISM's computational grid and examples of for-loops you will
 probably put between begin_access() and end_access().

 To ensure that ghost values are up to date add the following two lines
 before the code using ghosts:

 \code
 ierr = var.beginGhostComm(); CHKERRQ(ierr);
 ierr = var.endGhostComm(); CHKERRQ(ierr);
 \endcode

 \section imv_io Reading and writing variables

 PISM can read variables either from files with data on a grid matching the
 current grid (read()) or, using bilinear interpolation, from files
 containing data on a different (but compatible) grid (regrid()).

 To write a field to a prepared NetCDF file, use write(). (A file is prepared
 if it contains all the necessary dimensions, coordinate variables and global
 metadata.)

 If you need to "prepare" a file, do:
 \code
 PIO nc(grid.com, grid.rank, grid.config.get_string("output_format"));

 ierr = nc.open(filename, NC_WRITE); CHKERRQ(ierr);
 ierr = nc.def_time(config.get_string("time_dimension_name"),
                    config.get_string("calendar"),
                    grid.time->CF_units()); CHKERRQ(ierr);
 ierr = nc.append_time(grid.time->current()); CHKERRQ(ierr);
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

 \section imv_rev_counter Tracking if a field changed

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

 The state counter is \b not updated automatically. For the code snippet above
 to work, a bed deformation model has to call inc_state_counter() after an
 update.

 */
class IceModelVec {
public:
  IceModelVec();
  IceModelVec(const IceModelVec &other);
  virtual ~IceModelVec();

  virtual bool            was_created();
  virtual IceGrid*        get_grid() { return grid; }
  virtual int             get_ndims();
  //! \brief Returns the number of degrees of freedom per grid point.
  virtual int             get_dof() { return dof; }
  virtual int             get_stencil_width() { return da_stencil_width; }
  virtual int             get_nlevels() { return n_levels; }
  virtual vector<double>  get_levels() { return zlevels; }
  virtual bool            has_ghosts() { return localp; }

  virtual PetscErrorCode  range(PetscReal &min, PetscReal &max);
  virtual PetscErrorCode  norm(NormType n, PetscReal &out);
  virtual PetscErrorCode  add(PetscScalar alpha, IceModelVec &x);
  virtual PetscErrorCode  add(PetscScalar alpha, IceModelVec &x, IceModelVec &result);
  virtual PetscErrorCode  squareroot();
  virtual PetscErrorCode  shift(PetscScalar alpha);
  virtual PetscErrorCode  scale(PetscScalar alpha);
  virtual PetscErrorCode  multiply_by(IceModelVec &x, IceModelVec &result);
  virtual PetscErrorCode  multiply_by(IceModelVec &x);
  virtual PetscErrorCode  copy_to(Vec destination);
  virtual PetscErrorCode  copy_from(Vec source);
  virtual PetscErrorCode  copy_to(IceModelVec &destination);
  virtual PetscErrorCode  copy_from(IceModelVec &source);
  virtual PetscErrorCode  has_nan();
  virtual PetscErrorCode  set_name(string name, int component = 0);
  virtual PetscErrorCode  set_glaciological_units(string units);
  virtual PetscErrorCode  set_attr(string name, string value, int component = 0);
  virtual PetscErrorCode  set_attr(string name, double value, int component = 0);
  virtual PetscErrorCode  set_attr(string name, vector<double> values, int component = 0);
  virtual bool            has_attr(string name, int component = 0);
  virtual string          string_attr(string name, int component = 0);
  virtual double          double_attr(string name, int component = 0);
  virtual vector<double>  array_attr(string name, int component = 0);
  virtual PetscErrorCode  set_attrs(string my_pism_intent, string my_long_name,
				    string my_units, string my_standard_name, int component = 0);
  virtual PetscErrorCode  rename(const string &short_name, const string &long_name,
                                 const string &standard_name, int component = 0);
  virtual PetscErrorCode  set_intent(string pism_intent, int component = 0);
  virtual PetscErrorCode  read_attributes(string filename, int component = 0);
  virtual NCSpatialVariable get_metadata(int N = 0);
  virtual PetscErrorCode  set_metadata(NCSpatialVariable &var, int N);
  virtual bool            is_valid(PetscScalar a, int component = 0);
  virtual PetscErrorCode  define(const PIO &nc, PISM_IO_Type output_datatype);
  virtual PetscErrorCode  write(string filename);
  virtual PetscErrorCode  write(string filename, PISM_IO_Type nctype);
  virtual PetscErrorCode  dump(const char filename[]);
  virtual PetscErrorCode  read(string filename, unsigned int time);
  virtual PetscErrorCode  regrid(string filename, bool critical, int start = 0);
  virtual PetscErrorCode  regrid(string filename, PetscScalar default_value);

  virtual PetscErrorCode  begin_access();
  virtual PetscErrorCode  end_access();
  virtual PetscErrorCode  beginGhostComm();
  virtual PetscErrorCode  endGhostComm();
  virtual PetscErrorCode  beginGhostComm(IceModelVec &destination);
  virtual PetscErrorCode  endGhostComm(IceModelVec &destination);

  virtual PetscErrorCode  set(PetscScalar c);

  virtual int get_state_counter() const;
  virtual void inc_state_counter();

  bool   report_range;                 //!< If true, report range when regridding.
  bool   write_in_glaciological_units, //!< \brief If true, data is written to
  //!< a file in "human-friendly" units.
    time_independent;                  //!< \brief If true, corresponding
                                       //!< NetCDF variables do not depend on
                                       //!< the 't' dimension.
  PISM_IO_Type output_data_type;            //!< Corresponding NetCDF data type.
protected:
  vector<double> zlevels;
  int n_levels;                 //!< number of vertical levels

  bool shallow_copy;            //!< True if this IceModelVec is a shallow copy.
  Vec  v;                       //!< Internal storage
  string name;

  //! \brief NetCDF variable(s) corresponding to this IceModelVec; dof == 1
  //! vectors only have vars[0].
  vector<NCSpatialVariable> vars;

  IceGrid      *grid;
  int          dof,             //!< number of "degrees of freedom" per grid point
    da_stencil_width;           //!< stencil width supported by the DA
  DM           da;
  bool         localp;          //!< localp == true means "has ghosts"

  //! It is a map, because a temporary IceModelVec can be used to view
  //! different quantities, and a pointer because "shallow copies" should have
  //! the acces to the original map
  map<string,PetscViewer> *map_viewers;

  void         *array;  // will be PetscScalar** or PetscScalar*** in derived classes

  int access_counter;		// used in begin_access() and end_access()
  int state_counter;            //!< Internal IceModelVec "revision number"

  virtual PetscErrorCode create_2d_da(DM &result, PetscInt da_dof, PetscInt stencil_width);
  virtual PetscErrorCode destroy();
  virtual PetscErrorCode checkAllocated();
  virtual PetscErrorCode checkHaveArray();
  virtual PetscErrorCode checkCompatibility(const char*, IceModelVec &other);

  //! \brief Check the array indices and warn if they are out of range.
  void check_array_indices(int i, int j);
  virtual PetscErrorCode reset_attrs(int N);
  virtual PetscErrorCode get_interp_context(string filename, LocalInterpCtx* &lic);
};


enum PISM_Direction {North = 0, East, South, West};

//! \brief Star stencil points (in the map-plane).
template <typename T>
struct planeStar {
  T ij, e, w, n, s;
  void set(T input) {
    ij = e = w = n = s = input;
  }

  //! Get the element corresponding to a given direction.
  //! Use foo.ij to get the value at i,j.
  T& operator[](PISM_Direction direction) {
    switch (direction) {
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

//! Class for a 2d DA-based Vec.
class IceModelVec2 : public IceModelVec {
public:
  IceModelVec2() : IceModelVec() {}
  IceModelVec2(const IceModelVec2 &other) : IceModelVec(other) {};
  virtual PetscErrorCode view(PetscInt viewer_size);
  virtual PetscErrorCode view(PetscViewer v1, PetscViewer v2);
  using IceModelVec::write;
  virtual PetscErrorCode write(string filename, PISM_IO_Type nctype);
  virtual PetscErrorCode read(string filename, const unsigned int time);
  virtual PetscErrorCode regrid(string filename, bool critical, int start = 0);
  virtual PetscErrorCode regrid(string filename, PetscScalar default_value);
  // component-wise access:
  virtual PetscErrorCode get_component(int n, IceModelVec2S &result);
  virtual PetscErrorCode set_component(int n, IceModelVec2S &source);
protected:
  virtual PetscErrorCode create(IceGrid &my_grid, string my_short_name, bool has_ghosts,
                                int stencil_width, int dof);
  PetscErrorCode get_component(int n, Vec result);
  PetscErrorCode set_component(int n, Vec source);
};

//! A class for storing and accessing scalar 2D fields.
class IceModelVec2S : public IceModelVec2 {
  friend class IceModelVec2V;
  friend class IceModelVec2Stag;
public:
  IceModelVec2S() {}
  IceModelVec2S(const IceModelVec2S &other) : IceModelVec2(other) {}
  // does not need a copy constructor, because it does not add any new data members
  using IceModelVec2::create;
  virtual PetscErrorCode  create(IceGrid &my_grid, string my_name, bool has_ghosts, int width = 1);
  virtual PetscErrorCode  put_on_proc0(Vec onp0, VecScatter ctx, Vec g2, Vec g2natural);
  virtual PetscErrorCode  get_from_proc0(Vec onp0, VecScatter ctx, Vec g2, Vec g2natural);
  using IceModelVec::copy_to;
  using IceModelVec::copy_from;
  virtual PetscErrorCode  copy_to(IceModelVec &destination);
  virtual PetscErrorCode  copy_from(IceModelVec &source);
  PetscErrorCode  get_array(PetscScalar** &a);
  virtual PetscErrorCode set_to_magnitude(IceModelVec2S &v_x, IceModelVec2S &v_y);
  virtual PetscErrorCode mask_by(IceModelVec2S &M, PetscScalar fill = 0.0);
  virtual PetscErrorCode add(PetscScalar alpha, IceModelVec &x);
  virtual PetscErrorCode add(PetscScalar alpha, IceModelVec &x, IceModelVec &result);
  virtual PetscErrorCode multiply_by(IceModelVec &x, IceModelVec &result);
  virtual PetscErrorCode multiply_by(IceModelVec &x);
  virtual PetscErrorCode sum(PetscScalar &result);
  virtual PetscErrorCode min(PetscScalar &result);
  virtual PetscErrorCode max(PetscScalar &result);
  virtual PetscScalar diff_x(int i, int j);
  virtual PetscScalar diff_y(int i, int j);
  virtual PetscScalar diff_x_stagE(int i, int j);
  virtual PetscScalar diff_y_stagE(int i, int j);
  virtual PetscScalar diff_x_stagN(int i, int j);
  virtual PetscScalar diff_y_stagN(int i, int j);
  virtual PetscScalar diff_x_p(int i, int j);
  virtual PetscScalar diff_y_p(int i, int j);
  virtual PetscErrorCode view_matlab(PetscViewer my_viewer);
  virtual PetscErrorCode has_nan();

  //! Provides access (both read and write) to the internal PetscScalar array.
  /*!
    Note that i corresponds to the x direction and j to the y.
  */
  inline PetscScalar& operator() (int i, int j) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j);
#endif
    return static_cast<PetscScalar**>(array)[i][j];
  }

  inline planeStar<PetscScalar> star(int i, int j) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j);
    check_array_indices(i+1, j);
    check_array_indices(i-1, j);
    check_array_indices(i, j+1);
    check_array_indices(i, j-1);
#endif
    planeStar<PetscScalar> result;

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
  inline int as_int(int i, int j) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j);
#endif
    const PetscScalar **a = (const PetscScalar**) array;
    return static_cast<int>(floor(a[i][j] + 0.5));
  }

  inline planeStar<int> int_star(int i, int j) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j);
    check_array_indices(i+1, j);
    check_array_indices(i-1, j);
    check_array_indices(i, j+1);
    check_array_indices(i, j-1);
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
class PISMVector2 {
public:
  PISMVector2() : u(0), v(0) {}
  PISMVector2(PetscScalar a, PetscScalar b) : u(a), v(b) {}

  //! Magnitude squared.
  inline PetscScalar magnitude_squared() const {
    return u*u + v*v;
  }
  //! Magnitude.
  inline PetscScalar magnitude() const {
    return sqrt(magnitude_squared());
  }

  inline PISMVector2& operator=(const PISMVector2 &other) {
    // NOTE: we don't check for self-assignment because there is no memory
    // (de-)allocation here.
    u = other.u;
    v = other.v;
    return *this;
  }

  inline PISMVector2& operator+=(const PISMVector2 &other) {
    u += other.u;
    v += other.v;
    return *this;
  }

  inline PISMVector2& operator-=(const PISMVector2 &other) {
    u -= other.u;
    v -= other.v;
    return *this;
  }

  inline PISMVector2& operator*=(const PetscScalar &a) {
    u *= a;
    v *= a;
    return *this;
  }

  inline PISMVector2& operator/=(const PetscScalar &a) {
    u /= a;
    v /= a;
    return *this;
  }

  //! \brief Adds two vectors.
  inline PISMVector2 operator+(const PISMVector2 &other) const {
    return PISMVector2(u + other.u, v + other.v);
  }

  //! \brief Substracts two vectors.
  inline PISMVector2 operator-(const PISMVector2 &other) const {
    return PISMVector2(u - other.u, v - other.v);
  }

  //! \brief Scales a vector.
  inline PISMVector2 operator*(const PetscScalar &a) const {
    return PISMVector2(u * a, v * a);
  }

  //! \brief Scales a vector.
  inline PISMVector2 operator/(const PetscScalar &a) const {
    return PISMVector2(u / a, v / a);
  }

  PetscScalar u, v;
};

inline PISMVector2 operator*(const PetscScalar &a, const PISMVector2 &v1) {
  return v1 * a;
}

//! Class for storing and accessing 2D vector fields used in IceModel.
class IceModelVec2V : public IceModelVec2 {
public:
  IceModelVec2V();
  IceModelVec2V(const IceModelVec2V &other) : IceModelVec2(other) {}

  ~IceModelVec2V() {}

  using IceModelVec2::create;
  virtual PetscErrorCode create(IceGrid &my_grid, string my_short_name,
				bool has_ghosts, int stencil_width = 1);
  using IceModelVec::copy_to;
  using IceModelVec::copy_from;
  virtual PetscErrorCode copy_to(IceModelVec &destination);
  virtual PetscErrorCode copy_from(IceModelVec &source);
  virtual PetscErrorCode add(PetscScalar alpha, IceModelVec &x);
  virtual PetscErrorCode add(PetscScalar alpha, IceModelVec &x, IceModelVec &result);

  // I/O:
  using IceModelVec2::write;
  virtual PetscErrorCode get_array(PISMVector2 ** &a);
  virtual PetscErrorCode magnitude(IceModelVec2S &result);
  inline PISMVector2& operator()(int i, int j) {
#if (PISM_DEBUG == 1)
    check_array_indices(i, j);
#endif
    return static_cast<PISMVector2**>(array)[i][j];
  }

  inline planeStar<PISMVector2> star(int i, int j) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j);
    check_array_indices(i+1, j);
    check_array_indices(i-1, j);
    check_array_indices(i, j+1);
    check_array_indices(i, j-1);
#endif
    planeStar<PISMVector2> result;

    result.ij = operator()(i,j);
    result.e =  operator()(i+1,j);
    result.w =  operator()(i-1,j);
    result.n =  operator()(i,j+1);
    result.s =  operator()(i,j-1);

    return result;
  }

  // Metadata, etc:
  using IceModelVec2::is_valid;
  virtual bool           is_valid(PetscScalar u, PetscScalar v);
  virtual PetscErrorCode set_name(string name, int component = 0);
  virtual PetscErrorCode rename(const string &short_name, const string &long_name,
                                const string &standard_name, int component = 0);
  virtual PetscErrorCode rename(const string &short_name, const vector<string> &long_names,
                                const string &standard_name);
};

//! \brief A class for storing and accessing internal staggered-grid 2D fields.
//! Uses dof=2 storage. Does \b not support input and output. This class is
//! identical to IceModelVec2V, except that components are not called \c u and
//! \c v (to avoid confusion).
class IceModelVec2Stag : public IceModelVec2 {
public:
  IceModelVec2Stag() { dof = 2; vars.resize(dof); }
  IceModelVec2Stag(const IceModelVec2Stag &other) : IceModelVec2(other) {}
  using IceModelVec2::create;
  virtual PetscErrorCode create(IceGrid &my_grid, string my_name, bool has_ghosts, int width = 1);
  virtual PetscErrorCode get_array(PetscScalar*** &a);
  virtual PetscErrorCode begin_access();
  virtual PetscErrorCode end_access();
  inline PetscScalar& operator() (int i, int j, int k) {
#if (PISM_DEBUG == 1)
    check_array_indices(i, j);
#endif
    return static_cast<PetscScalar***>(array)[i][j][k];
  }
  virtual PetscErrorCode norm_all(NormType n, PetscReal &result0, PetscReal &result1);
  virtual PetscErrorCode staggered_to_regular(IceModelVec2S &result);
  virtual PetscErrorCode staggered_to_regular(IceModelVec2V &result);

  //! Returns the values at interfaces of the cell i,j using the staggered grid.
  /*! The ij member of the return value is set to 0, since it has no meaning in
    this context.
   */
  inline planeStar<PetscScalar> star(int i, int j) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j);
    check_array_indices(i+1, j);
    check_array_indices(i-1, j);
    check_array_indices(i, j+1);
    check_array_indices(i, j-1);
#endif
    planeStar<PetscScalar> result;

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
  IceModelVec3D(const IceModelVec3D &other);
  virtual ~IceModelVec3D();
public:

  virtual PetscErrorCode  begin_access();
  virtual PetscErrorCode  end_access();

  PetscErrorCode  setColumn(PetscInt i, PetscInt j, PetscScalar c);
  PetscErrorCode  setInternalColumn(PetscInt i, PetscInt j, PetscScalar *valsIN);
  PetscErrorCode  getInternalColumn(PetscInt i, PetscInt j, PetscScalar **valsOUT);

  PetscErrorCode  view_sounding(int i, int j, PetscInt viewer_size);
  PetscErrorCode  view_sounding(int i, int j, PetscViewer v);

  // note the IceModelVec3 with this method must be *local* while imv3_source must be *global*
  virtual PetscErrorCode beginGhostCommTransfer(IceModelVec3D &imv3_source);
  virtual PetscErrorCode endGhostCommTransfer(IceModelVec3D &imv3_source);
  virtual PetscScalar    getValZ(PetscInt i, PetscInt j, PetscScalar z);
  virtual PetscErrorCode isLegalLevel(PetscScalar z);
protected:
  virtual PetscErrorCode allocate(IceGrid &mygrid, string my_short_name,
                                  bool has_ghosts, vector<double> levels, int stencil_width = 1);
  virtual PetscErrorCode destroy();
  virtual PetscErrorCode has_nan();

  Vec sounding_buffer;
  map<string,PetscViewer> *sounding_viewers;
};


//! Class for a 3d DA-based Vec for ice scalar quantities.
class IceModelVec3 : public IceModelVec3D {
public:
  IceModelVec3() {}
  IceModelVec3(const IceModelVec3 &other) : IceModelVec3D(other) {}
  virtual ~IceModelVec3() {}

  virtual PetscErrorCode create(IceGrid &mygrid, string my_short_name,
                                bool has_ghosts, int stencil_width = 1);

  // need to call begin_access() before set...(i,j,...) or get...(i,j,...) *and* need call
  // end_access() afterward
  PetscErrorCode  getValColumn(PetscInt i, PetscInt j, PetscInt ks, PetscScalar *valsOUT);
  PetscErrorCode  getValColumnQUAD(PetscInt i, PetscInt j, PetscInt ks, PetscScalar *valsOUT);
  PetscErrorCode  getValColumnPL(PetscInt i, PetscInt j, PetscInt ks, PetscScalar *valsOUT);

  PetscErrorCode  setValColumnPL(PetscInt i, PetscInt j, PetscScalar *valsIN);

  PetscErrorCode  getPlaneStarZ(PetscInt i, PetscInt j, PetscScalar z,
                                planeStar<PetscScalar> *star);
  PetscErrorCode  getPlaneStar_fine(PetscInt i, PetscInt j, PetscInt k,
				    planeStar<PetscScalar> *star);
  PetscErrorCode  getPlaneStar(PetscInt i, PetscInt j, PetscInt k,
			       planeStar<PetscScalar> *star);

  PetscErrorCode  getHorSlice(Vec &gslice, PetscScalar z); // used in iMmatlab.cc
  PetscErrorCode  getHorSlice(IceModelVec2S &gslice, PetscScalar z);
  PetscErrorCode  getSurfaceValues(Vec &gsurf, IceModelVec2S &myH); // used in iMviewers.cc
  PetscErrorCode  getSurfaceValues(IceModelVec2S &gsurf, IceModelVec2S &myH);
  PetscErrorCode  getSurfaceValues(IceModelVec2S &gsurf, PetscScalar **H);
  PetscErrorCode  extend_vertically(int old_Mz, PetscScalar fill_value);
  PetscErrorCode  extend_vertically(int old_Mz, IceModelVec2S &fill_values);
protected:
  virtual PetscErrorCode  extend_vertically_private(int old_Mz);
};

#endif /* __IceModelVec_hh */

