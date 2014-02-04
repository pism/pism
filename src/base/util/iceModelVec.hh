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

class PIO;
class LocalInterpCtx;

//! Named constants for IceModelVec*::create.
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
  ierr = nc.open(filename, PISM_WRITE); CHKERRQ(ierr); // append == false
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

  virtual bool was_created();
  IceGrid* get_grid() { return grid; }
  unsigned int get_ndims();
  //! \brief Returns the number of degrees of freedom per grid point.
  unsigned int get_dof() { return m_dof; }
  unsigned int get_stencil_width() { return m_da_stencil_width; }
  int nlevels() { return m_n_levels; }
  std::vector<double>  get_levels() { return zlevels; }
  bool has_ghosts() { return m_has_ghosts; }

  virtual PetscErrorCode  range(double &min, double &max);
  virtual PetscErrorCode  norm(int n, double &out);
  virtual PetscErrorCode  norm_all(int n, std::vector<double> &result);
  virtual PetscErrorCode  add(double alpha, IceModelVec &x);
  virtual PetscErrorCode  squareroot();
  virtual PetscErrorCode  shift(double alpha);
  virtual PetscErrorCode  scale(double alpha);
  virtual PetscErrorCode  copy_to(Vec destination);
  virtual PetscErrorCode  copy_from(Vec source);
  virtual PetscErrorCode  copy_to(IceModelVec &destination);
  virtual PetscErrorCode  copy_from(IceModelVec &source);
  virtual Vec get_vec();
  virtual PetscErrorCode  has_nan();
  virtual PetscErrorCode  set_name(std::string name, int component = 0);
  virtual std::string name() const;
  virtual PetscErrorCode  set_glaciological_units(std::string units);
  virtual PetscErrorCode  set_attrs(std::string my_pism_intent, std::string my_long_name,
				    std::string my_units, std::string my_standard_name, int component = 0);
  virtual PetscErrorCode  rename(std::string short_name, std::string long_name,
                                 std::string standard_name, int component = 0);
  virtual PetscErrorCode  read_attributes(std::string filename, int component = 0);
  virtual PetscErrorCode  define(const PIO &nc, PISM_IO_Type output_datatype);

  virtual PetscErrorCode  write(std::string filename, PISM_IO_Type nctype = PISM_DOUBLE);
  virtual PetscErrorCode  read(std::string filename, unsigned int time);
  virtual PetscErrorCode  regrid(std::string filename, RegriddingFlag flag,
                                 double default_value = 0.0);

  virtual PetscErrorCode  write(const PIO &nc, PISM_IO_Type nctype = PISM_DOUBLE);
  virtual PetscErrorCode  read(const PIO &nc, unsigned int time);
  virtual PetscErrorCode  regrid(const PIO &nc, RegriddingFlag flag,
                                 double default_value = 0.0);

  virtual PetscErrorCode  begin_access();
  virtual PetscErrorCode  end_access();
  virtual PetscErrorCode  update_ghosts();
  virtual PetscErrorCode  update_ghosts(IceModelVec &destination);

  PetscErrorCode  set(double c);

  NCSpatialVariable& metadata(unsigned int N = 0);

  int get_state_counter() const;
  void inc_state_counter();
  void set_time_independent(bool flag);

  bool   m_report_range;                 //!< If true, report range when regridding.
  bool   write_in_glaciological_units;
  //!< \brief If true, data is written to a file in "human-friendly" units.

protected:
  std::vector<double> zlevels;
  unsigned int m_n_levels;                 //!< number of vertical levels

  Vec  v;                       //!< Internal storage
  std::string m_name;

  //! stores metadata (NetCDF variable attributes)
  std::vector<NCSpatialVariable> m_metadata;

  IceGrid *grid;

  unsigned int m_dof;                     //!< number of "degrees of freedom" per grid point
  unsigned int m_da_stencil_width;      //!< stencil width supported by the DA
  bool m_has_ghosts;            //!< m_has_ghosts == true means "has ghosts"
  DM   m_da;                    //!< DM; this IceModelVec does not own it!

  bool begin_end_access_use_dof;

  //! It is a map, because a temporary IceModelVec can be used to view
  //! different quantities
  std::map<std::string,PetscViewer> map_viewers;

  void *array;  // will be cast to double** or double*** in derived classes

  int access_counter;		// used in begin_access() and end_access()
  int state_counter;            //!< Internal IceModelVec "revision number"

  virtual PetscErrorCode destroy();
  virtual PetscErrorCode checkCompatibility(const char*, IceModelVec &other);

  //! \brief Check the array indices and warn if they are out of range.
  void check_array_indices(int i, int j, unsigned int k);
  virtual PetscErrorCode reset_attrs(unsigned int N);
  NormType int_to_normtype(int input);
private:
  // disable copy constructor and the assignment operator:
  IceModelVec(const IceModelVec &other);
  IceModelVec& operator=(const IceModelVec&);
  PetscErrorCode dump(const char filename[]);
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
  //! Use foo.ij to get the value at i,j (center of the star).
  inline T& operator[](PISM_Direction direction) {
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
  virtual PetscErrorCode view(int viewer_size);
  virtual PetscErrorCode view(PetscViewer v1, PetscViewer v2);
  using IceModelVec::write;
  using IceModelVec::read;
  using IceModelVec::regrid;
  virtual PetscErrorCode write(const PIO &nc, PISM_IO_Type nctype = PISM_DOUBLE);
  virtual PetscErrorCode read(const PIO &nc, const unsigned int time);
  virtual PetscErrorCode regrid(const PIO &nc, RegriddingFlag flag,
                                double default_value = 0.0);
  // component-wise access:
  virtual PetscErrorCode get_component(unsigned int n, IceModelVec2S &result);
  virtual PetscErrorCode set_component(unsigned int n, IceModelVec2S &source);
  inline double& operator() (int i, int j, int k) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, k);
#endif
    return static_cast<double***>(array)[i][j][k];
  }
  virtual PetscErrorCode create(IceGrid &my_grid, std::string my_short_name,
                                IceModelVecKind ghostedp, unsigned int stencil_width, int dof);
protected:
  PetscErrorCode get_component(unsigned int n, Vec result);
  PetscErrorCode set_component(unsigned int n, Vec source);
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
  virtual PetscErrorCode  create(IceGrid &my_grid, std::string my_name,
                                 IceModelVecKind ghostedp, int width = 1);
  virtual PetscErrorCode  put_on_proc0(Vec onp0, VecScatter ctx, Vec g2, Vec g2natural);
  virtual PetscErrorCode  get_from_proc0(Vec onp0, VecScatter ctx, Vec g2, Vec g2natural);
  using IceModelVec::copy_to;
  using IceModelVec::copy_from;
  virtual PetscErrorCode  copy_to(IceModelVec &destination);
  virtual PetscErrorCode  copy_from(IceModelVec &source);
  PetscErrorCode  get_array(double** &a);
  virtual PetscErrorCode set_to_magnitude(IceModelVec2S &v_x, IceModelVec2S &v_y);
  virtual PetscErrorCode mask_by(IceModelVec2S &M, double fill = 0.0);
  virtual PetscErrorCode add(double alpha, IceModelVec &x);
  virtual PetscErrorCode add(double alpha, IceModelVec &x, IceModelVec &result);
  virtual PetscErrorCode sum(double &result);
  virtual PetscErrorCode min(double &result);
  virtual PetscErrorCode max(double &result);
  virtual PetscErrorCode absmax(double &result);
  virtual double diff_x(int i, int j);
  virtual double diff_y(int i, int j);
  virtual double diff_x_stagE(int i, int j);
  virtual double diff_y_stagE(int i, int j);
  virtual double diff_x_stagN(int i, int j);
  virtual double diff_y_stagN(int i, int j);
  virtual double diff_x_p(int i, int j);
  virtual double diff_y_p(int i, int j);
  virtual PetscErrorCode view_matlab(PetscViewer my_viewer);
  virtual PetscErrorCode has_nan();

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
  inline int as_int(int i, int j) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, 0);
#endif
    const double **a = (const double**) array;
    return static_cast<int>(floor(a[i][j] + 0.5));
  }

  inline planeStar<int> int_star(int i, int j) {
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
class PISMVector2 {
public:
  PISMVector2() : u(0), v(0) {}
  PISMVector2(double a, double b) : u(a), v(b) {}

  //! Magnitude squared.
  inline double magnitude_squared() const {
    return u*u + v*v;
  }
  //! Magnitude.
  inline double magnitude() const {
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

  inline PISMVector2& operator*=(const double &a) {
    u *= a;
    v *= a;
    return *this;
  }

  inline PISMVector2& operator/=(const double &a) {
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
  inline PISMVector2 operator*(const double &a) const {
    return PISMVector2(u * a, v * a);
  }

  //! \brief Scales a vector.
  inline PISMVector2 operator/(const double &a) const {
    return PISMVector2(u / a, v / a);
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
  virtual PetscErrorCode create(IceGrid &my_grid, std::string my_short_name,
				IceModelVecKind ghostedp, unsigned int stencil_width = 1);
  using IceModelVec::copy_to;
  using IceModelVec::copy_from;
  virtual PetscErrorCode copy_to(IceModelVec &destination);
  virtual PetscErrorCode copy_from(IceModelVec &source);
  virtual PetscErrorCode add(double alpha, IceModelVec &x);
  virtual PetscErrorCode add(double alpha, IceModelVec &x, IceModelVec &result);

  // I/O:
  using IceModelVec2::write;
  virtual PetscErrorCode get_array(PISMVector2 ** &a);
  virtual PetscErrorCode magnitude(IceModelVec2S &result);
  inline PISMVector2& operator()(int i, int j) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, 0);
#endif
    return static_cast<PISMVector2**>(array)[i][j];
  }

  inline planeStar<PISMVector2> star(int i, int j) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, 0);
    check_array_indices(i+1, j, 0);
    check_array_indices(i-1, j, 0);
    check_array_indices(i, j+1, 0);
    check_array_indices(i, j-1, 0);
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
  virtual PetscErrorCode set_name(std::string name, int component = 0);
  virtual PetscErrorCode rename(std::string short_name, std::string long_name,
                                std::string standard_name, int component = 0);
  virtual PetscErrorCode rename(std::string short_name, std::vector<std::string> long_names,
                                std::string standard_name);
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
  virtual PetscErrorCode create(IceGrid &my_grid, std::string my_short_name, IceModelVecKind ghostedp,
                                unsigned int stencil_width = 1);
  virtual PetscErrorCode staggered_to_regular(IceModelVec2S &result);
  virtual PetscErrorCode staggered_to_regular(IceModelVec2V &result);
  virtual PetscErrorCode absmaxcomponents(double* z);

  //! Returns the values at interfaces of the cell i,j using the staggered grid.
  /*! The ij member of the return value is set to 0, since it has no meaning in
    this context.
   */
  inline planeStar<double> star(int i, int j) {
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

  PetscErrorCode  setColumn(int i, int j, double c);
  PetscErrorCode  setInternalColumn(int i, int j, double *valsIN);
  PetscErrorCode  getInternalColumn(int i, int j, double **valsOUT);

  virtual double    getValZ(int i, int j, double z);
  virtual PetscErrorCode isLegalLevel(double z);

  inline double& operator() (int i, int j, int k) {
#if (PISM_DEBUG==1)
    check_array_indices(i, j, k);
#endif
    return static_cast<double***>(array)[i][j][k];
  }

protected:
  virtual PetscErrorCode allocate(IceGrid &mygrid, std::string my_short_name,
                                  IceModelVecKind ghostedp, std::vector<double> levels,
                                  unsigned int stencil_width = 1);
  virtual PetscErrorCode has_nan();
};


//! Class for a 3d DA-based Vec for ice scalar quantities.
class IceModelVec3 : public IceModelVec3D {
public:
  IceModelVec3() {}
  virtual ~IceModelVec3() {}

  virtual PetscErrorCode create(IceGrid &mygrid, std::string my_short_name,
                                IceModelVecKind ghostedp,
                                unsigned int stencil_width = 1);

  // need to call begin_access() before set...(i,j,...) or get...(i,j,...) *and* need call
  // end_access() afterward
  PetscErrorCode  getValColumn(int i, int j, unsigned int ks, double *valsOUT);
  PetscErrorCode  getValColumnQUAD(int i, int j, unsigned int ks, double *valsOUT);
  PetscErrorCode  getValColumnPL(int i, int j, unsigned int ks, double *valsOUT);

  PetscErrorCode  setValColumnPL(int i, int j, double *valsIN);

  PetscErrorCode  getPlaneStarZ(int i, int j, double z,
                                planeStar<double> *star);
  PetscErrorCode  getPlaneStar_fine(int i, int j, unsigned int k,
				    planeStar<double> *star);
  PetscErrorCode  getPlaneStar(int i, int j, unsigned int k,
			       planeStar<double> *star);

  PetscErrorCode  getHorSlice(Vec &gslice, double z); // used in iMmatlab.cc
  PetscErrorCode  getHorSlice(IceModelVec2S &gslice, double z);
  PetscErrorCode  getSurfaceValues(Vec &gsurf, IceModelVec2S &myH); // used in iMviewers.cc
  PetscErrorCode  getSurfaceValues(IceModelVec2S &gsurf, IceModelVec2S &myH);
  PetscErrorCode  getSurfaceValues(IceModelVec2S &gsurf, double **H);
  PetscErrorCode  extend_vertically(int old_Mz, double fill_value);
  PetscErrorCode  extend_vertically(int old_Mz, IceModelVec2S &fill_values);
protected:
  virtual PetscErrorCode  extend_vertically_private(int old_Mz);
};

#endif /* __IceModelVec_hh */

