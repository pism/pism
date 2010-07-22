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

#ifndef __IceModelVec_hh
#define __IceModelVec_hh

#include <cstring>
#include <cstdlib>
#include <petscda.h>
#include "../udunits/udunits.h"
#include "nc_util.hh"
#include "NCSpatialVariable.hh"
#include "pism_const.hh"
#include "LocalInterpCtx.hh"

//! Abstract class for reading, writing, allocating, and accessing a DA-based PETSc Vec from within IceModel.
class IceModelVec {
public:
  IceModelVec();
  IceModelVec(const IceModelVec &other);
  virtual ~IceModelVec();

  virtual bool            was_created();
  virtual GridType        grid_type();

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
  virtual PetscErrorCode  set_name(const char name[], int component = 0);
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
  virtual bool            is_valid(PetscScalar a, int component = 0);
  virtual PetscErrorCode  write(const char filename[]);
  virtual PetscErrorCode  write(const char filename[], nc_type nctype);
  virtual PetscErrorCode  read(const char filename[], unsigned int time);
  virtual PetscErrorCode  regrid(const char filename[], LocalInterpCtx &lic,
                                 bool critical);
  virtual PetscErrorCode  regrid(const char filename[], LocalInterpCtx &lic, 
                                 PetscScalar default_value);

  virtual PetscErrorCode  begin_access();
  virtual PetscErrorCode  end_access();
  virtual PetscErrorCode  beginGhostComm();
  virtual PetscErrorCode  endGhostComm();
  virtual PetscErrorCode  beginGhostComm(IceModelVec &destination);
  virtual PetscErrorCode  endGhostComm(IceModelVec &destination);

  virtual PetscErrorCode  set(PetscScalar c);
 
  MaskInterp interpolation_mask;
  bool   use_interpolation_mask, write_in_glaciological_units, time_independent;
  nc_type output_data_type;
protected:

  bool shallow_copy;
  Vec  v;
  string name;

  //! NetCDF variable(s) corresponding to this IceModelVec; dof == 1 vectors
  //! only have vars[0].
  vector<NCSpatialVariable> vars;

  IceGrid      *grid;
  GridType     dims;
  int          dof, s_width;
  DA           da;
  bool         localp;
  //! it is a map, because a temporary IceModelVec can be used to view
  //! different quantities, and a pointer because "shallow copies" should have
  //! the acces to the original map
  map<string,PetscViewer> *map_viewers;

  void         *array;  // will be PetscScalar** or PetscScalar*** in derived classes

  int access_counter;		// used in begin_access() and end_access()

  virtual PetscErrorCode destroy();
  virtual PetscErrorCode checkAllocated();
  virtual PetscErrorCode checkHaveArray();
  virtual PetscErrorCode checkCompatibility(const char*, IceModelVec &other);
  virtual PetscErrorCode reset_attrs();
  virtual PetscErrorCode create_viewer(PetscInt viewer_size, string title, PetscViewer &viewer);
  virtual PetscErrorCode compute_viewer_size(int target, int &x, int &y);
  // FIXME: consider adding 
  //   virtual PetscErrorCode  checkSelfOwnsIt(const PetscInt i, const PetscInt j);
  //   virtual PetscErrorCode  checkSelfOwnsItGhosted(const PetscInt i, const PetscInt j);
};


struct planeStar {
  PetscScalar ij, ip1, im1, jp1, jm1;
};

struct planeBox {
  PetscScalar ij, ip1, im1, jp1, jm1, ip1jp1, im1jp1, ip1jm1, im1jm1;
};

//! Class for a 2d DA-based Vec for scalar quantities in IceModel.
class IceModelVec2 : public IceModelVec {
public:
  IceModelVec2() : IceModelVec() {}
  IceModelVec2(const IceModelVec2 &other) : IceModelVec(other) {};
  virtual PetscErrorCode create(IceGrid &my_grid, const char my_short_name[], bool local,
			 DAStencilType my_sten, int stencil_width, int dof);
  virtual PetscErrorCode view(PetscInt viewer_size);
protected:
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
  virtual PetscErrorCode  create(IceGrid &my_grid, const char my_name[], bool local, int width = 1);
  virtual PetscErrorCode  put_on_proc0(Vec onp0, VecScatter ctx, Vec g2, Vec g2natural);
  virtual PetscErrorCode  get_from_proc0(Vec onp0, VecScatter ctx, Vec g2, Vec g2natural);
  PetscErrorCode  get_array(PetscScalar** &a);
  virtual PetscErrorCode set_to_magnitude(IceModelVec2S &v_x, IceModelVec2S &v_y);
  virtual PetscErrorCode mask_by(IceModelVec2S &M, PetscScalar fill = 0.0);
  virtual PetscErrorCode sum(PetscScalar &result);
  virtual PetscScalar diff_x(int i, int j);
  virtual PetscScalar diff_y(int i, int j);
  virtual PetscScalar diff_x_stagE(int i, int j);
  virtual PetscScalar diff_y_stagE(int i, int j);
  virtual PetscScalar diff_x_stagN(int i, int j);
  virtual PetscScalar diff_y_stagN(int i, int j);
  virtual PetscScalar diff_x_p(int i, int j);
  virtual PetscScalar diff_y_p(int i, int j);
  virtual PetscErrorCode view_matlab(PetscViewer my_viewer);
  virtual PetscScalar& operator() (int i, int j);
  virtual PetscErrorCode has_nan();
};

//! \brief A class for storing and accessing internal staggered-grid 2D fields.
//! Uses dof=2 storage. Does \b not support input and output. This class is
//! identical to IceModelVec2V, except that components are not called \c u and
//! \c v (to avoid confusion).
class IceModelVec2Stag : public IceModelVec2 {
public:
  IceModelVec2Stag() { dof = 2; vars.resize(dof); }
  IceModelVec2Stag(const IceModelVec2S &other) : IceModelVec2(other) {}
  using IceModelVec2::create;
  virtual PetscErrorCode create(IceGrid &my_grid, const char my_name[], bool local, int width = 1);
  virtual PetscErrorCode get_array(PetscScalar*** &a);
  virtual PetscErrorCode begin_access();
  virtual PetscErrorCode end_access();
  virtual PetscScalar& operator() (int i, int j, int k);
  // component-wise access:
  virtual PetscErrorCode get_component(int n, IceModelVec2S &result);
  virtual PetscErrorCode set_component(int n, IceModelVec2S &source);
};

//! \brief A simple class "hiding" the fact that the mask is stored as
//! floating-point scalars (instead of integers).
class IceModelVec2Mask : public IceModelVec2 {
public:
  using IceModelVec2::create;
  virtual PetscErrorCode create(IceGrid &my_grid, const char my_name[], bool local, int width = 1);
  PetscErrorCode  get_array(PetscScalar** &a); // provides access to the storage (PetscScalar) array
  virtual PetscScalar& operator() (int i, int j);
  virtual PismMask value(int i, int j);	  // returns the mask value
  virtual bool is_grounded(int i, int j); // checks for MASK_SHEET || MASK_DRAGGING
  virtual bool is_floating(int i, int j); // checks for MASK_FLOATING || MASK_FLOATING_OCEAN0
  virtual PetscErrorCode fill_where_grounded(IceModelVec2S &v, const PetscScalar fillval);
  virtual PetscErrorCode fill_where_floating(IceModelVec2S &v, const PetscScalar fillval);
};

/// IceModeVec2V
struct PISMVector2 {
  PetscScalar u, v;
};

//! Class for storing and accessing 2D vector fields used in IceModel.
class IceModelVec2V : public IceModelVec2 {
public:
  IceModelVec2V();
  IceModelVec2V(const IceModelVec2V &original);
  ~IceModelVec2V();
  using IceModelVec2::create;
  virtual PetscErrorCode create(IceGrid &my_grid, const char my_short_name[],
				bool local, int stencil_width = 1); 

  // I/O:
  using IceModelVec2::write;
  virtual PetscErrorCode write(const char filename[], nc_type nctype); 
  virtual PetscErrorCode read(const char filename[], const unsigned int time); 
  virtual PetscErrorCode regrid(const char filename[], LocalInterpCtx &lic, bool critical); 
  virtual PetscErrorCode regrid(const char filename[], LocalInterpCtx &lic, PetscScalar default_value); 
  virtual PetscErrorCode get_array(PISMVector2 ** &a); 
  virtual PetscErrorCode magnitude(IceModelVec2S &result); 
  virtual PISMVector2&   operator()(int i, int j);
  virtual PetscErrorCode view(PetscInt viewer_size);
  // component-wise access:
  virtual PetscErrorCode get_component(int n, IceModelVec2S &result);
  virtual PetscErrorCode set_component(int n, IceModelVec2S &source);
  // Metadata, etc:
  using IceModelVec2::is_valid;
  virtual bool           is_valid(PetscScalar u, PetscScalar v); 
protected:
  DA component_da;		//!< same as \c da, but for one component only
  PetscErrorCode destroy();	
};

//! Class for a 3d DA-based Vec for bedrock (lithosphere) scalar quantities in IceModel.
class IceModelVec3Bedrock : public IceModelVec {
public:
  IceModelVec3Bedrock();
  IceModelVec3Bedrock(const IceModelVec3Bedrock &other);
  virtual ~IceModelVec3Bedrock();
  virtual PetscErrorCode create(IceGrid &mygrid, const char my_short_name[], bool local);

  PetscErrorCode  setInternalColumn(PetscInt i, PetscInt j, PetscScalar *valsIN);
  PetscErrorCode  setColumn(PetscInt i, PetscInt j, PetscScalar c);
  PetscErrorCode  getInternalColumn(PetscInt i, PetscInt j, PetscScalar **valsOUT);
  PetscErrorCode  setValColumnPL(PetscInt i, PetscInt j, PetscScalar *valsIN);
  PetscErrorCode  getValColumnPL(PetscInt i, PetscInt j, PetscScalar *valsOUT);
  PetscErrorCode  view_sounding(int i, int j, PetscInt viewer_size);

protected:  
  map<string,PetscViewer> *sounding_viewers;
  Vec sounding_buffer;
  PetscErrorCode  isLegalLevel(PetscScalar z);
  virtual PetscErrorCode destroy();
};


//! Class for a 3d DA-based Vec for ice scalar quantities in IceModel.
class IceModelVec3 : public IceModelVec {
public:
  IceModelVec3();
  IceModelVec3(const IceModelVec3 &other);
  virtual ~IceModelVec3();
  virtual PetscErrorCode  create(IceGrid &mygrid, const char my_short_name[],
				 bool local, int stencil_width = 1);

  // note the IceModelVec3 with this method must be *local* while imv3_source must be *global*
  virtual PetscErrorCode beginGhostCommTransfer(IceModelVec3 &imv3_source);
  virtual PetscErrorCode endGhostCommTransfer(IceModelVec3 &imv3_source);

  // need call begin_access() before set...() or get...() *and* need call end_access() afterward
  PetscErrorCode  setValColumnPL(PetscInt i, PetscInt j, PetscScalar *valsIN);
  PetscErrorCode  setColumn(PetscInt i, PetscInt j, PetscScalar c);
  PetscErrorCode  setInternalColumn(PetscInt i, PetscInt j, PetscScalar *valsIN);

  PetscScalar     getValZ(PetscInt i, PetscInt j, PetscScalar z);

  PetscErrorCode  getPlaneStarZ(PetscInt i, PetscInt j, PetscScalar z,
                                planeStar *star);
  PetscErrorCode  getPlaneStar_fine(PetscInt i, PetscInt j, PetscInt k,
				    planeStar *star);
  PetscErrorCode  getPlaneStar(PetscInt i, PetscInt j, PetscInt k,
			       planeStar *star);

//Idea; see IceModel::horizontalVelocitySIARegular(), for example; would only work for 4 neighbors
//   if used with two IceModelVec3:
//  PetscErrorCode  averageStagNeighCol(PetscInt i, PetscInt j, 
//                                      PetscInt nlevels, PetscScalar *levelsIN, 
//                                      PetscScalar *avOUT);

  PetscErrorCode  getValColumnPL(PetscInt i, PetscInt j, PetscInt ks, PetscScalar *valsOUT);
  PetscErrorCode  getValColumnQUAD(PetscInt i, PetscInt j, PetscInt ks, PetscScalar *valsOUT);
  PetscErrorCode  getValColumn(PetscInt i, PetscInt j, PetscInt ks, PetscScalar *valsOUT);
  PetscErrorCode  getInternalColumn(PetscInt i, PetscInt j, PetscScalar **valsPTR);

  PetscErrorCode  getHorSlice(Vec &gslice, PetscScalar z); // used in iMmatlab.cc
  PetscErrorCode  getHorSlice(IceModelVec2S &gslice, PetscScalar z);
  PetscErrorCode  getSurfaceValues(Vec &gsurf, IceModelVec2S &myH); // used in iMviewers.cc
  PetscErrorCode  getSurfaceValues(IceModelVec2S &gsurf, IceModelVec2S &myH);
  PetscErrorCode  getSurfaceValues(IceModelVec2S &gsurf, PetscScalar **H);
  PetscErrorCode  extend_vertically(int old_Mz, PetscScalar fill_value);
  PetscErrorCode  extend_vertically(int old_Mz, IceModelVec2S &fill_values);

  PetscErrorCode view_surface(IceModelVec2S &thickness, PetscInt viewer_size);
  PetscErrorCode view_horizontal_slice(PetscScalar level, PetscInt viewer_size);
  PetscErrorCode view_sounding(int i, int j, PetscInt viewer_size);
  virtual PetscErrorCode  has_nan();

protected:  
  virtual PetscErrorCode  destroy();
  PetscErrorCode  isLegalLevel(PetscScalar z);
  virtual PetscErrorCode  extend_vertically_private(int old_Mz);
  virtual PetscErrorCode  create_da(DA &result, PetscInt Mz);
  map<string,PetscViewer> *slice_viewers, *sounding_viewers;
  Vec sounding_buffer;
};

#endif /* __IceModelVec_hh */

