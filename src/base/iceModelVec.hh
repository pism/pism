// Copyright (C) 2008, 2009 Ed Bueler and Constantine Khroulev
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

#include <cstring>
#include <cstdlib>
#include <petscda.h>
#include <netcdf.h>
#include "../udunits/udunits.h"
#include "nc_util.hh"
#include "NCVariable.hh"
#include "pism_const.hh"

#ifndef __IceModelVec_hh
#define __IceModelVec_hh

//! Abstract class for reading, writing, allocating, and accessing a DA-based PETSc Vec from within IceModel.
class IceModelVec {
public:
  IceModelVec();
  IceModelVec(const IceModelVec &other);
  virtual ~IceModelVec();

  virtual PetscErrorCode  create(IceGrid &mygrid, const char my_short_name[], bool local);
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
  virtual PetscErrorCode  copy_to_global(Vec destination);
  virtual PetscErrorCode  copy_to(IceModelVec &destination);
  virtual PetscErrorCode  copy_from(IceModelVec &source);
  virtual PetscErrorCode  set_name(const char name[]);
  virtual PetscErrorCode  set_glaciological_units(string units);
  virtual PetscErrorCode  set_attr(string name, string value);
  virtual PetscErrorCode  set_attr(string name, double value);
  virtual PetscErrorCode  set_attr(string name, vector<double> values);
  virtual bool            has_attr(string name);
  virtual string          string_attr(string name);
  virtual double          double_attr(string name);
  virtual vector<double>  array_attr(string name);
  virtual PetscErrorCode  set_attrs(string my_pism_intent, string my_long_name,
				    string my_units, string my_standard_name);
  virtual bool            is_valid(PetscScalar a);
  virtual PetscErrorCode  write(const char filename[]);
  virtual PetscErrorCode  write(const char filename[], nc_type nctype);
  virtual PetscErrorCode  read(const char filename[], unsigned int time);
  virtual PetscErrorCode  regrid(const char filename[], LocalInterpCtx &lic, bool critical);
  virtual PetscErrorCode  regrid(const char filename[], LocalInterpCtx &lic, PetscScalar default_value);

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

  NCSpatialVariable var1;	//!< a NetCDF variable corresponding to this
				//!IceModelVec; called var1 because some
				//!IceModelVecs will have more: var2, etc.

  IceGrid      *grid;
  GridType     dims;
  int dof;
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


//! Class for a 2d DA-based Vec for ice scalar quantities in IceModel.
class IceModelVec2 : public IceModelVec {
public:
  IceModelVec2();
  IceModelVec2(const IceModelVec2 &other);
  // does not need a copy constructor, because it does not add any new data members
  virtual PetscErrorCode  create(IceGrid &my_grid, const char my_short_name[], bool local);
  PetscErrorCode  create(IceGrid &my_grid, const char my_short_name[], bool local,
			 DAStencilType my_sten, int stencil_width);
  virtual PetscErrorCode  put_on_proc0(Vec onp0, VecScatter ctx, Vec g2, Vec g2natural);
  virtual PetscErrorCode  get_from_proc0(Vec onp0, VecScatter ctx, Vec g2, Vec g2natural);
  PetscErrorCode  get_array(PetscScalar** &a);
  virtual PetscErrorCode set_to_magnitude(IceModelVec2 &v_x, IceModelVec2 &v_y);
  virtual PetscErrorCode mask_by(IceModelVec2 &M, PetscScalar fill = 0.0);
  virtual PetscErrorCode view(Vec g2, PetscInt viewer_size);
  virtual PetscScalar& operator() (int i, int j);
};

//! \brief A simple class "hiding" the fact that the mask is stored as
//! floating-point scalars (instead of integers).
class IceModelVec2Mask : public IceModelVec2 {
public:
  virtual PismMask value(int i, int j);	  // returns the mask value
  virtual bool is_grounded(int i, int j); // checks for MASK_SHEET || MASK_DRAGGING
  virtual bool is_floating(int i, int j); // checks for MASK_FLOATING || MASK_FLOATING_OCEAN0
  virtual PetscErrorCode fill_where_grounded(IceModelVec2 &v, const PetscScalar fillval);
  virtual PetscErrorCode fill_where_floating(IceModelVec2 &v, const PetscScalar fillval);
};

//! Class for a 3d DA-based Vec for bedrock (lithosphere) scalar quantities in IceModel.
class IceModelVec3Bedrock : public IceModelVec {
public:
  IceModelVec3Bedrock();
  IceModelVec3Bedrock(const IceModelVec3Bedrock &other);
  virtual PetscErrorCode create(IceGrid &mygrid, const char my_short_name[], bool local);

  PetscErrorCode  setInternalColumn(PetscInt i, PetscInt j, PetscScalar *valsIN);
  PetscErrorCode  setColumn(PetscInt i, PetscInt j, PetscScalar c);
  PetscErrorCode  getInternalColumn(PetscInt i, PetscInt j, PetscScalar **valsOUT);

  PetscErrorCode  setValColumnPL(PetscInt i, PetscInt j, PetscInt nlevels, 
				 PetscScalar *levelsIN, PetscScalar *valsIN);
  PetscErrorCode  getValColumnPL(PetscInt i, PetscInt j, PetscInt nlevels, 
				 PetscScalar *levelsIN, PetscScalar *valsOUT);
  PetscErrorCode  getValColumnQUAD(PetscInt i, PetscInt j, PetscInt nlevels, 
				   PetscScalar *levelsIN, PetscScalar *valsOUT);
  PetscErrorCode view_sounding(int i, int j, PetscInt viewer_size);

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
  virtual PetscErrorCode  create(IceGrid &mygrid, const char my_short_name[], bool local);

  // note the IceModelVec3 with this method must be *local* while imv3_source must be *global*
  virtual PetscErrorCode  beginGhostCommTransfer(IceModelVec3 &imv3_source);
  virtual PetscErrorCode  endGhostCommTransfer(IceModelVec3 &imv3_source);

  // need call begin_access() before set...() or get...() *and* need call end_access() afterward
  PetscErrorCode  setValColumnPL(PetscInt i, PetscInt j, PetscInt nlevels, 
                               PetscScalar *levelsIN, PetscScalar *valsIN);
  PetscErrorCode  setColumn(PetscInt i, PetscInt j, PetscScalar c);
  PetscErrorCode  setInternalColumn(PetscInt i, PetscInt j, PetscScalar *valsIN);

  PetscScalar     getValZ(PetscInt i, PetscInt j, PetscScalar z);

  PetscErrorCode  getPlaneStarZ(PetscInt i, PetscInt j, PetscScalar z,
                                planeStar *star);

//Idea; see IceModel::horizontalVelocitySIARegular(), for example; would only work for 4 neighbors
//   if used with two IceModelVec3:
//  PetscErrorCode  averageStagNeighCol(PetscInt i, PetscInt j, 
//                                      PetscInt nlevels, PetscScalar *levelsIN, 
//                                      PetscScalar *avOUT);

  PetscErrorCode  getValColumnPL(PetscInt i, PetscInt j, PetscInt nlevels, 
                                 PetscScalar *levelsIN, PetscScalar *valsOUT);
  PetscErrorCode  getValColumnQUAD(PetscInt i, PetscInt j, PetscInt nlevels, 
                                   PetscScalar *levelsIN, PetscScalar *valsOUT);
  PetscErrorCode  getInternalColumn(PetscInt i, PetscInt j, PetscScalar **valsPTR);

  PetscErrorCode  getHorSlice(Vec &gslice, PetscScalar z); // used in iMmatlab.cc
  PetscErrorCode  getHorSlice(IceModelVec2 &gslice, PetscScalar z);
  PetscErrorCode  getSurfaceValues(Vec &gsurf, IceModelVec2 &myH); // used in iMviewers.cc
  PetscErrorCode  getSurfaceValues(IceModelVec2 &gsurf, IceModelVec2 &myH);
  PetscErrorCode  getSurfaceValues(IceModelVec2 &gsurf, PetscScalar **H);
  PetscErrorCode  extend_vertically(int old_Mz, PetscScalar fill_value);
  PetscErrorCode  extend_vertically(int old_Mz, IceModelVec2 &fill_values);

  PetscErrorCode view_surface(IceModelVec2 &thickness, Vec g2, PetscInt viewer_size);
  PetscErrorCode view_horizontal_slice(PetscScalar level, Vec g2, PetscInt viewer_size);
  PetscErrorCode view_sounding(int i, int j, PetscInt viewer_size);

protected:  
  virtual PetscErrorCode  destroy();
  PetscErrorCode  isLegalLevel(PetscScalar z);
  virtual PetscErrorCode  extend_vertically_private(int old_Mz);
  map<string,PetscViewer> *slice_viewers, *sounding_viewers;
  Vec sounding_buffer;
};

#endif /* __IceModelVec_hh */

