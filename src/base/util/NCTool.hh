// Copyright (C) 2007--2011 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of Pism.
//
// Pism is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// Pism is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with Pism; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef __NCTool_hh
#define __NCTool_hh

#include <petsc.h>

// The following is a stupid kludge necessary to make NetCDF 4.x work in
// serial mode in an MPI program:
#ifndef MPI_INCLUDED
#define MPI_INCLUDED 1
#endif
#include <netcdf.h>		// nc_type

// Note: as far as I (CK) can tell, MPI_INCLUDED is a MPICH invention.

#include "udunits.h"	// utUnit
#include <string>
#include <vector>
#include <map>

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

enum AxisType {X_AXIS, Y_AXIS, Z_AXIS, T_AXIS, UNKNOWN_AXIS};

int nc_check(int stat);
int check_err(const int stat, const int line, const char *file);

//! A fairly low-level NetCDF wrapper used in various PISM classes.
/*!
  \section nctool_overview PISM's NetCDF wrapper

  This class does not contain any grid-specific code to make it possible to
  use it to access NetCDF files that do not contain gridded data.

  All the work is done on processor 0. Other processors get results matching
  ones on processor zero, unless stated otherwise.

  Most direct uses of this class are related to determining grid parameters if
  an input file. In the majority of other cases one can use higher-level
  classes, such as IceModelVec.

  Note: "non-const" methods of this this class are the only ones that can make
  a particular instance NCTool refer to a different file. In other words,
  passing "const NCTool &" as an argument ensures that after a call the file is
  still open.
 */
class NCTool {
public:
  NCTool(MPI_Comm c, PetscMPIInt r);
  virtual ~NCTool();
  virtual PetscErrorCode open_for_reading(string filename);
  virtual PetscErrorCode open_for_writing(string filename);
  virtual PetscErrorCode move_if_exists(string filename);
  virtual PetscErrorCode close();

  virtual PetscErrorCode find_variable(string short_name, string standard_name,
			       int *varid, bool &exists, bool &found_by_standard_name) const;
  virtual PetscErrorCode find_variable(string short_name, string standard_name,
			       int *varid, bool &exists) const;
  virtual PetscErrorCode find_variable(string short_name, int *varid, bool &exists) const;
  virtual PetscErrorCode find_dimension(string short_name, int *dimid, bool &exists) const;
  virtual PetscErrorCode append_time(string name, PetscReal time) const;
  virtual PetscErrorCode append_time_bounds(string name, PetscReal t0, PetscReal t1) const;
  virtual PetscErrorCode write_history(string history, bool overwrite = false) const;
  virtual bool check_dimension(string dim, int len) const;

  virtual PetscErrorCode get_dim_length(string name, unsigned int *len) const;
  virtual PetscErrorCode get_dim_limits(string name, double *min, double *max) const;

  virtual PetscErrorCode get_dimension(string name, vector<double> &result) const;
  virtual PetscErrorCode put_dimension(int varid, const vector<double> &vals) const;

  virtual PetscErrorCode inq_unlimdim(int &unlimdimid) const;
  virtual PetscErrorCode inq_dimname(int dimid, string &name) const;
  virtual PetscErrorCode inq_dimtype(string name, AxisType &result) const;
  virtual PetscErrorCode inq_dimids(int varid, vector<int> &dimids) const;
  virtual PetscErrorCode inq_nattrs(int varid, int &N) const;
  virtual PetscErrorCode inq_att_name(int varid, int n, string &name) const;
  virtual PetscErrorCode inq_att_type(int varid, string name, nc_type &typep) const;
  virtual PetscErrorCode get_att_text(int varid, string name, string &result) const;
  virtual PetscErrorCode get_att_double(int varid, string name, vector<double> &result) const;
  virtual PetscErrorCode get_units(int varid, bool &has_units, utUnit &units) const;
  virtual PetscErrorCode get_nrecords(int &nrecords) const;
  virtual PetscErrorCode get_nrecords(string short_name, string std_name,
                                      unsigned int &nrecords) const;
  virtual PetscErrorCode set_attrs(int varid, map<string,string> attrs) const;
  virtual PetscErrorCode create_dimension(string name, int length, map<string,string> attrs,
                                          int &dimid, int &varid) const;

  virtual int get_ncid() const;
  virtual PetscErrorCode define_mode() const;
  virtual PetscErrorCode data_mode() const;
protected:
  int ncid;
  MPI_Comm com;
  PetscMPIInt rank;
  mutable bool def_mode;             // note: only processor 0 should use this
};

#endif // __NCTool_hh

