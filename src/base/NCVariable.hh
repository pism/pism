// Copyright (C) 2009 Constantine Khroulev
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

#ifndef __NCVariable_hh
#define __NCVariable_hh

#include <map>
#include <vector>
#include <string>
#include <petscda.h>
#include <netcdf.h>
#include "../udunits/udunits.h"
#include "nc_util.hh"

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

struct NCVariable {
public:
  NCVariable();
  void init(const char name[], IceGrid &g, GridType d);
  PetscErrorCode set_units(const char*);
  PetscErrorCode set_glaciological_units(const char*);
  PetscErrorCode write(const char filename[], nc_type nctype,
		       bool write_in_glaciological_units, Vec v);
  PetscErrorCode read(const char filename[], unsigned int time, Vec v);
  PetscErrorCode regrid(const char filename[], LocalInterpCtx &lic,
			bool critical, bool set_default_value,
			PetscScalar default_value, MaskInterp *, Vec v);
  PetscErrorCode change_units(Vec v, utUnit *from, utUnit *to);
  PetscErrorCode reset();
  void set(const char*, double);
  double get(const char*);
  bool has(const char name[]);
			
  string short_name;

  // Attributes:
  /*!
    \li long_name
    \li standard_name
    \li pism_intent
    \li units
    \li glaciological_units
   */
  map<string, string> strings;
  map<string, vector<double> > doubles;

protected:
  PetscErrorCode write_attributes(int ncid, int varid, nc_type nctype,
				  bool write_in_glaciological_units);
  PetscErrorCode report_range(Vec v);
  PetscErrorCode check_range(Vec v);
  PetscErrorCode read_valid_range(int ncid, int varid);
  PetscErrorCode define(int ncid, nc_type nctype, int &varid);
  IceGrid *grid;
  GridType dims;
  utUnit units,		      //!< internal (PISM) units
         glaciological_units; //!< for diagnostic variables: units to use when writing
			      //!< to a NetCDF file and for standard out reports
};

#endif
