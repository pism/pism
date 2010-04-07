// Copyright (C) 2009, 2010 Constantine Khroulev
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


//! \brief A class for handling variable metadata, reading, writing and converting
//! from input units and to output units.
/*! A NetCDF variable can have any number of attributes, but some of them are
  treated differently:

  \li units: specifies internal units. When read, a variable is converted to
  these units. When written, it is converted from these to glaciological_units
  if write_in_glaciological_units is true.

  \li glaciological_units: is never written to a file; replaces 'units' in the
  output if write_in_glaciological_units is true.

  \li valid_min, valid_max: specify the valid range of a variable. Are read
  from an input file \b only if not specified previously. If both are set, then
  valid_range is used in the output instead.

  Also:

  \li empty string attributes are ignored (they are not written to the output;
  file and has("foo") returns false if "foo" is absent or equal to an empty
  string).
 */
class NCVariable {
public:
  NCVariable();
  virtual ~NCVariable() {};
  void init(string name, MPI_Comm c, PetscMPIInt r);
  virtual PetscErrorCode set_units(string);
  virtual PetscErrorCode set_glaciological_units(string);
  virtual PetscErrorCode reset();
  virtual void set(string, double);
  virtual double get(string) const;
  virtual void set_string(string name, string value);
  virtual string get_string(string) const;
  virtual bool has(string) const;
  virtual bool is_valid(PetscScalar a) const;
			
  string short_name;

  // Attributes:
  /*! Typical attributes stored here:

    \li long_name
    \li standard_name
    \li pism_intent
    \li units
    \li glaciological_units
   */
  map<string, vector<double> > doubles;

protected:
  virtual PetscErrorCode write_attributes(const NCTool &nc, int varid, nc_type nctype,
					  bool write_in_glaciological_units) const;
  virtual PetscErrorCode read_valid_range(const NCTool &nc, int varid);
  MPI_Comm com;
  PetscMPIInt rank;
  map<string, string> strings;
  utUnit units,		      //!< internal (PISM) units
         glaciological_units; //!< \brief for diagnostic variables: units to
			      //!use when writing to a NetCDF file and for
			      //!standard out reports
};

//! A class for reading, writing and accessing PISM configuration flags and parameters.
class NCConfigVariable : public NCVariable {
public:
  virtual PetscErrorCode print() const;
  virtual PetscErrorCode read(const char filename[]);
  virtual PetscErrorCode write(const char filename[]) const;
  virtual double get(string) const;
  virtual bool   get_flag(string) const;
  virtual string get_string(string name) const;
  virtual void   set_flag(string, bool);
  virtual PetscErrorCode flag_from_option(string, string);
  virtual PetscErrorCode scalar_from_option(string, string);
  virtual void import_from(const NCConfigVariable &other);
  virtual void update_from(const NCConfigVariable &other);
protected:
  string config_filename;
  virtual PetscErrorCode write_attributes(const NCTool &nc, int varid, nc_type nctype,
					  bool write_in_glaciological_units) const;
  virtual PetscErrorCode define(int ncid, int &varid) const;
};

//! A class for reading and writing NetCDF global attributes.
/*! This is not a variable, because it has no value, but it is similar to
  NCConfigVariable, because uses of these attributes are similar.
*/
class NCGlobalAttributes : public NCConfigVariable {
public:
  virtual PetscErrorCode read(const char filename[]);
  virtual PetscErrorCode write(const char filename[]) const;
  virtual void prepend_history(string message);
protected:
  virtual PetscErrorCode write_attributes(const NCTool &nc, int, nc_type, bool) const;
};

//! An internal class for reading, writing and converting time-series.
class NCTimeseries : public NCVariable {
public:
  string dimension_name;
  void    init(string name, string dim_name, MPI_Comm c, PetscMPIInt r);
  virtual PetscErrorCode read(const char filename[], vector<double> &data);
  virtual PetscErrorCode write(const char filename[], size_t start, vector<double> &data, nc_type nctype = NC_DOUBLE);
  virtual PetscErrorCode write(const char filename[], size_t start, double data, nc_type nctype = NC_DOUBLE);
  virtual PetscErrorCode change_units(vector<double> &data, utUnit *from, utUnit *to);
  virtual PetscErrorCode report_range(vector<double> &data);
protected:
  virtual PetscErrorCode define(int ncid, int dimid, int &varid, nc_type nctype);
  PetscErrorCode define_dimension(int ncid, int &dimid);
};

#endif
