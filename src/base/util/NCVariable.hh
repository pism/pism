// Copyright (C) 2009--2012 Constantine Khroulev
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

#include <set>
#include <map>
#include <vector>
#include <string>
#include <petscsys.h>
#include "udunits.h"
#include "PIO.hh"

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond


//! \brief A class for handling variable metadata, reading, writing and converting
//! from input units and to output units.
/*! A NetCDF variable can have any number of attributes, but some of them get
    special treatment:

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
  virtual int get_ndims() const;
  string short_name;

  // Attributes:
  /*! Typical attributes stored here:

    \li long_name
    \li standard_name
    \li pism_intent
    \li units
    \li glaciological_units
   */
  map<string, vector<double> > doubles; //!< scalar and array attributes

  virtual PetscErrorCode read_attributes(const PIO &nc);

  virtual PetscErrorCode define(const PIO &nc, PISM_IO_Type nctype,
                                bool write_in_glaciological_units = true) = 0;
protected:
  virtual PetscErrorCode write_attributes(const PIO &nc, PISM_IO_Type nctype,
					  bool write_in_glaciological_units) const;
  virtual PetscErrorCode read_valid_range(const PIO &nc, string name);
  MPI_Comm com;
  PetscMPIInt rank;
  map<string, string> strings;  //!< string and boolean attributes
  utUnit units,		      //!< internal (PISM) units
         glaciological_units; //!< \brief for diagnostic variables: units to
			      //!< use when writing to a NetCDF file and for
			      //!< standard out reports
  int ndims;
};

//! A class for reading, writing and accessing PISM configuration flags and parameters.
class NCConfigVariable : public NCVariable {
public:
  NCConfigVariable();
  ~NCConfigVariable();
  virtual PetscErrorCode print(PetscInt verbosity_threshhold) const;
  virtual PetscErrorCode print() const { print(4); return 0; };
  virtual PetscErrorCode warn_about_unused_parameters() const;
  virtual PetscErrorCode read(const PIO &nc);
  virtual PetscErrorCode write(const PIO &nc);

  virtual PetscErrorCode read(string filename);
  virtual PetscErrorCode write(string filename);

  virtual string get_config_filename() const;
  virtual double get(string) const;
  virtual double get(string name, string u1, string u2) const;
  virtual bool   get_flag(string) const;
  virtual string get_string(string name) const;
  // Set a flag (overriding the default in pism_config.cdl). Should not be used
  // in pismr code.
  virtual void   set_flag(string, bool);
  // Set parameters and remember that they were set using a command-line option
  virtual PetscErrorCode set_flag_from_option(string name, bool value);
  virtual PetscErrorCode set_scalar_from_option(string name, double value);
  virtual PetscErrorCode set_string_from_option(string name, string value);
  virtual PetscErrorCode set_keyword_from_option(string name, string value);
  // Set parameters by ptocessing a command-line option
  virtual PetscErrorCode flag_from_option(string, string);
  virtual PetscErrorCode scalar_from_option(string, string);
  virtual PetscErrorCode string_from_option(string, string);
  virtual PetscErrorCode keyword_from_option(string, string, string);
  // Import settings from an override file
  virtual void import_from(const NCConfigVariable &other);
  virtual void update_from(const NCConfigVariable &other);

  virtual PetscErrorCode define(const PIO &nc, PISM_IO_Type nctype,
                                bool write_in_glaciological_units = true);
protected:
  string config_filename;       //!< \brief the name of the file this config database
                                //!< was initialized from 
  virtual PetscErrorCode write_attributes(const PIO &nc, PISM_IO_Type nctype,
					  bool write_in_glaciological_units) const;

  double get_quiet(string) const;
  string get_string_quiet(string) const;
  bool   get_flag_quiet(string) const;

  std::set<string> parameters_set;
  mutable std::set<string> parameters_used;
  bool options_left_set;
};

//! \brief A class for reading and writing NetCDF global attributes.
/*! This is not a variable, because it has no value, but it is similar to
  NCConfigVariable, because uses of these attributes are similar.
*/
class NCGlobalAttributes : public NCConfigVariable {
public:
  using NCConfigVariable::read;
  using NCConfigVariable::write;
  virtual PetscErrorCode read(const PIO &nc);
  virtual PetscErrorCode write(const PIO &nc);

  virtual void prepend_history(string message);
  virtual void set_from_config(const NCConfigVariable &input);
protected:
  virtual PetscErrorCode write_attributes(const PIO &nc, PISM_IO_Type, bool) const;
};

//! An internal class for reading, writing and converting time-series.
class NCTimeseries : public NCVariable {
public:
  string dimension_name;        //!< the name of the NetCDF dimension this timeseries depends on
  void    init(string name, string dim_name, MPI_Comm c, PetscMPIInt r);

  virtual PetscErrorCode read(const PIO &nc, bool use_reference_date, vector<double> &data);
  virtual PetscErrorCode write(const PIO &nc, size_t start, vector<double> &data, PISM_IO_Type nctype = PISM_DOUBLE);
  virtual PetscErrorCode write(const PIO &nc, size_t start, double data, PISM_IO_Type nctype = PISM_DOUBLE);

  // virtual PetscErrorCode read(string filename, bool use_reference_date, vector<double> &data);
  // virtual PetscErrorCode write(string filename, size_t start, vector<double> &data, PISM_IO_Type nctype = PISM_DOUBLE);
  // virtual PetscErrorCode write(string filename, size_t start, double data, PISM_IO_Type nctype = PISM_DOUBLE);

  virtual PetscErrorCode change_units(vector<double> &data, utUnit *from, utUnit *to);
  virtual PetscErrorCode get_bounds_name(const PIO &nc, string &result);
  virtual PetscErrorCode report_range(vector<double> &data);

  virtual PetscErrorCode define(const PIO &nc, PISM_IO_Type nctype, bool);
};

class NCTimeBounds : public NCVariable
{
public:
  void init(string var_name, string dim_name, MPI_Comm c, PetscMPIInt r);
  virtual PetscErrorCode read(const PIO &nc, bool use_reference_date, vector<double> &data);
  virtual PetscErrorCode write(const PIO &nc, size_t start, vector<double> &data, PISM_IO_Type nctype = PISM_DOUBLE);
  virtual PetscErrorCode write(const PIO &nc, size_t start, double a, double b, PISM_IO_Type nctype = PISM_DOUBLE);

  // virtual PetscErrorCode read(string filename, bool use_reference_date, vector<double> &data);
  // virtual PetscErrorCode write(string filename, size_t start, vector<double> &data, PISM_IO_Type nctype = PISM_DOUBLE);
  // virtual PetscErrorCode write(string filename, size_t start, double a, double b, PISM_IO_Type nctype = PISM_DOUBLE);

  virtual PetscErrorCode change_units(vector<double> &data, utUnit *from, utUnit *to);

  virtual PetscErrorCode define(const PIO &nc, PISM_IO_Type nctype, bool);
protected:
  string dimension_name, bounds_name;
};

#endif  // __NCVariable_hh
