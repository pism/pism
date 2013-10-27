// Copyright (C) 2009--2013 Constantine Khroulev
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

#ifndef __NCVariable_hh
#define __NCVariable_hh

#include <set>
#include <map>
#include <vector>
#include <string>
#include <petscsys.h>
#include "PISMUnits.hh"
#include "PIO.hh"

class PISMTime;

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
  NCVariable(PISMUnitSystem system);
  virtual ~NCVariable() {}
  void init(std::string name, MPI_Comm c, PetscMPIInt r);
  virtual PetscErrorCode set_units(std::string);
  virtual PetscErrorCode set_glaciological_units(std::string);
  virtual PetscErrorCode reset();
  virtual void set(std::string, double);
  virtual double get(std::string) const;
  virtual void set_string(std::string name, std::string value);
  virtual std::string get_string(std::string) const;
  virtual bool has(std::string) const;
  virtual bool is_valid(PetscScalar a) const;
  virtual int get_ndims() const;
  std::string short_name;

  // Attributes:
  /*! Typical attributes stored here:

    \li long_name
    \li standard_name
    \li pism_intent
    \li units
    \li glaciological_units
   */
  std::map<std::string, std::vector<double> > doubles; //!< scalar and array attributes

  virtual PetscErrorCode read_attributes(const PIO &nc);

  virtual PetscErrorCode define(const PIO &nc, PISM_IO_Type nctype,
                                bool write_in_glaciological_units = true) = 0;
protected:
  virtual PetscErrorCode write_attributes(const PIO &nc, PISM_IO_Type nctype,
					  bool write_in_glaciological_units) const;
  virtual PetscErrorCode read_valid_range(const PIO &nc, std::string name);
  MPI_Comm com;
  PetscMPIInt rank;
  std::map<std::string, std::string> strings;  //!< string and boolean attributes
  PISMUnit m_units,		      //!< internal (PISM) units
    m_glaciological_units; //!< \brief for diagnostic variables: units
  //!< to use when writing to a NetCDF file and for standard out reports
  int ndims;
};

//! A class for reading, writing and accessing PISM configuration flags and parameters.
class NCConfigVariable : public NCVariable {
public:
  NCConfigVariable(PISMUnitSystem system);
  ~NCConfigVariable();
  virtual PetscErrorCode print(PetscInt verbosity_threshhold = 4) const;
  virtual PetscErrorCode warn_about_unused_parameters() const;
  virtual PetscErrorCode read(const PIO &nc);
  virtual PetscErrorCode write(const PIO &nc);

  virtual PetscErrorCode read(std::string filename);
  virtual PetscErrorCode write(std::string filename);

  virtual std::string get_config_filename() const;
  virtual PISMUnitSystem get_unit_system() const;
  virtual double get(std::string) const;
  virtual double get(std::string name, std::string u1, std::string u2) const;
  virtual bool   get_flag(std::string) const;
  virtual std::string get_string(std::string name) const;
  // Set a flag (overriding the default in pism_config.cdl). Should not be used
  // in pismr code.
  virtual void   set_flag(std::string, bool);
  // Set parameters and remember that they were set using a command-line option
  virtual PetscErrorCode set_flag_from_option(std::string name, bool value);
  virtual PetscErrorCode set_scalar_from_option(std::string name, double value);
  virtual PetscErrorCode set_string_from_option(std::string name, std::string value);
  virtual PetscErrorCode set_keyword_from_option(std::string name, std::string value);
  // Set parameters by ptocessing a command-line option
  virtual PetscErrorCode flag_from_option(std::string, std::string);
  virtual PetscErrorCode scalar_from_option(std::string, std::string);
  virtual PetscErrorCode string_from_option(std::string, std::string);
  virtual PetscErrorCode keyword_from_option(std::string, std::string, std::string);
  // Import settings from an override file
  virtual void import_from(const NCConfigVariable &other);
  virtual void update_from(const NCConfigVariable &other);

  virtual PetscErrorCode define(const PIO &nc, PISM_IO_Type nctype,
                                bool write_in_glaciological_units = true);
protected:
  std::string config_filename;       //!< \brief the name of the file this config database
                                //!< was initialized from 
  virtual PetscErrorCode write_attributes(const PIO &nc, PISM_IO_Type nctype,
					  bool write_in_glaciological_units) const;

  double get_quiet(std::string) const;
  std::string get_string_quiet(std::string) const;
  bool   get_flag_quiet(std::string) const;

  std::set<std::string> parameters_set;
  mutable std::set<std::string> parameters_used;
  bool options_left_set;
  PISMUnitSystem m_unit_system;
};

//! \brief A class for reading and writing NetCDF global attributes.
/*! This is not a variable, because it has no value, but it is similar to
  NCConfigVariable, because uses of these attributes are similar.
*/
class NCGlobalAttributes : public NCConfigVariable {
public:
  NCGlobalAttributes(PISMUnitSystem system);
  using NCConfigVariable::read;
  using NCConfigVariable::write;
  virtual PetscErrorCode read(const PIO &nc);
  virtual PetscErrorCode write(const PIO &nc);

  virtual void prepend_history(std::string message);
  virtual void set_from_config(const NCConfigVariable &input);
protected:
  virtual PetscErrorCode write_attributes(const PIO &nc, PISM_IO_Type, bool) const;
};

//! An internal class for reading, writing and converting time-series.
class NCTimeseries : public NCVariable {
public:
  NCTimeseries(PISMUnitSystem system);
  std::string dimension_name;        //!< the name of the NetCDF dimension this timeseries depends on
  void    init(std::string name, std::string dim_name, MPI_Comm c, PetscMPIInt r);

  virtual PetscErrorCode read(const PIO &nc, PISMTime *time, std::vector<double> &data);
  virtual PetscErrorCode write(const PIO &nc, size_t start, std::vector<double> &data, PISM_IO_Type nctype = PISM_DOUBLE);
  virtual PetscErrorCode write(const PIO &nc, size_t start, double data, PISM_IO_Type nctype = PISM_DOUBLE);

  virtual PetscErrorCode change_units(std::vector<double> &data, PISMUnit &from, PISMUnit &to);
  virtual PetscErrorCode get_bounds_name(const PIO &nc, std::string &result);
  virtual PetscErrorCode report_range(std::vector<double> &data);

  virtual PetscErrorCode define(const PIO &nc, PISM_IO_Type nctype, bool);
};

class NCTimeBounds : public NCVariable
{
public:
  NCTimeBounds(PISMUnitSystem system);
  void init(std::string var_name, std::string dim_name, MPI_Comm c, PetscMPIInt r);
  virtual PetscErrorCode read(const PIO &nc, PISMTime *time, std::vector<double> &data);
  virtual PetscErrorCode write(const PIO &nc, size_t start, std::vector<double> &data, PISM_IO_Type nctype = PISM_DOUBLE);
  virtual PetscErrorCode write(const PIO &nc, size_t start, double a, double b, PISM_IO_Type nctype = PISM_DOUBLE);

  virtual PetscErrorCode change_units(std::vector<double> &data, PISMUnit &from, PISMUnit &to);

  virtual PetscErrorCode define(const PIO &nc, PISM_IO_Type nctype, bool);
protected:
  std::string dimension_name, bounds_name;
};

#endif  // __NCVariable_hh
