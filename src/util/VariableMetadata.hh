// Copyright (C) 2009--2018 Constantine Khroulev
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

#ifndef __VariableMetadata_hh
#define __VariableMetadata_hh

#include <set>
#include <map>
#include <vector>
#include <string>

#include "Units.hh"

#include "pism/util/io/IO_Flags.hh" // IO_Type

namespace pism {

class Logger;

//! @brief A class for handling variable metadata, reading, writing and converting
//! from input units and to output units.
/*! A NetCDF variable can have any number of attributes, but some of them get
  special treatment:

  - units: specifies internal units. When read, a variable is
  converted to these units. When written, it is converted from these
  to glaciological_units.

  - glaciological_units: is never written to a file; replaces 'units'
  in the output file.

  - valid_min, valid_max: specify the valid range of a variable. Are
  read from an input file *only* if not specified previously. If
  both are set, then valid_range is used in the output instead.

  Also:

  - empty string attributes are ignored (they are not written to the
  output file and has_attribute("foo") returns false if "foo" is
  absent or equal to an empty string).

  Typical attributes stored here:

  - long_name
  - standard_name
  - pism_intent
  - units
  - glaciological_units (saved to files as "units")

  Use the `name` of "PISM_GLOBAL" to read and write global attributes.
  (See also PIO.)

*/

class VariableMetadata {
public:
  VariableMetadata(const std::string &name, units::System::Ptr system, unsigned int ndims = 0);
  virtual ~VariableMetadata();

  // setters
  void set_double(const std::string &name, double value);
  void set_doubles(const std::string &name, const std::vector<double> &values);
  void set_name(const std::string &name);
  void set_string(const std::string &name, const std::string &value);

  void set_time_independent(bool flag);
  void set_output_type(IO_Type type);

  void clear_all_doubles();
  void clear_all_strings();

  // getters
  units::System::Ptr unit_system() const;

  double get_double(const std::string &name) const;
  std::vector<double> get_doubles(const std::string &name) const;
  std::string get_name() const;
  std::string get_string(const std::string &name) const;

  unsigned int get_n_spatial_dimensions() const;

  bool has_attribute(const std::string &name) const;
  bool has_attributes() const;
  bool get_time_independent() const;
  IO_Type get_output_type() const;

  typedef std::map<std::string,std::string> StringAttrs;
  const StringAttrs& get_all_strings() const;

  typedef std::map<std::string,std::vector<double> > DoubleAttrs;
  const DoubleAttrs& get_all_doubles() const;

  void report_to_stdout(const Logger &log, int verbosity_threshold) const;
  void check_range(const std::string &filename, double min, double max);
  void report_range(const Logger &log, double min, double max, bool found_by_standard_name);

protected:
  unsigned int m_n_spatial_dims;

private:
  //! @brief The unit system to use.
  units::System::Ptr m_unit_system;

  //! string and boolean attributes
  std::map<std::string, std::string> m_strings;

  //! scalar and array attributes
  std::map<std::string, std::vector<double> > m_doubles;
  std::string m_short_name;
  bool m_time_independent;

  IO_Type m_output_type;
};

bool set_contains(const std::set<std::string> &S, const VariableMetadata &field);

//! Spatial NetCDF variable (corresponding to a 2D or 3D scalar field).
class SpatialVariableMetadata : public VariableMetadata {
public:
  SpatialVariableMetadata(units::System::Ptr system, const std::string &name);
  SpatialVariableMetadata(units::System::Ptr system, const std::string &name,
                          const std::vector<double> &zlevels);
  SpatialVariableMetadata(const SpatialVariableMetadata &other);
  virtual ~SpatialVariableMetadata();

  void set_levels(const std::vector<double> &levels);
  const std::vector<double>& get_levels() const;

  VariableMetadata& get_x();
  VariableMetadata& get_y();
  VariableMetadata& get_z();

  const VariableMetadata& get_x() const;
  const VariableMetadata& get_y() const;
  const VariableMetadata& get_z() const;

private:
  VariableMetadata m_x, m_y, m_z;
  std::vector<double> m_zlevels;
  void init_internal(const std::string &name,
                     const std::vector<double> &z_levels);
};

//! An internal class for reading, writing and converting time-series.
class TimeseriesMetadata : public VariableMetadata {
public:
  TimeseriesMetadata(const std::string &name, const std::string &dimension_name,
                     units::System::Ptr system);
  virtual ~TimeseriesMetadata();

  std::string get_dimension_name() const;
private:
  //! the name of the NetCDF dimension this timeseries depends on
  std::string m_dimension_name;
};

class TimeBoundsMetadata : public TimeseriesMetadata
{
public:
  TimeBoundsMetadata(const std::string &name, const std::string &dimension_name,
                     units::System::Ptr system);
  virtual ~TimeBoundsMetadata();
  std::string get_bounds_name() const;
private:
  std::string m_bounds_name;
};

} // end of namespace pism

#endif  // __VariableMetadata_hh
