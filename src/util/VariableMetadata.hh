// Copyright (C) 2009--2018, 2020, 2021, 2022, 2023 Constantine Khroulev
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

#include "pism/util/Units.hh"

namespace pism {
namespace io {
enum Type : int;
}

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
  (See also File.)

*/

class VariableMetadata;

class ConstAttribute {
public:
  friend class VariableMetadata;
  ConstAttribute(const ConstAttribute&) = delete;
  ConstAttribute& operator=(const ConstAttribute&) = delete;

  operator std::string() const;
  operator double() const;
  operator std::vector<double> () const;
protected:
  ConstAttribute(const VariableMetadata *var, const std::string &name);
  ConstAttribute(ConstAttribute&& a) noexcept;

  std::string m_name;
  VariableMetadata* m_var;
};

class Attribute : public ConstAttribute {
public:
  friend class VariableMetadata;
  void operator=(const std::string &value);
  void operator=(const std::initializer_list<double> &value);
  void operator=(const std::vector<double> &value);
private:
  using ConstAttribute::ConstAttribute;
};

class VariableMetadata {
public:
  VariableMetadata(const std::string &name, units::System::Ptr system, unsigned int ndims = 0);
  virtual ~VariableMetadata() = default;

  Attribute operator[](const std::string &name) {
    return Attribute(this, name);
  }

  const ConstAttribute operator[](const std::string &name) const {
    return ConstAttribute(this, name);
  }

  // getters and setters
  double get_number(const std::string &name) const;
  void set_number(const std::string &name, double value);

  std::vector<double> get_numbers(const std::string &name) const;
  void set_numbers(const std::string &name, const std::vector<double> &values);

  std::string get_name() const;
  void set_name(const std::string &name);

  std::string get_string(const std::string &name) const;
  void set_string(const std::string &name, const std::string &value);
  void set_units_without_validation(const std::string &value);

  bool get_time_independent() const;
  void set_time_independent(bool flag);

  io::Type get_output_type() const;
  void set_output_type(io::Type type);

  void clear_all_doubles();
  void clear_all_strings();

  // more getters
  units::System::Ptr unit_system() const;

  unsigned int n_spatial_dimensions() const;

  bool has_attribute(const std::string &name) const;
  bool has_attributes() const;

  typedef std::map<std::string,std::string> StringAttrs;
  const StringAttrs& all_strings() const;

  typedef std::map<std::string,std::vector<double> > DoubleAttrs;
  const DoubleAttrs& all_doubles() const;

  void report_to_stdout(const Logger &log, int verbosity_threshold) const;
  void check_range(const std::string &filename, double min, double max) const;
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

  io::Type m_output_type;
};

//! Spatial NetCDF variable (corresponding to a 2D or 3D scalar field).
class SpatialVariableMetadata : public VariableMetadata {
public:
  SpatialVariableMetadata(units::System::Ptr system, const std::string &name,
                          const std::vector<double> &zlevels = {0.0});
  virtual ~SpatialVariableMetadata() = default;

  const std::vector<double>& levels() const;

  VariableMetadata& x();
  VariableMetadata& y();
  VariableMetadata& z();

  const VariableMetadata& x() const;
  const VariableMetadata& y() const;
  const VariableMetadata& z() const;

private:
  VariableMetadata m_x, m_y, m_z;
  std::vector<double> m_zlevels;
};

} // end of namespace pism

#endif  // __VariableMetadata_hh
