// Copyright (C) 2009--2018, 2020, 2021, 2022, 2023, 2024, 2024, 2025 Constantine Khroulev
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

#ifndef PISM_VARIABLEMETADATA_H
#define PISM_VARIABLEMETADATA_H

#include <map>
#include <vector>
#include <string>
#include <memory>

namespace pism {
namespace units {
class System;
}

namespace io {
enum Type : int;
}

class Grid;
class Logger;

//! @brief A class for handling variable metadata, reading, writing and converting
//! from input units and to output units.
/*! A NetCDF variable can have any number of attributes, but some of them get
  special treatment:

  - units: specifies internal units. When read, a variable is
  converted to these units. When written, it is converted from these
  to output_units.

  - output_units: is never written to a file; replaces 'units'
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
  - units
  - output_units (saved to files as "units")

  Use the `name` of "PISM_GLOBAL" to read and write global attributes.
  (See also File.)

*/

class VariableMetadata;

/*!
 * Syntactic sugar used to make it easier to get attributes.
 *
 * This class makes it possible to get both string and numeric attributes using code that
 * looks like `variable = metadata["attribute"]`. It tries to convert to the type of
 * `variable` and throws an error if there is a type mismatch.
 */
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

/*!
 * Syntactic sugar used to make it easier to set attributes.
 *
 * This class makes it possible to set both string and numeric attributes using code that
 * looks like `metadata["attribute"] = value`.
 */
class Attribute : public ConstAttribute {
public:
  friend class VariableMetadata;
  void operator=(const std::string &value);
  void operator=(const std::initializer_list<double> &value);
  void operator=(const std::vector<double> &value);
private:
  using ConstAttribute::ConstAttribute;
};

class DimensionMetadata;

class VariableAttributes {
public:
  //! string and boolean attributes
  std::map<std::string, std::string> strings;

  //! scalar and array attributes
  std::map<std::string, std::vector<double> > numbers;

  //! @brief The unit system to use.
  std::shared_ptr<units::System> unit_system;

  bool is_set(const std::string &name) const;
};

class VariableMetadata {
public:
  VariableMetadata(const std::string &name, std::shared_ptr<units::System> system,
                   unsigned int ndims = 0);
  virtual ~VariableMetadata() = default;

  Attribute operator[](const std::string &name) {
    return Attribute(this, name);
  }

  ConstAttribute operator[](const std::string &name) const {
    return ConstAttribute(this, name);
  }

  // setters for common attributes

  VariableMetadata &long_name(const std::string &input) {
    return set_string("long_name", input);
  }

  VariableMetadata &standard_name(const std::string &input) {
    return set_string("standard_name", input);
  }

  VariableMetadata &units(const std::string &input) {
    return set_string("units", input);
  }

  VariableMetadata &output_units(const std::string &input) {
    return set_string("output_units", input);
  }

  // getters and setters
  double get_number(const std::string &name) const;
  VariableMetadata &set_number(const std::string &name, double value);

  std::vector<double> get_numbers(const std::string &name) const;
  VariableMetadata &set_numbers(const std::string &name, const std::vector<double> &values);

  std::string get_name() const;
  VariableMetadata &set_name(const std::string &name);

  std::string get_string(const std::string &name) const;
  VariableMetadata &set_string(const std::string &name, const std::string &value);
  VariableMetadata &set_units_without_validation(const std::string &value);

  bool get_time_independent() const;
  VariableMetadata &set_time_independent(bool flag);

  io::Type get_output_type() const;
  VariableMetadata &set_output_type(io::Type type);

  //! Clear all attributes
  VariableMetadata &clear();

  // more getters
  std::shared_ptr<units::System> unit_system() const;

  unsigned int n_spatial_dimensions() const;

  std::vector<DimensionMetadata> dimensions() const;

  bool has_attribute(const std::string &name) const;
  bool has_attributes() const;

  const std::map<std::string, std::string> &all_strings() const;
  const std::map<std::string, std::vector<double> > &all_doubles() const;
  const VariableAttributes &attributes() const;

  void report_to_stdout(const Logger &log, int verbosity_threshold) const;
  void check_range(const std::string &filename, double min, double max) const;
  void report_range(const Logger &log, double min, double max) const;

protected:
  unsigned int m_n_spatial_dims;

  virtual std::vector<DimensionMetadata> dimensions_impl() const;

private:
  VariableAttributes m_attributes;

  std::string m_short_name;
  bool m_time_independent;

  io::Type m_output_type;
};

class DimensionMetadata : public VariableMetadata {
public:
  DimensionMetadata(const std::string &name, std::shared_ptr<units::System> system,
                    int length);
  int length() const;
private:
  std::vector<DimensionMetadata> dimensions_impl() const;

  int m_length;
};

class OtherMetadata : public VariableMetadata {
public:
  OtherMetadata(const std::string &name,
                const std::vector<DimensionMetadata> &dimensions,
                std::shared_ptr<units::System> system);
  OtherMetadata(const std::string &name,
                const std::vector<std::tuple<std::string, int>> &dimensions,
                std::shared_ptr<units::System> system);
private:
  std::vector<DimensionMetadata> m_dimensions;
  std::vector<DimensionMetadata> dimensions_impl() const;
};

//! Spatial NetCDF variable (corresponding to a 2D or 3D scalar field).
class SpatialVariableMetadata : public VariableMetadata {
public:
  SpatialVariableMetadata(std::shared_ptr<units::System> system, const std::string &name,
                          const Grid &grid,
                          const std::vector<double> &levels = { 0.0 });
  SpatialVariableMetadata(std::shared_ptr<units::System> system, const std::string &name,
                          unsigned int Mx, double dx, unsigned int My, double dy,
                          const std::vector<double> &levels = { 0.0 });
  virtual ~SpatialVariableMetadata() = default;

  const std::vector<double>& levels() const;

  DimensionMetadata& x();
  DimensionMetadata& y();
  DimensionMetadata& z();

  const DimensionMetadata& x() const;
  const DimensionMetadata& y() const;
  const DimensionMetadata& z() const;

private:
  std::vector<DimensionMetadata> dimensions_impl() const;
  DimensionMetadata m_x, m_y, m_z;
  std::vector<double> m_zlevels;
};

// Comparison operator for VariableMetadata (we need it to store VariableMetadata in
// sorted containers)
inline bool operator<(const VariableMetadata &a, const VariableMetadata &b) {
  return a.get_name() < b.get_name();
}

} // end of namespace pism

#endif  // PISM_VARIABLEMETADATA_H
