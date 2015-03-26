// Copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015 Constantine Khroulev and Ed Bueler
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

#include <sstream>
#include <set>
#include <cassert>
#include <algorithm>

#include "VariableMetadata.hh"
#include "PIO.hh"
#include "pism_options.hh"
#include "IceGrid.hh"

#include "PISMConfigInterface.hh"
#include "error_handling.hh"

namespace pism {

static void compute_range(MPI_Comm com, double *data, size_t data_size, double *min, double *max) {
  double
    min_result = data[0],
    max_result = data[0];
  for (size_t k = 0; k < data_size; ++k) {
    min_result = std::min(min_result, data[k]);
    max_result = std::max(max_result, data[k]);
  }

  if (min) {
    *min = GlobalMin(com, min_result);
  }

  if (max) {
    *max = GlobalMax(com, max_result);
  }
}


VariableMetadata::VariableMetadata(const std::string &name, const UnitSystem &system, unsigned int ndims)
  : m_n_spatial_dims(ndims),
    m_unit_system(system),
    m_short_name(name) {

  clear_all_strings();
  clear_all_doubles();

  // long_name is unset
  // standard_name is unset
  // pism_intent is unset
  // coordinates is unset

  // valid_min and valid_max are unset
}

VariableMetadata::~VariableMetadata() {
  // empty
}

/** Get the number of spatial dimensions.
 */
unsigned int VariableMetadata::get_n_spatial_dimensions() const {
  return m_n_spatial_dims;
}

UnitSystem VariableMetadata::unit_system() const {
  return m_unit_system;
}

//! @brief Check if the range `[min, max]` is a subset of `[valid_min, valid_max]`.
/*! Throws an exception if this check failed.
 */
void VariableMetadata::check_range(const std::string &filename, double min, double max) {

  const std::string &units_string = get_string("units");
  const char
    *units = units_string.c_str(),
    *name = get_name().c_str(),
    *file = filename.c_str();

  if (has_attribute("valid_min") and has_attribute("valid_max")) {
    double
      valid_min = get_double("valid_min"),
      valid_max = get_double("valid_max");
    if ((min < valid_min) or (max > valid_max)) {
      throw RuntimeError::formatted("some values of '%s' in '%s' are outside the valid range [%e, %e] (%s).\n"
                                    "computed min = %e %s, computed max = %e %s",
                                    name, file,
                                    valid_min, valid_max, units, min, units, max, units);
    }
  } else if (has_attribute("valid_min")) {
    double valid_min = get_double("valid_min");
    if (min < valid_min) {
      throw RuntimeError::formatted("some values of '%s' in '%s' are less than the valid minimum (%e %s).\n"
                                    "computed min = %e %s, computed max = %e %s",
                                    name, file,
                                    valid_min, units, min, units, max, units);
    }
  } else if (has_attribute("valid_max")) {
    double valid_max = get_double("valid_max");
    if (max > valid_max) {
      throw RuntimeError::formatted("some values of '%s' in '%s' are greater than the valid maximum (%e %s).\n"
                                    "computed min = %e %s, computed max = %e %s",
                                    name, file,
                                    valid_max, units, min, units, max, units);
    }
  }
}

//! 3D version
SpatialVariableMetadata::SpatialVariableMetadata(const UnitSystem &system, const std::string &name,
                                     const IceGrid &g, const std::vector<double> &zlevels)
  : VariableMetadata("unnamed", system),
    m_x("x", system),
    m_y("y", system),
    m_z("z", system) {

  init_internal(name, g, zlevels);
}

//! 2D version
SpatialVariableMetadata::SpatialVariableMetadata(const UnitSystem &system, const std::string &name,
                                     const IceGrid &g)
  : VariableMetadata("unnamed", system),
    m_x("x", system),
    m_y("y", system),
    m_z("z", system) {

  std::vector<double> z(1, 0.0);
  init_internal(name, g, z);
}

void SpatialVariableMetadata::init_internal(const std::string &name, const IceGrid &g,
                                      const std::vector<double> &z_levels) {
  m_time_dimension_name = "t";        // will be overridden later

  m_x.set_string("axis", "X");
  m_x.set_string("long_name", "X-coordinate in Cartesian system");
  m_x.set_string("standard_name", "projection_x_coordinate");
  m_x.set_string("units", "m");

  m_y.set_string("axis", "Y");
  m_y.set_string("long_name", "Y-coordinate in Cartesian system");
  m_y.set_string("standard_name", "projection_y_coordinate");
  m_y.set_string("units", "m");

  m_z.set_string("axis", "Z");
  m_z.set_string("long_name", "Z-coordinate in Cartesian system");
  m_z.set_string("units", "m");
  m_z.set_string("positive", "up");

  set_string("coordinates", "lat lon");
  set_string("grid_mapping", "mapping");

  m_variable_order = "xyz";

  set_name(name);
  m_grid = &g;

  m_zlevels = z_levels;

  this->set_time_independent(false);

  if (m_zlevels.size() > 1) {
    get_z().set_name("z");      // default; can be overridden easily
    m_n_spatial_dims = 3;
  } else {
    get_z().set_name("");
    m_n_spatial_dims = 2;
  }

  m_variable_order = m_grid->config.get_string("output_variable_order");
}

SpatialVariableMetadata::SpatialVariableMetadata(const SpatialVariableMetadata &other)
  : VariableMetadata(other), m_x(other.m_x), m_y(other.m_y), m_z(other.m_z) {
  m_time_dimension_name = other.m_time_dimension_name;
  m_variable_order      = other.m_variable_order;
  m_zlevels             = other.m_zlevels;
  m_grid                = other.m_grid;
}

SpatialVariableMetadata::~SpatialVariableMetadata() {
  // empty
}

void SpatialVariableMetadata::set_levels(const std::vector<double> &levels) {
  if (levels.size() < 1) {
    throw RuntimeError("argument \"levels\" has to have length 1 or greater");
  }
  m_zlevels = levels;
}

const std::vector<double>& SpatialVariableMetadata::get_levels() const {
  return m_zlevels;
}

void SpatialVariableMetadata::set_time_independent(bool flag) {
  if (flag == true) {
    m_time_dimension_name = "";
  } else {
    m_time_dimension_name = m_grid->config.get_string("time_dimension_name");
  }
}

const IceGrid& SpatialVariableMetadata::grid() const {
  return *m_grid;
}


//! Read a variable from a file into an array `output`.
/*! This also converts the data from input units to internal units if needed.
 */
void SpatialVariableMetadata::read(const PIO &nc, unsigned int time, double *output) {

  assert(m_grid != NULL);

  // Find the variable:
  std::string name_found;
  bool found_by_standard_name = false, variable_exists = false;
  nc.inq_var(get_name(), get_string("standard_name"),
             variable_exists, name_found, found_by_standard_name);

  if (not variable_exists) {
    throw RuntimeError::formatted("Can't find '%s' (%s) in '%s'.",
                                  get_name().c_str(),
                                  get_string("standard_name").c_str(), nc.inq_filename().c_str());
  }

  // Sanity check: the variable in an input file should have the expected
  // number of spatial dimensions.
  {
    // Set of spatial dimensions this field has.
    std::set<int> axes;
    axes.insert(X_AXIS);
    axes.insert(Y_AXIS);
    if (get_z().get_name().empty() == false) {
      axes.insert(Z_AXIS);
    }

    std::vector<std::string> input_dims;
    int input_ndims = 0;                 // number of spatial dimensions (input file)
    size_t matching_dim_count = 0; // number of matching dimensions

    input_dims = nc.inq_vardims(name_found);
    std::vector<std::string>::iterator j = input_dims.begin();
    while (j != input_dims.end()) {
      AxisType tmp = nc.inq_dimtype(*j);

      if (tmp != T_AXIS) {
        ++input_ndims;
      }

      if (axes.find(tmp) != axes.end()) {
        ++matching_dim_count;
      }

      ++j;
    }


    if (axes.size() != matching_dim_count) {

      // Print the error message and stop:
      throw RuntimeError::formatted("found the %dD variable %s (%s) in '%s' while trying to read\n"
                                    "'%s' ('%s'), which is %d-dimensional.",
                                    input_ndims, name_found.c_str(),
                                    join(input_dims, ",").c_str(),
                                    nc.inq_filename().c_str(),
                                    get_name().c_str(), get_string("long_name").c_str(),
                                    static_cast<int>(axes.size()));
    }
  }

  // make sure we have at least one level
  unsigned int nlevels = std::max(m_zlevels.size(), (size_t)1);

  nc.get_vec(*m_grid, name_found, nlevels, time, output);

  std::string input_units = nc.get_att_text(name_found, "units");

  if (has_attribute("units") && input_units.empty()) {
    const std::string &units_string = get_string("units"),
      &long_name = get_string("long_name");
    verbPrintf(2, nc.com(),
               "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
               "              Assuming that it is in '%s'.\n",
               get_name().c_str(), long_name.c_str(),
               units_string.c_str());
    input_units = units_string;
  }

  // Convert data:
  size_t size = m_grid->xm() * m_grid->ym() * nlevels;

  UnitConverter(m_unit_system,
                input_units,
                get_string("units")).convert_doubles(output, size);
}

//! \brief Write a \b global Vec `v` to a variable.
/*!
  Defines a variable and converts the units if needed.
 */
void SpatialVariableMetadata::write(const PIO &nc, IO_Type nctype,
                              bool write_in_glaciological_units,
                              const double *input) const {

  // find or define the variable
  std::string name_found;
  bool exists, found_by_standard_name;
  nc.inq_var(get_name(), get_string("standard_name"),
             exists, name_found, found_by_standard_name);

  if (not exists) {
    define(nc, nctype, write_in_glaciological_units);
    name_found = get_name();
  }

  // make sure we have at least one level
  unsigned int nlevels = std::max(m_zlevels.size(), (size_t)1);

  if (write_in_glaciological_units) {
    size_t data_size = m_grid->xm() * m_grid->ym() * nlevels;

    // create a temporary array, convert to glaciological units, and
    // save
    std::vector<double> tmp(data_size);
    for (size_t k = 0; k < data_size; ++k) {
      tmp[k] = input[k];
    }

    UnitConverter(m_unit_system,
                  get_string("units"),
                  get_string("glaciological_units")).convert_doubles(&tmp[0], tmp.size());

    nc.put_vec(*m_grid, name_found, nlevels, &tmp[0]);
  } else {
    nc.put_vec(*m_grid, name_found, nlevels, input);
  }
}

//! \brief Regrid from a NetCDF file into a distributed array `output`.
/*!
  - if `flag` is `CRITICAL` or `CRITICAL_FILL_MISSING`, stops if the
    variable was not found in the input file
  - if `flag` is one of `CRITICAL_FILL_MISSING` and
    `OPTIONAL_FILL_MISSING`, replace _FillValue with `default_value`.
  - sets `v` to `default_value` if `flag` is `OPTIONAL` and the
    variable was not found in the input file
  - uses the last record in the file
 */
void SpatialVariableMetadata::regrid(const PIO &nc, RegriddingFlag flag, bool do_report_range,
                               double default_value, double *output) {
  unsigned int t_length = nc.inq_nrecords(get_name(), get_string("standard_name"));

  this->regrid(nc, t_length - 1, flag, do_report_range, default_value, output);
}

void SpatialVariableMetadata::regrid(const PIO &nc, unsigned int t_start,
                               RegriddingFlag flag, bool do_report_range,
                               double default_value, double *output) {
  assert(m_grid != NULL);

  const size_t data_size = m_grid->xm() * m_grid->ym() * m_zlevels.size();

  // Find the variable
  bool exists, found_by_standard_name;
  std::string name_found;
  nc.inq_var(get_name(), get_string("standard_name"),
             exists, name_found, found_by_standard_name);

  if (exists) {                      // the variable was found successfully

    if (flag == OPTIONAL_FILL_MISSING or flag == CRITICAL_FILL_MISSING) {
      verbPrintf(2, nc.com(),
                 "PISM WARNING: Replacing missing values with %f [%s] in variable '%s' read from '%s'.\n",
                 default_value, get_string("units").c_str(), get_name().c_str(),
                 nc.inq_filename().c_str());

      nc.regrid_vec_fill_missing(*m_grid, name_found, m_zlevels,
                                 t_start, default_value, output);
    } else {
      nc.regrid_vec(*m_grid, name_found, m_zlevels, t_start, output);
    }

    // Now we need to get the units string from the file and convert
    // the units, because check_range and report_range expect data to
    // be in PISM (MKS) units.

    std::string input_units = nc.get_att_text(name_found, "units");

    if (input_units.empty()) {
      std::string internal_units = get_string("units");
      input_units = internal_units;
      if (not internal_units.empty()) {
        verbPrintf(2, nc.com(),
                   "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
                   "              Assuming that it is in '%s'.\n",
                   get_name().c_str(),
                   get_string("long_name").c_str(),
                   internal_units.c_str());
      }
    }

    // Convert data:
    UnitConverter(m_unit_system,
                  input_units,
                  get_string("units")).convert_doubles(output, data_size);

    // Read the valid range info:
    nc.read_valid_range(name_found, *this);

    double min = 0.0, max = 0.0;
    compute_range(nc.com(), output, data_size, &min, &max);
    
    // Check the range and warn the user if needed:
    check_range(nc.inq_filename(), min, max);

    if (do_report_range) {
      // We can report the success, and the range now:
      verbPrintf(2, nc.com(), "  FOUND ");

      this->report_range(nc.com(), min, max, found_by_standard_name);
    }
  } else {                // couldn't find the variable
    if (flag == CRITICAL or flag == CRITICAL_FILL_MISSING) {
      // if it's critical, print an error message and stop
      throw RuntimeError::formatted("Can't find '%s' in the regridding file '%s'.",
                                    get_name().c_str(), nc.inq_filename().c_str());
    }

    // If it is optional, fill with the provided default value.
    // UnitConverter constructor will make sure that units are compatible.
    UnitConverter c(m_unit_system,
                    this->get_string("units"),
                    this->get_string("glaciological_units"));

    std::string spacer(get_name().size(), ' ');
    verbPrintf(2, nc.com(),
               "  absent %s / %-10s\n"
               "         %s \\ not found; using default constant %7.2f (%s)\n",
               get_name().c_str(),
               get_string("long_name").c_str(),
               spacer.c_str(), c(default_value),
               get_string("glaciological_units").c_str());

    for (size_t k = 0; k < data_size; ++k) {
      output[k] = default_value;
    }
  } // end of if (exists)
}


//! Report the range of a \b global Vec `v`.
void SpatialVariableMetadata::report_range(MPI_Comm com, double min, double max,
                                     bool found_by_standard_name) {

  // UnitConverter constructor will make sure that units are compatible.
  UnitConverter c(m_unit_system,
                  this->get_string("units"),
                  this->get_string("glaciological_units"));
  min = c(min);
  max = c(max);

  std::string spacer(get_name().size(), ' ');

  if (has_attribute("standard_name")) {

    if (found_by_standard_name) {
      verbPrintf(2, com,
                 " %s / standard_name=%-10s\n"
                 "         %s \\ min,max = %9.3f,%9.3f (%s)\n",
                 get_name().c_str(),
                 get_string("standard_name").c_str(), spacer.c_str(), min, max,
                 get_string("glaciological_units").c_str());
    } else {
      verbPrintf(2, com,
                 " %s / WARNING! standard_name=%s is missing, found by short_name\n"
                 "         %s \\ min,max = %9.3f,%9.3f (%s)\n",
                 get_name().c_str(),
                 get_string("standard_name").c_str(), spacer.c_str(), min, max,
                 get_string("glaciological_units").c_str());
    }

  } else {
    verbPrintf(2, com,
               " %s / %-10s\n"
               "         %s \\ min,max = %9.3f,%9.3f (%s)\n",
               get_name().c_str(),
               get_string("long_name").c_str(), spacer.c_str(), min, max,
               get_string("glaciological_units").c_str());
  }
}

//! \brief Define dimensions a variable depends on.
void SpatialVariableMetadata::define_dimensions(const PIO &nc) const {
  bool exists;

  // x
  exists = nc.inq_dim(get_x().get_name());
  if (not exists) {
    nc.def_dim(m_grid->Mx(), m_x);
    nc.put_dim(get_x().get_name(), m_grid->x());
  }

  // y
  exists = nc.inq_dim(get_y().get_name());
  if (not exists) {
    nc.def_dim(m_grid->My(), m_y);
    nc.put_dim(get_y().get_name(), m_grid->y());
  }

  // z
  std::string z_name = get_z().get_name();
  if (z_name.empty() == false) {
    exists = nc.inq_dim(z_name);
    if (not exists) {
      unsigned int nlevels = std::max(m_zlevels.size(), (size_t)1); // make sure we have at least one level
      nc.def_dim(nlevels, m_z);
      nc.put_dim(z_name, m_zlevels);
    }
  }
}

//! Define a NetCDF variable corresponding to a VariableMetadata object.
void SpatialVariableMetadata::define(const PIO &nc, IO_Type nctype,
                               bool write_in_glaciological_units) const {
  std::vector<std::string> dims;

  bool exists = nc.inq_var(get_name());
  if (exists) {
    return;
  }

  define_dimensions(nc);
  std::string variable_order = m_variable_order;
  // "..._bounds" should be stored with grid corners (corresponding to
  // the "z" dimension here) last, so we override the variable storage
  // order here
  if (ends_with(get_name(), "_bounds") and variable_order == "zyx") {
    variable_order = "yxz";
  }

  std::string
    x = get_x().get_name(),
    y = get_y().get_name(),
    z = get_z().get_name(),
    t = m_time_dimension_name;

  nc.redef();

  if (not t.empty()) {
    dims.push_back(t);
  }

  // Use t,x,y,z(zb) variable order: it is weird, but matches the in-memory
  // storage order and so is *a lot* faster.
  if (variable_order == "xyz") {
    dims.push_back(x);
    dims.push_back(y);
  }

  // Use the t,y,x,z variable order: also weird, somewhat slower, but 2D fields
  // are stored in the "natural" order.
  if (variable_order == "yxz") {
    dims.push_back(y);
    dims.push_back(x);
  }

  if (not z.empty()) {
    dims.push_back(z);
  }

  // Use the t,z(zb),y,x variables order: more natural for plotting and post-processing,
  // but requires transposing data while writing and is *a lot* slower.
  if (variable_order == "zyx") {
    dims.push_back(y);
    dims.push_back(x);
  }

  nc.def_var(get_name(), nctype, dims);

  nc.write_attributes(*this, nctype, write_in_glaciological_units);
}

VariableMetadata& SpatialVariableMetadata::get_x() {
  return m_x;
}

VariableMetadata& SpatialVariableMetadata::get_y() {
  return m_y;
}

VariableMetadata& SpatialVariableMetadata::get_z() {
  return m_z;
}

const VariableMetadata& SpatialVariableMetadata::get_x() const {
  return m_x;
}

const VariableMetadata& SpatialVariableMetadata::get_y() const {
  return m_y;
}

const VariableMetadata& SpatialVariableMetadata::get_z() const {
  return m_z;
}

//! Checks if an attribute is present. Ignores empty strings, except
//! for the "units" attribute.
bool VariableMetadata::has_attribute(const std::string &name) const {

  std::map<std::string,std::string>::const_iterator j = m_strings.find(name);
  if (j != m_strings.end()) {
    if (name != "units" && (j->second).empty()) {
      return false;
    }

    return true;
  }

  if (m_doubles.find(name) != m_doubles.end()) {
    return true;
  }

  return false;
}

void VariableMetadata::set_name(const std::string &name) {
  m_short_name = name;
}

//! Set a scalar attribute to a single (scalar) value.
void VariableMetadata::set_double(const std::string &name, double value) {
  m_doubles[name] = std::vector<double>(1, value);
}

//! Set a scalar attribute to a single (scalar) value.
void VariableMetadata::set_doubles(const std::string &name, const std::vector<double> &values) {
  m_doubles[name] = values;
}

void VariableMetadata::clear_all_doubles() {
  m_doubles.clear();
}

void VariableMetadata::clear_all_strings() {
  m_strings.clear();
}

std::string VariableMetadata::get_name() const {
  return m_short_name;
}

//! Get a single-valued scalar attribute.
double VariableMetadata::get_double(const std::string &name) const {
  std::map<std::string,std::vector<double> >::const_iterator j = m_doubles.find(name);
  if (j != m_doubles.end()) {
    return (j->second)[0];
  } else {
    throw RuntimeError::formatted("variable \"%s\" does not have a double attribute \"%s\"",
                                  get_name().c_str(), name.c_str());
  }
}

//! Get an array-of-doubles attribute.
std::vector<double> VariableMetadata::get_doubles(const std::string &name) const {
  std::map<std::string,std::vector<double> >::const_iterator j = m_doubles.find(name);
  if (j != m_doubles.end()) {
    return j->second;
  } else {
    return std::vector<double>();
  }
}

const VariableMetadata::StringAttrs& VariableMetadata::get_all_strings() const {
  return m_strings;
}

const VariableMetadata::DoubleAttrs& VariableMetadata::get_all_doubles() const {
  return m_doubles;
}

//! Set a string attribute.
void VariableMetadata::set_string(const std::string &name, const std::string &value) {

  if (name == "units") {
    // create a dummy object to validate the units string
    Unit tmp(m_unit_system, value);

    m_strings[name] = value;
    m_strings["glaciological_units"] = value;
  } else if (name == "glaciological_units") {
    m_strings[name] = value;

    Unit internal(m_unit_system, get_string("units")),
      glaciological(m_unit_system, value);

    if (not UnitConverter::are_convertible(internal, glaciological)) {
      throw RuntimeError::formatted("units \"%s\" and \"%s\" are not compatible",
                                    get_string("units").c_str(), value.c_str());
    }
  } else if (name == "short_name") {
    set_name(name);
  } else {
    m_strings[name] = value;
  }
}

//! Get a string attribute.
/*!
 * Returns an empty string if an attribute is not set.
 */
std::string VariableMetadata::get_string(const std::string &name) const {
  if (name == "short_name") {
     return get_name();
  }

  std::map<std::string,std::string>::const_iterator j = m_strings.find(name);
  if (j != m_strings.end()) {
    return j->second;
  } else {
    return std::string();
  }
}

void VariableMetadata::report_to_stdout(MPI_Comm com, int verbosity_threshold) const {

  // Print text attributes:
  const VariableMetadata::StringAttrs &strings = this->get_all_strings();
  VariableMetadata::StringAttrs::const_iterator i;
  for (i = strings.begin(); i != strings.end(); ++i) {
    std::string name  = i->first;
    std::string value = i->second;

    if (value.empty()) {
      continue;
    }

    verbPrintf(verbosity_threshold, com, "  %s = \"%s\"\n",
               name.c_str(), value.c_str());
  }

  // Print double attributes:
  const VariableMetadata::DoubleAttrs &doubles = this->get_all_doubles();
  VariableMetadata::DoubleAttrs::const_iterator j;
  for (j = doubles.begin(); j != doubles.end(); ++j) {
    std::string name  = j->first;
    std::vector<double> values = j->second;

    if (values.empty()) {
      continue;
    }

    if ((fabs(values[0]) >= 1.0e7) || (fabs(values[0]) <= 1.0e-4)) {
      // use scientific notation if a number is big or small
      verbPrintf(verbosity_threshold, com, "  %s = %12.3e\n",
                 name.c_str(), values[0]);
    } else {
      verbPrintf(verbosity_threshold, com, "  %s = %12.5f\n",
                 name.c_str(), values[0]);
    }

  }
}

TimeseriesMetadata::TimeseriesMetadata(const std::string &name, const std::string &dimension_name,
                           const UnitSystem &system)
  : VariableMetadata(name, system, 0) {
  m_dimension_name = dimension_name;
}

TimeseriesMetadata::~TimeseriesMetadata()
{
  // empty
}

std::string TimeseriesMetadata::get_dimension_name() const {
  return m_dimension_name;
}

//! Define a NetCDF variable corresponding to a time-series.
void TimeseriesMetadata::define(const PIO &nc, IO_Type nctype, bool) const {

  bool exists = nc.inq_var(get_name());
  if (exists) {
    return;
  }

  nc.redef();

  exists = nc.inq_dim(m_dimension_name);
  if (not exists) {
    nc.def_dim(PISM_UNLIMITED, VariableMetadata(m_dimension_name, m_unit_system));
  }

  exists = nc.inq_var(get_name());
  if (not exists) {
    std::vector<std::string> dims(1, m_dimension_name);
    nc.redef();
    nc.def_var(get_name(), nctype, dims);
  }

  nc.write_attributes(*this, PISM_FLOAT, true);
}

/// TimeBoundsMetadata

TimeBoundsMetadata::TimeBoundsMetadata(const std::string &var_name, const std::string &dim_name,
                           const UnitSystem &system)
  : TimeseriesMetadata(var_name, dim_name, system) {
  m_bounds_name    = "nv";      // number of vertexes
}

TimeBoundsMetadata::~TimeBoundsMetadata() {
  // empty
}

void TimeBoundsMetadata::define(const PIO &nc, IO_Type nctype, bool) const {
  std::vector<std::string> dims;
  bool exists = false;
  
  std::string dimension_name = get_dimension_name();

  exists = nc.inq_var(get_name());
  if (exists) {
    return;
  }

  nc.redef();

  exists = nc.inq_dim(dimension_name);
  if (not exists) {
    nc.def_dim(PISM_UNLIMITED, VariableMetadata(dimension_name, m_unit_system));
  }

  exists = nc.inq_dim(m_bounds_name);
  if (not exists) {
    nc.def_dim(2, VariableMetadata(m_bounds_name, m_unit_system));
  }

  dims.push_back(dimension_name);
  dims.push_back(m_bounds_name);

  nc.redef();

  nc.def_var(get_name(), nctype, dims);

  nc.write_attributes(*this, nctype, true);
}

} // end of namespace pism
