// Copyright (C) 2009, 2010, 2011, 2012, 2013, 2014 Constantine Khroulev and Ed Bueler
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
#include <assert.h>
#include <gsl/gsl_math.h>

#include "NCVariable.hh"
#include "PIO.hh"
#include "pism_options.hh"
#include "IceGrid.hh"
#include "LocalInterpCtx.hh"
#include "PISMTime.hh"

#include "PISMConfig.hh"
#include "error_handling.hh"

namespace pism {

NCVariable::NCVariable(const std::string &name, const UnitSystem &system, unsigned int ndims)
  : m_n_spatial_dims(ndims),
    m_units(system, "1"),
    m_glaciological_units(system, "1"),
    m_short_name(name) {

  clear_all_strings();
  clear_all_doubles();

  // long_name is unset
  // standard_name is unset
  // pism_intent is unset
  // coordinates is unset

  // valid_min and valid_max are unset
}

NCVariable::~NCVariable() {
  // empty
}

/** Get the number of spatial dimensions.
 */
unsigned int NCVariable::get_n_spatial_dimensions() const {
  return m_n_spatial_dims;
}

//! Set the internal units.
PetscErrorCode NCVariable::set_units(const std::string &new_units) {

  // Do not use NCVariable::set_string here, because it is written in
  // a way that forces users to use set_units() to set units.
  m_strings["units"] = new_units;
  m_strings["glaciological_units"] = new_units;

  m_units = Unit(m_units.get_system(), new_units);

  // Set the glaciological units too (ensures internal consistency):
  m_glaciological_units = m_units;

  return 0;
}

//! Set the glaciological (output) units.
/*! These units are used for output (if write_in_glaciological_units is set)
  and for standard out reports.
 */
PetscErrorCode NCVariable::set_glaciological_units(const std::string &new_units) {
  // Save the human-friendly version of the string; this is to avoid getting
  // things like '3.16887646408185e-08 meter second-1' instead of 'm year-1'
  // (and thus violating the CF conventions).

  // Do not use NCVariable::set_string here, because it is written in
  // a way that forces users to use set_glaciological_units() to set
  // "glaciological" units.
  m_strings["glaciological_units"] = new_units;

  m_glaciological_units = Unit(m_units.get_system(), new_units);

  assert(UnitConverter::are_convertible(m_units, m_glaciological_units) == true);

  return 0;
}

Unit NCVariable::get_units() const {
  return m_units;
}

Unit NCVariable::get_glaciological_units() const {
  return m_glaciological_units;
}

//! 3D version
NCSpatialVariable::NCSpatialVariable(const UnitSystem &system, const std::string &name,
                                     IceGrid &g, std::vector<double> &zlevels)
  : NCVariable("unnamed", system),
    m_x("x", system),
    m_y("y", system),
    m_z("z", system) {

  init_internal(name, g, zlevels);
}

//! 2D version
NCSpatialVariable::NCSpatialVariable(const UnitSystem &system, const std::string &name,
                                     IceGrid &g)
  : NCVariable("unnamed", system),
    m_x("x", system),
    m_y("y", system),
    m_z("z", system) {

  std::vector<double> z(1, 0.0);
  init_internal(name, g, z);
}

void NCSpatialVariable::init_internal(const std::string &name, IceGrid &g,
                                      std::vector<double> &z_levels) {
  m_time_dimension_name = "t";        // will be overriden later

  m_x.set_string("axis", "X");
  m_x.set_string("long_name", "X-coordinate in Cartesian system");
  m_x.set_string("standard_name", "projection_x_coordinate");
  m_x.set_units("m");

  m_y.set_string("axis", "Y");
  m_y.set_string("long_name", "Y-coordinate in Cartesian system");
  m_y.set_string("standard_name", "projection_y_coordinate");
  m_y.set_units("m");

  m_z.set_string("axis", "Z");
  m_z.set_string("long_name", "Z-coordinate in Cartesian system");
  m_z.set_units("m");
  m_z.set_string("positive", "up");

  set_string("coordinates", "lat lon");
  set_string("grid_mapping", "mapping");

  m_variable_order = "xyz";

  set_name(name);
  m_grid = &g;
  m_com = g.com;

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

NCSpatialVariable::NCSpatialVariable(const NCSpatialVariable &other)
  : NCVariable(other), m_x(other.m_x), m_y(other.m_y), m_z(other.m_z) {
  m_time_dimension_name = other.m_time_dimension_name;
  m_variable_order      = other.m_variable_order;
  m_zlevels             = other.m_zlevels;
  m_grid                = other.m_grid;
  m_com                 = other.m_com;
}

NCSpatialVariable::~NCSpatialVariable() {
  // empty
}

void NCSpatialVariable::set_levels(const std::vector<double> &levels) {
  assert(levels.size() >= 1);
  m_zlevels = levels;
}

void NCSpatialVariable::set_time_independent(bool flag) {
  if (flag == true) {
    m_time_dimension_name = "";
  } else {
    m_time_dimension_name = m_grid->config.get_string("time_dimension_name");
  }
}

//! Read a variable from a file into a \b global Vec v.
/*! This also converts the data from input units to internal units if needed.
 */
PetscErrorCode NCSpatialVariable::read(const PIO &nc, unsigned int time, Vec v) {

  assert(m_grid != NULL);

  // Find the variable:
  std::string name_found;
  bool found_by_standard_name = false, variable_exists = false;
  nc.inq_var(get_name(), get_string("standard_name"),
             variable_exists, name_found, found_by_standard_name);

  if (!variable_exists) {
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

      // Join input dimension names:
      j = input_dims.begin();
      std::string tmp = *j++;
      while (j != input_dims.end()) {
        tmp += std::string(", ") + *j++;
      }

      // Print the error message and stop:
      throw RuntimeError::formatted("found the %dD variable %s (%s) in '%s' while trying to read\n"
                                    "'%s' ('%s'), which is %d-dimensional.",
                                    input_ndims, name_found.c_str(), tmp.c_str(), nc.inq_filename().c_str(),
                                    get_name().c_str(), get_string("long_name").c_str(),
                                    static_cast<int>(axes.size()));
    }
  }

  unsigned int nlevels = PetscMax(m_zlevels.size(), 1); // make sure we have at least one level
  nc.get_vec(m_grid, name_found, nlevels, time, v);

  bool input_has_units;
  Unit input_units(get_units().get_system(), "1");

  nc.inq_units(name_found, input_has_units, input_units);

  if (has_attribute("units") && (!input_has_units)) {
    const std::string &units_string = get_string("units"),
      &long_name = get_string("long_name");
    verbPrintf(2, m_com,
               "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
               "              Assuming that it is in '%s'.\n",
               get_name().c_str(), long_name.c_str(),
               units_string.c_str());
    input_units = get_units();
  }

  // Convert data:
  convert_vec(v, input_units, get_units());

  return 0;
}

//! \brief Write a \b global Vec `v` to a variable.
/*!
  Defines a variable and converts the units if needed.
 */
PetscErrorCode NCSpatialVariable::write(const PIO &nc, IO_Type nctype,
                                        bool write_in_glaciological_units, Vec v) const {

  // find or define the variable
  std::string name_found;
  bool exists, found_by_standard_name;
  nc.inq_var(get_name(), get_string("standard_name"),
             exists, name_found, found_by_standard_name);

  if (!exists) {
    define(nc, nctype, write_in_glaciological_units);
    name_found = get_name();
  }

  if (write_in_glaciological_units) {
    convert_vec(v, get_units(), get_glaciological_units());
  }

  // Actually write data:
  unsigned int nlevels = PetscMax(m_zlevels.size(), 1); // make sure we have at least one level
  nc.put_vec(m_grid, name_found, nlevels, v);

  if (write_in_glaciological_units) {
    convert_vec(v, get_glaciological_units(), get_units()); // restore the units
  }

  return 0;
}

//! \brief Regrid from a NetCDF file into a \b global Vec `v`.
/*!
  - if `flag` is `CRITICAL` or `CRITICAL_FILL_MISSING`, stops if the
    variable was not found in the input file
  - if `flag` is one of `CRITICAL_FILL_MISSING` and
    `OPTIONAL_FILL_MISSING`, replace _FillValue with `default_value`.
  - sets `v` to `default_value` if `flag` is `OPTIONAL` and the
    variable was not found in the input file
  - uses the last record in the file
 */
PetscErrorCode NCSpatialVariable::regrid(const PIO &nc, RegriddingFlag flag, bool do_report_range,
                                         double default_value, Vec v) {
  unsigned int t_length = nc.inq_nrecords(get_name(), get_string("standard_name"));

  this->regrid(nc, t_length - 1, flag, do_report_range, default_value, v);

  return 0;
}

PetscErrorCode NCSpatialVariable::regrid(const PIO &nc, unsigned int t_start,
                                         RegriddingFlag flag, bool do_report_range,
                                         double default_value, Vec v) {
  PetscErrorCode ierr;

  assert(m_grid != NULL);

  // Find the variable
  bool exists, found_by_standard_name;
  std::string name_found;
  nc.inq_var(get_name(), get_string("standard_name"),
             exists, name_found, found_by_standard_name);

  if (exists == true) {                      // the variable was found successfully

    if (flag == OPTIONAL_FILL_MISSING ||
        flag == CRITICAL_FILL_MISSING) {
      verbPrintf(2, m_com,
                 "PISM WARNING: Replacing missing values with %f [%s] in variable '%s' read from '%s'.\n",
                 default_value, get_string("units").c_str(), get_name().c_str(),
                 nc.inq_filename().c_str());

      nc.regrid_vec_fill_missing(m_grid, name_found, m_zlevels,
                                 t_start, default_value, v);
    } else {
      nc.regrid_vec(m_grid, name_found, m_zlevels, t_start, v);
    }

    // Now we need to get the units string from the file and convert the units,
    // because check_range and report_range expect the data to be in PISM (MKS)
    // units.

    bool input_has_units;
    Unit input_units(get_units().get_system(), "1");

    nc.inq_units(name_found, input_has_units, input_units);

    if (input_has_units == false) {
      input_units = get_units();
      if (get_string("units") != "") {
        verbPrintf(2, m_com,
                   "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
                   "              Assuming that it is in '%s'.\n",
                   get_name().c_str(),
                   get_string("long_name").c_str(),
                   get_string("units").c_str());
      }
    }

    // Convert data:
    convert_vec(v, input_units, get_units());

    // Read the valid range info:
    nc.read_valid_range(name_found, *this);

    // Check the range and warn the user if needed:
    check_range(nc.inq_filename(), v);

    if (do_report_range == true) {
      // We can report the success, and the range now:
      verbPrintf(2, m_com, "  FOUND ");
      this->report_range(v, found_by_standard_name);
    }
  } else {                // couldn't find the variable
    if (flag == CRITICAL ||
        flag == CRITICAL_FILL_MISSING) { // if it's critical, print an error message and stop
      throw RuntimeError::formatted("Can't find '%s' in the regridding file '%s'.",
                                    get_name().c_str(), nc.inq_filename().c_str());
    }

    // If it is optional, fill with the provided default value.
    assert(UnitConverter::are_convertible(get_units(), get_glaciological_units()) == true);
    UnitConverter c(this->get_units(), this->get_glaciological_units());

    std::string spacer(get_name().size(), ' ');
    verbPrintf(2, m_com,
               "  absent %s / %-10s\n"
               "         %s \\ not found; using default constant %7.2f (%s)\n",
               get_name().c_str(),
               get_string("long_name").c_str(),
               spacer.c_str(), c(default_value),
               get_string("glaciological_units").c_str());
    ierr = VecSet(v, default_value);
    PISM_PETSC_CHK(ierr, "VecSet");
  } // end of if (exists)

  return 0;
}


//! Report the range of a \b global Vec `v`.
PetscErrorCode NCSpatialVariable::report_range(Vec v, bool found_by_standard_name) {
  PetscErrorCode ierr;
  double min, max;

  ierr = VecMin(v, NULL, &min);
  PISM_PETSC_CHK(ierr, "VecMin");
  ierr = VecMax(v, NULL, &max);
  PISM_PETSC_CHK(ierr, "VecMax");

  assert(UnitConverter::are_convertible(get_units(), get_glaciological_units()) == true);
  UnitConverter c(this->get_units(), this->get_glaciological_units());
  min = c(min);
  max = c(max);

  std::string spacer(get_name().size(), ' ');

  if (has_attribute("standard_name")) {

    if (found_by_standard_name) {
      verbPrintf(2, m_com,
                 " %s / standard_name=%-10s\n"
                 "         %s \\ min,max = %9.3f,%9.3f (%s)\n",
                 get_name().c_str(),
                 get_string("standard_name").c_str(), spacer.c_str(), min, max,
                 get_string("glaciological_units").c_str());
    } else {
      verbPrintf(2, m_com,
                 " %s / WARNING! standard_name=%s is missing, found by short_name\n"
                 "         %s \\ min,max = %9.3f,%9.3f (%s)\n",
                 get_name().c_str(),
                 get_string("standard_name").c_str(), spacer.c_str(), min, max,
                 get_string("glaciological_units").c_str());
    }

  } else {

    verbPrintf(2, m_com,
               " %s / %-10s\n"
               "         %s \\ min,max = %9.3f,%9.3f (%s)\n",
               get_name().c_str(),
               get_string("long_name").c_str(), spacer.c_str(), min, max,
               get_string("glaciological_units").c_str());
  }

  return 0;
}

//! Check if the range of a \b global Vec `v` is in the range specified by valid_min and valid_max attributes.
PetscErrorCode NCSpatialVariable::check_range(const std::string &filename, Vec v) {
  PetscReal min = 0.0, max = 0.0;
  PetscErrorCode ierr;

  assert(m_grid != NULL);

  // Vec v is always global here (so VecMin and VecMax work as expected)
  ierr = VecMin(v, NULL, &min);
  PISM_PETSC_CHK(ierr, "VecMin");
  ierr = VecMax(v, NULL, &max);
  PISM_PETSC_CHK(ierr, "VecMax");

  const std::string &units_string = get_string("units");

  if (has_attribute("valid_min") && has_attribute("valid_max")) {
    double
      valid_min = get_double("valid_min"),
      valid_max = get_double("valid_max");
    if ((min < valid_min) || (max > valid_max)) {
      throw RuntimeError::formatted("some values of '%s' in '%s' are outside the valid range [%f, %f] (%s).",
                                    get_name().c_str(), filename.c_str(),
                                    valid_min, valid_max, units_string.c_str());
    }
  } else if (has_attribute("valid_min")) {
    double valid_min = get_double("valid_min");
    if (min < valid_min) {
      throw RuntimeError::formatted("some values of '%s' in '%s' are less than the valid minimum (%f %s).",
                                    get_name().c_str(), filename.c_str(),
                                    valid_min, units_string.c_str());
    }
  } else if (has_attribute("valid_max")) {
    double valid_max = get_double("valid_max");
    if (max > valid_max) {
      throw RuntimeError::formatted("some values of '%s' in '%s' are greater than the valid maximum (%f %s).\n",
                                    get_name().c_str(), filename.c_str(),
                                    valid_max, units_string.c_str());
    }
  }

  return 0;
}

//! \brief Define dimensions a variable depends on.
PetscErrorCode NCSpatialVariable::define_dimensions(const PIO &nc) const {
  bool exists;

  // x
  exists = nc.inq_dim(get_x().get_name());
  if (!exists) {
    nc.def_dim(m_grid->Mx(), m_x);
    nc.put_dim(get_x().get_name(), m_grid->x);
  }

  // y
  exists = nc.inq_dim(get_y().get_name());
  if (!exists) {
    nc.def_dim(m_grid->My, m_y);
    nc.put_dim(get_y().get_name(), m_grid->y);
  }

  // z
  std::string z_name = get_z().get_name();
  if (z_name.empty() == false) {
    exists = nc.inq_dim(z_name);
    if (!exists) {
      unsigned int nlevels = PetscMax(m_zlevels.size(), 1); // make sure we have at least one level
      nc.def_dim(nlevels, m_z);
      nc.put_dim(z_name, m_zlevels);
    }
  }

  return 0;
}

//! Define a NetCDF variable corresponding to a NCVariable object.
PetscErrorCode NCSpatialVariable::define(const PIO &nc, IO_Type nctype,
                                         bool write_in_glaciological_units) const {
  std::vector<std::string> dims;

  bool exists = nc.inq_var(get_name());
  if (exists) {
    return 0;
  }

  define_dimensions(nc);
  std::string variable_order = m_variable_order;
  // "..._bounds" should be stored with grid corners (corresponding to
  // the "z" dimension here) last, so we override the variable storage
  // order here
  if (ends_with(get_name(), "_bounds") && variable_order == "zyx") {
    variable_order = "yxz";
  }

  std::string x = get_x().get_name(),
    y = get_y().get_name(),
    z = get_z().get_name(),
    t = m_time_dimension_name;

  nc.redef();

  if (t.empty() == false) {
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

  if (z.empty() == false) {
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

  return 0;
}

NCVariable& NCSpatialVariable::get_x() {
  return m_x;
}

NCVariable& NCSpatialVariable::get_y() {
  return m_y;
}

NCVariable& NCSpatialVariable::get_z() {
  return m_z;
}

const NCVariable& NCSpatialVariable::get_x() const {
  return m_x;
}

const NCVariable& NCSpatialVariable::get_y() const {
  return m_y;
}

const NCVariable& NCSpatialVariable::get_z() const {
  return m_z;
}

//! Checks if an attribute is present. Ignores empty strings, except
//! for the "units" attribute.
bool NCVariable::has_attribute(const std::string &name) const {

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

void NCVariable::set_name(const std::string &name) {
  m_short_name = name;
}

//! Set a scalar attribute to a single (scalar) value.
void NCVariable::set_double(const std::string &name, double value) {
  m_doubles[name] = std::vector<double>(1, value);
}

//! Set a scalar attribute to a single (scalar) value.
void NCVariable::set_doubles(const std::string &name, const std::vector<double> &values) {
  m_doubles[name] = values;
}

void NCVariable::clear_all_doubles() {
  m_doubles.clear();
}

void NCVariable::clear_all_strings() {
  m_strings.clear();
}

std::string NCVariable::get_name() const {
  return m_short_name;
}

//! Get a single-valued scalar attribute.
double NCVariable::get_double(const std::string &name) const {
  std::map<std::string,std::vector<double> >::const_iterator j = m_doubles.find(name);
  if (j != m_doubles.end()) {
    return (j->second)[0];
  } else {
    return GSL_NAN;
  }
}

//! Get an array-of-doubles attribute.
std::vector<double> NCVariable::get_doubles(const std::string &name) const {
  std::map<std::string,std::vector<double> >::const_iterator j = m_doubles.find(name);
  if (j != m_doubles.end()) {
    return j->second;
  } else {
    return std::vector<double>();
  }
}

const std::map<std::string,std::string>& NCVariable::get_all_strings() const {
  return m_strings;
}

const std::map<std::string,std::vector<double> >& NCVariable::get_all_doubles() const {
  return m_doubles;
}

//! Set a string attribute.
void NCVariable::set_string(const std::string &name, const std::string &value) {

  assert(name != "units");
  assert(name != "glaciological_units");

  if (name == "short_name") {
    set_name(name);
  } else {
    m_strings[name] = value;
  }
}

//! Get a string attribute.
/*!
 * Returns an empty string if an attribute is not set.
 */
std::string NCVariable::get_string(const std::string &name) const {
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

PetscErrorCode NCVariable::report_to_stdout(MPI_Comm com, int verbosity_threshold) const {

  // Print text attributes:
  const NCVariable::StringAttrs &strings = this->get_all_strings();
  NCVariable::StringAttrs::const_iterator i;
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
  const NCVariable::DoubleAttrs &doubles = this->get_all_doubles();
  NCVariable::DoubleAttrs::const_iterator j;
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
  return 0;
}

NCTimeseries::NCTimeseries(const std::string &name, const std::string &dimension_name,
                           const UnitSystem &system)
  : NCVariable(name, system, 0) {
  m_dimension_name = dimension_name;
}

NCTimeseries::~NCTimeseries()
{
  // empty
}

std::string NCTimeseries::get_dimension_name() const {
  return m_dimension_name;
}

//! Define a NetCDF variable corresponding to a time-series.
PetscErrorCode NCTimeseries::define(const PIO &nc, IO_Type nctype, bool) const {

  bool exists = nc.inq_var(get_name());
  if (exists) {
    return 0;
  }

  nc.redef();

  exists = nc.inq_dim(m_dimension_name);
  if (exists == false) {
    NCVariable tmp(m_dimension_name, get_units().get_system());
    nc.def_dim(PISM_UNLIMITED, tmp);
  }

  exists = nc.inq_var(get_name());
  if (exists == false) {
    std::vector<std::string> dims(1);
    dims[0] = m_dimension_name;
    nc.redef();
    nc.def_var(get_name(), nctype, dims);
  }

  nc.write_attributes(*this, PISM_FLOAT, true);

  return 0;
}

/// NCTimeBounds

NCTimeBounds::NCTimeBounds(const std::string &var_name, const std::string &dim_name,
                           const UnitSystem &system)
  : NCTimeseries(var_name, dim_name, system) {
  m_bounds_name    = "nv";      // number of vertexes
}

NCTimeBounds::~NCTimeBounds() {
  // empty
}

PetscErrorCode NCTimeBounds::define(const PIO &nc, IO_Type nctype, bool) const {
  std::vector<std::string> dims;
  bool exists = false;
  
  std::string dimension_name = get_dimension_name();

  UnitSystem system = get_units().get_system();

  exists = nc.inq_var(get_name());
  if (exists) {
    return 0;
  }

  nc.redef();

  exists = nc.inq_dim(dimension_name);
  if (exists == false) {
    NCVariable tmp(dimension_name, system);
    nc.def_dim(PISM_UNLIMITED, tmp);
  }

  exists = nc.inq_dim(m_bounds_name);
  if (exists == false) {
    NCVariable tmp(m_bounds_name, system);
    nc.def_dim(2, tmp);
  }

  dims.push_back(dimension_name);
  dims.push_back(m_bounds_name);

  nc.redef();

  nc.def_var(get_name(), nctype, dims);

  nc.write_attributes(*this, nctype, true);

  return 0;
}

} // end of namespace pism
