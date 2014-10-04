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

namespace pism {

NCVariable::NCVariable(const std::string &name, const UnitSystem &system, unsigned int ndims)
  : m_n_spatial_dims(ndims), m_units(system), m_glaciological_units(system), m_short_name(name) {

  m_units.reset();
  m_glaciological_units.reset();

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

  int errcode = m_units.parse(new_units);

  assert(errcode == 0 && "invalid units specification");

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

  int errcode = m_glaciological_units.parse(new_units);

  assert(errcode == 0 && "invalid units specification");

  assert(units_are_convertible(m_units, m_glaciological_units) == true);

  return 0;
}

Unit NCVariable::get_units() const {
  return m_units;
}

Unit NCVariable::get_glaciological_units() const {
  return m_glaciological_units;
}

NCSpatialVariable::NCSpatialVariable(const UnitSystem &system)
  : NCVariable("unnamed", system),
    m_x("x", system),
    m_y("y", system),
    m_z("z", system) {

  m_grid = NULL;
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

void NCSpatialVariable::init_2d(const std::string &name, IceGrid &g) {
  std::vector<double> z(1);

  init_3d(name, g, z);
}

//! \brief Initialize a NCSpatialVariable instance.
void NCSpatialVariable::init_3d(const std::string &name, IceGrid &g, std::vector<double> &z_levels) {
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
  PetscErrorCode ierr;

  assert(m_grid != NULL);

  // Find the variable:
  std::string name_found;
  bool found_by_standard_name = false, variable_exists = false;
  ierr = nc.inq_var(get_name(), get_string("standard_name"),
                    variable_exists, name_found, found_by_standard_name); CHKERRQ(ierr);

  if (!variable_exists) {
    ierr = PetscPrintf(m_com,
                       "PISM ERROR: Can't find '%s' (%s) in '%s'.\n",
                       get_name().c_str(),
                       get_string("standard_name").c_str(), nc.inq_filename().c_str());
    CHKERRQ(ierr);
    PISMEnd();
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

    ierr = nc.inq_vardims(name_found, input_dims);
    std::vector<std::string>::iterator j = input_dims.begin();
    while (j != input_dims.end()) {
      AxisType tmp;
      ierr = nc.inq_dimtype(*j, tmp); CHKERRQ(ierr);

      if (tmp != T_AXIS)
        ++input_ndims;

      if (axes.find(tmp) != axes.end())
        ++matching_dim_count;

      ++j;
    }


    if (axes.size() != matching_dim_count) {

      // Join input dimension names:
      j = input_dims.begin();
      std::string tmp = *j++;
      while (j != input_dims.end())
        tmp += std::string(", ") + *j++;

      // Print the error message and stop:
      PetscPrintf(m_com,
                  "PISM ERROR: found the %dD variable %s(%s) in '%s' while trying to read\n"
                  "            '%s' ('%s'), which is %d-dimensional.\n",
                  input_ndims, name_found.c_str(), tmp.c_str(), nc.inq_filename().c_str(),
                  get_name().c_str(), get_string("long_name").c_str(),
                  static_cast<int>(axes.size()));
      PISMEnd();

    }
  }

  unsigned int nlevels = PetscMax(m_zlevels.size(), 1); // make sure we have at least one level
  ierr = nc.get_vec(m_grid, name_found, nlevels, time, v); CHKERRQ(ierr);

  bool input_has_units;
  Unit input_units(get_units().get_system());

  ierr = nc.inq_units(name_found, input_has_units, input_units); CHKERRQ(ierr);

  if (has_attribute("units") && (!input_has_units)) {
    const std::string &units_string = get_string("units"),
      &long_name = get_string("long_name");
    ierr = verbPrintf(2, m_com,
                      "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
                      "              Assuming that it is in '%s'.\n",
                      get_name().c_str(), long_name.c_str(),
                      units_string.c_str()); CHKERRQ(ierr);
    input_units = get_units();
  }

  // Convert data:
  ierr = units_check(get_name(), input_units, get_units()); CHKERRQ(ierr);
  ierr = convert_vec(v, input_units, get_units()); CHKERRQ(ierr);

  return 0;
}

//! \brief Write a \b global Vec `v` to a variable.
/*!
  Defines a variable and converts the units if needed.
 */
PetscErrorCode NCSpatialVariable::write(const PIO &nc, IO_Type nctype,
                                        bool write_in_glaciological_units, Vec v) const {
  PetscErrorCode ierr;

  // find or define the variable
  std::string name_found;
  bool exists, found_by_standard_name;
  ierr = nc.inq_var(get_name(), get_string("standard_name"),
                    exists, name_found, found_by_standard_name); CHKERRQ(ierr);

  if (!exists) {
    ierr = define(nc, nctype, write_in_glaciological_units); CHKERRQ(ierr);
    name_found = get_name();
  }

  if (write_in_glaciological_units) {
    ierr = convert_vec(v, get_units(), get_glaciological_units()); CHKERRQ(ierr);
  }

  // Actually write data:
  unsigned int nlevels = PetscMax(m_zlevels.size(), 1); // make sure we have at least one level
  ierr = nc.put_vec(m_grid, name_found, nlevels, v); CHKERRQ(ierr);

  if (write_in_glaciological_units) {
    ierr = convert_vec(v, get_glaciological_units(), get_units()); CHKERRQ(ierr); // restore the units
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
  PetscErrorCode ierr;
  unsigned int t_length = 0;

  ierr = nc.inq_nrecords(get_name(), get_string("standard_name"),
                         t_length); CHKERRQ(ierr);

  ierr = this->regrid(nc, t_length - 1, flag, do_report_range, default_value, v); CHKERRQ(ierr);

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
  ierr = nc.inq_var(get_name(), get_string("standard_name"),
                    exists, name_found, found_by_standard_name); CHKERRQ(ierr);

  if (exists == true) {                      // the variable was found successfully

    if (flag == OPTIONAL_FILL_MISSING ||
        flag == CRITICAL_FILL_MISSING) {
      ierr = verbPrintf(2, m_com,
                        "PISM WARNING: Replacing missing values with %f [%s] in variable '%s' read from '%s'.\n",
                        default_value, get_string("units").c_str(), get_name().c_str(),
                        nc.inq_filename().c_str()); CHKERRQ(ierr);

      ierr = nc.regrid_vec_fill_missing(m_grid, name_found, m_zlevels,
                                        t_start, default_value,
                                        v); CHKERRQ(ierr);
    } else {
      ierr = nc.regrid_vec(m_grid, name_found, m_zlevels, t_start, v); CHKERRQ(ierr);
    }

    // Now we need to get the units string from the file and convert the units,
    // because check_range and report_range expect the data to be in PISM (MKS)
    // units.

    bool input_has_units;
    Unit input_units(get_units().get_system());

    ierr = nc.inq_units(name_found, input_has_units, input_units); CHKERRQ(ierr);

    if (input_has_units == false) {
      input_units = get_units();
      if (get_string("units") != "") {
        ierr = verbPrintf(2, m_com,
                          "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
                          "              Assuming that it is in '%s'.\n",
                          get_name().c_str(),
                          get_string("long_name").c_str(),
                          get_string("units").c_str()); CHKERRQ(ierr);
      }
    }

    // Convert data:
    ierr = units_check(get_name(), input_units, get_units()); CHKERRQ(ierr);
    ierr = convert_vec(v, input_units, get_units()); CHKERRQ(ierr);

    // Read the valid range info:
    ierr = nc.read_valid_range(name_found, *this); CHKERRQ(ierr);

    // Check the range and warn the user if needed:
    ierr = check_range(nc.inq_filename(), v); CHKERRQ(ierr);

    if (do_report_range == true) {
      // We can report the success, and the range now:
      ierr = verbPrintf(2, m_com, "  FOUND ");
      ierr = this->report_range(v, found_by_standard_name); CHKERRQ(ierr);
    }
  } else {                // couldn't find the variable
    if (flag == CRITICAL ||
        flag == CRITICAL_FILL_MISSING) { // if it's critical, print an error message and stop
      ierr = PetscPrintf(m_com,
                         "PISM ERROR: Can't find '%s' in the regridding file '%s'.\n",
                         get_name().c_str(), nc.inq_filename().c_str());
      CHKERRQ(ierr);
      PISMEnd();
    }

    // If it is optional, fill with the provided default value.
    cv_converter *c = get_glaciological_units().get_converter_from(get_units());
    assert(c != NULL);
    double tmp = cv_convert_double(c, default_value);
    cv_free(c);

    std::string spacer(get_name().size(), ' ');
    ierr = verbPrintf(2, m_com,
                      "  absent %s / %-10s\n"
                      "         %s \\ not found; using default constant %7.2f (%s)\n",
                      get_name().c_str(),
                      get_string("long_name").c_str(),
                      spacer.c_str(), tmp,
                      get_string("glaciological_units").c_str());
    CHKERRQ(ierr);
    ierr = VecSet(v, default_value); CHKERRQ(ierr);
  } // end of if (exists)

  return 0;
}


//! Report the range of a \b global Vec `v`.
PetscErrorCode NCSpatialVariable::report_range(Vec v, bool found_by_standard_name) {
  PetscErrorCode ierr;
  double min, max;

  ierr = VecMin(v, NULL, &min); CHKERRQ(ierr);
  ierr = VecMax(v, NULL, &max); CHKERRQ(ierr);

  cv_converter *c = get_glaciological_units().get_converter_from(get_units());
  assert(c != NULL);
  min = cv_convert_double(c, min);
  max = cv_convert_double(c, max);
  cv_free(c);

  std::string spacer(get_name().size(), ' ');

  if (has_attribute("standard_name")) {

    if (found_by_standard_name) {
      ierr = verbPrintf(2, m_com,
                        " %s / standard_name=%-10s\n"
                        "         %s \\ min,max = %9.3f,%9.3f (%s)\n",
                        get_name().c_str(),
                        get_string("standard_name").c_str(), spacer.c_str(), min, max,
                        get_string("glaciological_units").c_str()); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, m_com,
                        " %s / WARNING! standard_name=%s is missing, found by short_name\n"
                        "         %s \\ min,max = %9.3f,%9.3f (%s)\n",
                        get_name().c_str(),
                        get_string("standard_name").c_str(), spacer.c_str(), min, max,
                        get_string("glaciological_units").c_str()); CHKERRQ(ierr);
    }

  } else {

    ierr = verbPrintf(2, m_com,
                      " %s / %-10s\n"
                      "         %s \\ min,max = %9.3f,%9.3f (%s)\n",
                      get_name().c_str(),
                      get_string("long_name").c_str(), spacer.c_str(), min, max,
                      get_string("glaciological_units").c_str()); CHKERRQ(ierr);
  }

  return 0;
}

//! Check if the range of a \b global Vec `v` is in the range specified by valid_min and valid_max attributes.
PetscErrorCode NCSpatialVariable::check_range(const std::string &filename, Vec v) {
  PetscReal min = 0.0, max = 0.0;
  PetscErrorCode ierr;
  bool failed = false;

  assert(m_grid != NULL);

  // Vec v is always global here (so VecMin and VecMax work as expected)
  ierr = VecMin(v, NULL, &min); CHKERRQ(ierr);
  ierr = VecMax(v, NULL, &max); CHKERRQ(ierr);

  const std::string &units_string = get_string("units");

  if (has_attribute("valid_min") && has_attribute("valid_max")) {
    double valid_min = get_double("valid_min"),
      valid_max = get_double("valid_max");
    if ((min < valid_min) || (max > valid_max)) {
      ierr = PetscPrintf(m_com,
                         "PISM ERROR: some values of '%s' are outside the valid range [%f, %f] (%s).\n",
                         get_name().c_str(), valid_min, valid_max, units_string.c_str()); CHKERRQ(ierr);
      failed = true;
    }

  } else if (has_attribute("valid_min")) {
    double valid_min = get_double("valid_min");
    if (min < valid_min) {
      ierr = PetscPrintf(m_com,
                         "PISM ERROR: some values of '%s' are less than the valid minimum (%f %s).\n",
                         get_name().c_str(), valid_min, units_string.c_str()); CHKERRQ(ierr);
      failed = true;
    }

  } else if (has_attribute("valid_max")) {
    double valid_max = get_double("valid_max");
    if (max > valid_max) {
      ierr = PetscPrintf(m_com,
                         "PISM ERROR: some values of '%s' are greater than the valid maximum (%f %s).\n",
                         get_name().c_str(), valid_max, units_string.c_str()); CHKERRQ(ierr);
      failed = true;
    }
  }

  if (failed == true) {
    PetscPrintf(m_com,
                "            Please inspect variable '%s' in '%s' and replace offending entries.\n"
                "            PISM will stop now.\n",
                get_name().c_str(), filename.c_str());
    PISMEnd();
  }

  return 0;
}

//! \brief Define dimensions a variable depends on.
PetscErrorCode NCSpatialVariable::define_dimensions(const PIO &nc) const {
  PetscErrorCode ierr;
  bool exists;

  // x
  ierr = nc.inq_dim(get_x().get_name(), exists); CHKERRQ(ierr);
  if (!exists) {
    ierr = nc.def_dim(m_grid->Mx, m_x); CHKERRQ(ierr);
    ierr = nc.put_dim(get_x().get_name(), m_grid->x); CHKERRQ(ierr);
  }

  // y
  ierr = nc.inq_dim(get_y().get_name(), exists); CHKERRQ(ierr);
  if (!exists) {
    ierr = nc.def_dim(m_grid->My, m_y); CHKERRQ(ierr);
    ierr = nc.put_dim(get_y().get_name(), m_grid->y); CHKERRQ(ierr);
  }

  // z
  std::string z_name = get_z().get_name();
  if (z_name.empty() == false) {
    ierr = nc.inq_dim(z_name, exists); CHKERRQ(ierr);
    if (!exists) {
      unsigned int nlevels = PetscMax(m_zlevels.size(), 1); // make sure we have at least one level
      ierr = nc.def_dim(nlevels, m_z); CHKERRQ(ierr);
      ierr = nc.put_dim(z_name, m_zlevels); CHKERRQ(ierr);
    }
  }

  return 0;
}

//! Define a NetCDF variable corresponding to a NCVariable object.
PetscErrorCode NCSpatialVariable::define(const PIO &nc, IO_Type nctype,
                                         bool write_in_glaciological_units) const {
  int ierr;
  std::vector<std::string> dims;
  bool exists;

  ierr = nc.inq_var(get_name(), exists); CHKERRQ(ierr);
  if (exists)
    return 0;

  ierr = define_dimensions(nc); CHKERRQ(ierr);
  std::string variable_order = m_variable_order;
  // "..._bounds" should be stored with grid corners (corresponding to
  // the "z" dimension here) last, so we override the variable storage
  // order here
  if (ends_with(get_name(), "_bounds") && variable_order == "zyx")
    variable_order = "yxz";

  std::string x = get_x().get_name(),
    y = get_y().get_name(),
    z = get_z().get_name(),
    t = m_time_dimension_name;

  ierr = nc.redef(); CHKERRQ(ierr);

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

  ierr = nc.def_var(get_name(), nctype, dims); CHKERRQ(ierr);

  ierr = nc.write_attributes(*this, nctype, write_in_glaciological_units); CHKERRQ(ierr);

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
    if (name != "units" && (j->second).empty())
      return false;

    return true;
  }

  if (m_doubles.find(name) != m_doubles.end())
    return true;

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
  if (j != m_doubles.end())
    return (j->second)[0];
  else
    return GSL_NAN;
}

//! Get an array-of-doubles attribute.
std::vector<double> NCVariable::get_doubles(const std::string &name) const {
  std::map<std::string,std::vector<double> >::const_iterator j = m_doubles.find(name);
  if (j != m_doubles.end())
    return j->second;
  else
    return std::vector<double>();
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
  if (name == "short_name") return get_name();

  std::map<std::string,std::string>::const_iterator j = m_strings.find(name);
  if (j != m_strings.end())
    return j->second;
  else
    return std::string();
}

PetscErrorCode NCVariable::report_to_stdout(MPI_Comm com, int verbosity_threshold) const {
  PetscErrorCode ierr;

  // Print text attributes:
  const NCVariable::StringAttrs &strings = this->get_all_strings();
  NCVariable::StringAttrs::const_iterator i;
  for (i = strings.begin(); i != strings.end(); ++i) {
    std::string name  = i->first;
    std::string value = i->second;

    if (value.empty()) continue;

    ierr = verbPrintf(verbosity_threshold, com, "  %s = \"%s\"\n",
                      name.c_str(), value.c_str()); CHKERRQ(ierr);
  }

  // Print double attributes:
  const NCVariable::DoubleAttrs &doubles = this->get_all_doubles();
  NCVariable::DoubleAttrs::const_iterator j;
  for (j = doubles.begin(); j != doubles.end(); ++j) {
    std::string name  = j->first;
    std::vector<double> values = j->second;

    if (values.empty()) continue;

    if ((fabs(values[0]) >= 1.0e7) || (fabs(values[0]) <= 1.0e-4)) {
      // use scientific notation if a number is big or small
      ierr = verbPrintf(verbosity_threshold, com, "  %s = %12.3e\n",
                        name.c_str(), values[0]); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(verbosity_threshold, com, "  %s = %12.5f\n",
                        name.c_str(), values[0]); CHKERRQ(ierr);
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
  PetscErrorCode ierr;

  bool exists;
  ierr = nc.inq_var(get_name(), exists); CHKERRQ(ierr);
  if (exists)
    return 0;

  ierr = nc.redef(); CHKERRQ(ierr);

  ierr = nc.inq_dim(m_dimension_name, exists); CHKERRQ(ierr);
  if (exists == false) {
    NCVariable tmp(m_dimension_name, get_units().get_system());
    ierr = nc.def_dim(PISM_UNLIMITED, tmp); CHKERRQ(ierr);
  }

  ierr = nc.inq_var(get_name(), exists); CHKERRQ(ierr);
  if (exists == false) {
    std::vector<std::string> dims(1);
    dims[0] = m_dimension_name;
    ierr = nc.redef(); CHKERRQ(ierr);
    ierr = nc.def_var(get_name(), nctype, dims); CHKERRQ(ierr);
  }

  ierr = nc.write_attributes(*this, PISM_FLOAT, true); CHKERRQ(ierr);

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
  PetscErrorCode ierr;
  std::vector<std::string> dims;
  bool exists = false;
  
  std::string dimension_name = get_dimension_name();

  UnitSystem system = get_units().get_system();

  ierr = nc.inq_var(get_name(), exists); CHKERRQ(ierr);
  if (exists)
    return 0;

  ierr = nc.redef(); CHKERRQ(ierr);

  ierr = nc.inq_dim(dimension_name, exists); CHKERRQ(ierr);
  if (exists == false) {
    NCVariable tmp(dimension_name, system);
    ierr = nc.def_dim(PISM_UNLIMITED, tmp); CHKERRQ(ierr);
  }

  ierr = nc.inq_dim(m_bounds_name, exists); CHKERRQ(ierr);
  if (exists == false) {
    NCVariable tmp(m_bounds_name, system);
    ierr = nc.def_dim(2, tmp); CHKERRQ(ierr);
  }

  dims.push_back(dimension_name);
  dims.push_back(m_bounds_name);

  ierr = nc.redef(); CHKERRQ(ierr);

  ierr = nc.def_var(get_name(), nctype, dims); CHKERRQ(ierr);

  ierr = nc.write_attributes(*this, nctype, true); CHKERRQ(ierr);

  return 0;
}

} // end of namespace pism
