// Copyright (C) 2009, 2010, 2011 Constantine Khroulev and Ed Bueler
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

#include "NCVariable.hh"
#include "NCSpatialVariable.hh"
#include "PISMIO.hh"
#include "pism_const.hh"
#include <algorithm>

NCVariable::NCVariable() {
  reset();
}

//! Initialize a NCVariable instance.
void NCVariable::init(string name, MPI_Comm c, PetscMPIInt r) {
  short_name = name;

  com = c;
  rank = r;
}

//! Set the internal units.
/*! Units should not be set by accessing the \c strings member directly. This
  method also checks if \c new_units are valid and initializes the \c units structure.
 */
PetscErrorCode NCVariable::set_units(string new_units) {
  strings["units"] = new_units;

  if (utScan(new_units.c_str(), &units) != 0) {
    SETERRQ2(1, "PISM ERROR: NCVariable '%s': unknown or invalid units specification '%s'.",
	     short_name.c_str(), new_units.c_str());
  }

  // Set the glaciological units too:
  utCopy(&units, &glaciological_units);
  strings["glaciological_units"] = new_units;

  return 0;
}

//! Set the glaciological (output) units.
/*! These units are used for output (if write_in_glaciological_units is set)
  and for standard out reports.

  \c glaciological_units should not be set by accessing the \c strings member
  directly. This method also checks if \c new_units are valid and compatible
  with the internal units.
 */
PetscErrorCode NCVariable::set_glaciological_units(string new_units) {
  double a, b;			// dummy variables
  string &units_string = strings["units"];

  if (utScan(new_units.c_str(), &glaciological_units) != 0) {
    SETERRQ2(1, "PISM ERROR: NCVariable '%s': unknown or invalid units specification '%s'.",
	     short_name.c_str(), new_units.c_str());
  }
  
  if (utConvert(&units, &glaciological_units, &a, &b) == UT_ECONVERT) {
    SETERRQ3(1, "PISM ERROR: NCVariable '%s': attempted to set glaciological units to '%s', which is not compatible with '%s'.\n",
	     short_name.c_str(), new_units.c_str(), units_string.c_str());
  }

  // Save the human-friendly version of the string; this is to avoid getting
  // things like '3.16887646408185e-08 meter second-1' instead of 'm year-1'
  // (and thus violating the CF conventions).
  strings["glaciological_units"] = new_units;
  return 0;
}

PetscErrorCode NCSpatialVariable::reset() {
  NCVariable::reset();

  time_independent = false;
  strings["coordinates"] = "lat lon";
  strings["grid_mapping"] = "mapping";
  return 0;
}

NCSpatialVariable::NCSpatialVariable() {
  grid = NULL;
  reset();
}

void NCSpatialVariable::init(string name, IceGrid &g, GridType d) {
  NCVariable::init(name, g.com, g.rank);
  short_name = name;
  grid = &g;
  dims = d;
}

//! Read a variable from a file into a \b global Vec v.
/*! This also converts the data from input units to internal units if needed.
 */
PetscErrorCode NCSpatialVariable::read(const char filename[], unsigned int time, Vec v) {
  PetscErrorCode ierr;
  bool variable_exists;
  int varid;
  PISMIO nc(grid);

  if (grid == NULL)
    SETERRQ(1, "NCVariable::read: grid is NULL.");

  if (grid->da2 == PETSC_NULL)
    SETERRQ(1, "NCVariable::read: grid.da2 is NULL.");

  // Open the file:
  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
  
  // Find the variable:
  ierr = nc.find_variable(short_name, strings["standard_name"],
			  &varid, variable_exists); CHKERRQ(ierr);

  if (!variable_exists) {
    ierr = PetscPrintf(com,
		      "PISM ERROR: Can't find '%s' (%s) in '%s'.\n",
		       short_name.c_str(),
		       strings["standard_name"].c_str(), filename);
    CHKERRQ(ierr);
    PISMEnd();
  }

  ierr = nc.get_var(varid, v, dims, time); CHKERRQ(ierr);  

  bool input_has_units;
  utUnit input_units;

  ierr = nc.get_units(varid,
		      input_has_units, input_units); CHKERRQ(ierr);

  if ( has("units") && (!input_has_units) ) {
    string &units_string = strings["units"],
      &long_name = strings["long_name"];
    ierr = verbPrintf(2, com,
		      "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
		      "              Assuming that it is in '%s'.\n",
		      short_name.c_str(), long_name.c_str(),
		      units_string.c_str()); CHKERRQ(ierr);
    utCopy(&units, &input_units);
  }

  // Convert data:
  ierr = change_units(v, &input_units, &units); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);
  return 0;
}

//! Write a \b global Vec \c v to a variable.
/*!
  Defines a variable and converts the units if needed.
 */
PetscErrorCode NCSpatialVariable::write(const char filename[], nc_type nctype,
					bool write_in_glaciological_units, Vec v) {
  PetscErrorCode ierr;
  bool exists;
  PISMIO nc(grid);
  int varid;

  ierr = nc.open_for_writing(filename, true, true); CHKERRQ(ierr);
  // append == true and check_dims == true

  // find or define the variable
  ierr = nc.find_variable(short_name, strings["standard_name"],
			  &varid, exists); CHKERRQ(ierr);

  if (!exists) {
    ierr = define(nc, varid, nctype, write_in_glaciological_units); CHKERRQ(ierr);
  }

  if (write_in_glaciological_units) {
    ierr = change_units(v, &units, &glaciological_units); CHKERRQ(ierr);
  }

  // Actually write data:
  ierr = nc.put_var(varid, v, dims); CHKERRQ(ierr);  
  
  if (write_in_glaciological_units) {
    ierr = change_units(v, &glaciological_units, &units); CHKERRQ(ierr); // restore the units
  }

  ierr = nc.close(); CHKERRQ(ierr);
  return 0;
}

//! Regrid from a NetCDF file into a \b global Vec \c v.
/*!
  \li stops if critical == true and the variable was not found
  \li sets \c v to \c default_value if \c set_default_value == true and the variable was not found
  \li interpolation mask can be NULL if it is not used.
 */
PetscErrorCode NCSpatialVariable::regrid(const char filename[], LocalInterpCtx &lic,
					 bool critical, bool set_default_value,
					 PetscScalar default_value,
					 Vec v) {
  int varid;
  bool exists;
  PetscErrorCode ierr;
  PISMIO nc(grid);

  if (grid == NULL)
    SETERRQ(1, "NCVariable::regrid: grid is NULL.");

  if (grid->da2 == PETSC_NULL)
    SETERRQ(1, "NCVariable::regrid: grid.da2 is NULL.");

  // Open the file
  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);

  // Find the variable
  bool found_by_standard_name;
  ierr = nc.find_variable(short_name, strings["standard_name"],
			  &varid, exists,
			  found_by_standard_name); CHKERRQ(ierr);

  if (!exists) {		// couldn't find the variable
    if (critical) {		// if it's critical, print an error message and stop
      // SETERRQ1(1, "Variable '%s' was not found.\n", short_name);
      ierr = PetscPrintf(com,
			"PISM ERROR: Can't find '%s' in the regridding file '%s'.\n",
			 short_name.c_str(), filename);
      CHKERRQ(ierr);
      PISMEnd();
    }

    if (set_default_value) {	// if it's not and we have a default value, set it
      double slope, intercept, tmp;
      utConvert(&units, &glaciological_units, &slope, &intercept);
      tmp = intercept + slope*default_value;
      
      ierr = verbPrintf(2, com, 
			"  absent %-10s/ %-60s\n   %-16s\\ not found; using default constant %7.2f (%s)\n",
			short_name.c_str(),
			strings["long_name"].c_str(),
			"", tmp,
			strings["glaciological_units"].c_str());
      CHKERRQ(ierr);
      ierr = VecSet(v, default_value); CHKERRQ(ierr);
    } else {			// otherwise leave it alone
      ierr = verbPrintf(2, com, 
			"  absent %-10s/ %-60s\n   %-16s\\ not found; continuing without setting it\n",
			short_name.c_str(),
			strings["long_name"].c_str(), "");
      CHKERRQ(ierr);
    }
  } else {			// the variable was found successfully
    ierr = nc.regrid_var(varid, dims, lic, v); CHKERRQ(ierr);

    // Now we need to get the units string from the file and convert the units,
    // because check_range and report_range expect the data to be in PISM (SI)
    // units.

    bool input_has_units;
    utUnit input_units;

    ierr = nc.get_units(varid, input_has_units, input_units); CHKERRQ(ierr);

    if ( has("units") && (!input_has_units) ) {
      ierr = verbPrintf(2, com,
			"PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
			"              Assuming that it is in '%s'.\n",
			short_name.c_str(),
			strings["long_name"].c_str(),
			strings["units"].c_str()); CHKERRQ(ierr);
      utCopy(&units, &input_units);
    }

    // Convert data:
    ierr = change_units(v, &input_units, &units); CHKERRQ(ierr);

    // Read the valid range info:
    ierr = read_valid_range(nc, varid); CHKERRQ(ierr);

    // Check the range and warn the user if needed:
    ierr = check_range(v); CHKERRQ(ierr);

    if (lic.report_range) {
      // We can report the success, and the range now:
      ierr = verbPrintf(2, com, "  FOUND ");
      ierr = report_range(v, found_by_standard_name); CHKERRQ(ierr);
    }
  } // end of if(exists)

  ierr = nc.close(); CHKERRQ(ierr);
  return 0;
}

//! Read the valid range information from a file.
/*! Reads \c valid_min, \c valid_max and \c valid_range attributes; if \c
    valid_range is found, sets the pair \c valid_min and \c valid_max instead.
 */
PetscErrorCode NCVariable::read_valid_range(const NCTool &nc, int varid) {
  string input_units_string;
  utUnit input_units;
  vector<double> bounds;
  double slope, intercept;
  int stat;

  // Never reset valid_min/max if any of them was set internally.
  if (has("valid_min") || has("valid_max"))
    return 0;

  // Read the units: The following code ignores the units in the input file if
  // a) they are absent :-) b) they are invalid c) they are not compatible with
  // internal units.
  stat = nc.get_att_text(varid, "units", input_units_string); CHKERRQ(stat);
  if (input_units_string != "") {
    stat = utScan(input_units_string.c_str(), &input_units);
    if (stat != 0)
      utCopy(&units, &input_units);
  }

  if (utConvert(&input_units, &units, &slope, &intercept) != 0) {
    slope = 1;
    intercept = 0;
  }

  stat = nc.get_att_double(varid, "valid_range", bounds); CHKERRQ(stat);
  if (bounds.size() == 2) {		// valid_range is present
    set("valid_min", intercept + slope*bounds[0]);
    set("valid_max", intercept + slope*bounds[1]);
  } else {			// valid_range has the wrong length or is missing
    stat = nc.get_att_double(varid, "valid_min", bounds); CHKERRQ(stat);
    if (bounds.size() == 1) {		// valid_min is present
      set("valid_min", intercept + slope*bounds[0]);
    }

    stat = nc.get_att_double(varid, "valid_max", bounds); CHKERRQ(stat);
    if (bounds.size() == 1) {		// valid_max is present
      set("valid_max", intercept + slope*bounds[0]);
    }
  }

  return 0;
}

//! Converts \c v from internal to glaciological units.
PetscErrorCode NCSpatialVariable::to_glaciological_units(Vec v) {
  return change_units(v, &units, &glaciological_units);
}

//! Converts \c v from the units corresponding to \c from to the ones corresponding to \c to.
/*!
  Does nothing if this transformation is trivial.
 */
PetscErrorCode NCSpatialVariable::change_units(Vec v, utUnit *from, utUnit *to) {
  PetscErrorCode ierr;
  double slope, intercept;
  string from_name, to_name;
  char *tmp;
  bool use_slope, use_intercept;

  // Get string representations of units:
  utPrint(from, &tmp);
  from_name = tmp;
  utPrint(to, &tmp);
  to_name = tmp;

  // Get the slope and the intercept of the linear transformation.
  ierr = utConvert(from, to, &slope, &intercept);

  if (ierr != 0) { 		// can't convert
    if (ierr == UT_ECONVERT) {	// because units are incompatible
      ierr = PetscPrintf(com,
			 "PISM ERROR: IceModelVec '%s': attempted to convert data from '%s' to '%s'.\n",
			 short_name.c_str(), from_name.c_str(), to_name.c_str());
      PISMEnd();
    } else {			// some other error
      return 2;
    }
  }

  use_slope     = PetscAbsReal(slope - 1.0) > 1e-16;
  use_intercept = PetscAbsReal(intercept)   > 1e-16;

  if (use_slope && use_intercept) {
    ierr = VecScale(v, slope); CHKERRQ(ierr);
    ierr = VecShift(v, intercept); CHKERRQ(ierr);
  } else if (use_slope && !use_intercept) {
    ierr = VecScale(v, slope); CHKERRQ(ierr);
  } else if (!use_slope && use_intercept) {
    ierr = VecShift(v, intercept); CHKERRQ(ierr);
  }

  return 0;
}

//! Write variable attributes to a NetCDF file.
/*!

  \li if write_in_glaciological_units == true, "glaciological_units" are
  written under the name "units" plus the valid range is written in
  glaciological units.

  \li if both valid_min and valid_max are set, then valid_range is written
  instead of the valid_min, valid_max pair.
 */
PetscErrorCode NCVariable::write_attributes(const NCTool &nc, int varid, nc_type nctype,
					    bool write_in_glaciological_units) const {
  int ierr;

  if (rank != 0) return 0;

  int ncid = nc.get_ncid();

  ierr = nc.define_mode(); CHKERRQ(ierr);

  // units, valid_min, valid_max and valid_range need special treatment:
  if (has("units")) {
    string output_units = get_string("units");
    if (write_in_glaciological_units)
      output_units = get_string("glaciological_units");
    ierr = nc_put_att_text(ncid, varid,
			   "units",
			   output_units.size(), output_units.c_str());
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  double bounds[2], fill_value = 0.0;

  if (has("_FillValue")) {
    fill_value = get("_FillValue");
  }

  // We need to save valid_min, valid_max and valid_range in the units
  // matching the ones in the output.
  if (write_in_glaciological_units) {
    double slope, intercept;

    ierr = utConvert(&units, &glaciological_units, &slope, &intercept); CHKERRQ(ierr);

    bounds[0]  = intercept + slope*get("valid_min");
    bounds[1]  = intercept + slope*get("valid_max");
    fill_value = intercept + slope*fill_value;
  } else {
    bounds[0] = get("valid_min");
    bounds[1] = get("valid_max");
  }

  if (has("_FillValue")) {
    ierr = nc_put_att_double(ncid, varid, "_FillValue", nctype, 1, &fill_value);
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  if (has("valid_min") && has("valid_max")) {
    ierr = nc_put_att_double(ncid, varid, "valid_range", nctype, 2, bounds);
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  } else if (has("valid_min")) {
    ierr = nc_put_att_double(ncid, varid, "valid_min", nctype, 1, &bounds[0]);
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  } else if (has("valid_max")) {
    ierr = nc_put_att_double(ncid, varid, "valid_max", nctype, 1, &bounds[1]);
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  // Write text attributes:
  map<string, string>::const_iterator i;
  for (i = strings.begin(); i != strings.end(); ++i) {
    string name  = i->first;
    string value = i->second;

    if (name == "units" ||
	name == "glaciological_units" ||
	value.empty())
      continue;

    ierr = nc_put_att_text(ncid, varid, name.c_str(), value.size(), value.c_str()); 
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  // Write double attributes:
  map<string, vector<double> >::const_iterator j;
  for (j = doubles.begin(); j != doubles.end(); ++j) {
    string name  = j->first;
    vector<double> values = j->second;

    if (name == "valid_min" ||
	name == "valid_max" ||
	name == "valid_range" ||
	name == "_FillValue" ||
	values.empty())
      continue;

    ierr = nc_put_att_double(ncid, varid, name.c_str(), nctype, values.size(), &values[0]);
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  return 0;
}


//! Report the range of a \b global Vec \c v.
PetscErrorCode NCSpatialVariable::report_range(Vec v, bool found_by_standard_name) {
  double slope, intercept;
  PetscErrorCode ierr;
  PetscReal min, max;

  // Get the conversion coefficients:
  utConvert(&units, &glaciological_units, &slope, &intercept);

  ierr = VecMin(v, PETSC_NULL, &min); CHKERRQ(ierr);
  ierr = VecMax(v, PETSC_NULL, &max); CHKERRQ(ierr);

  // Note that in some cases the following conversion does nothing
  min = min * slope + intercept;
  max = max * slope + intercept;

  if (has("standard_name")) {

    if (found_by_standard_name) {
      ierr = verbPrintf(2, com, 
			" %-10s/ standard_name=%-60s\n   %-16s\\ min,max = %9.3f,%9.3f (%s)\n",
			short_name.c_str(),
			strings["standard_name"].c_str(), "", min, max,
			strings["glaciological_units"].c_str()); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, com, 
			" %-10s/ WARNING! standard_name=%s is missing, found by short_name\n   %-16s\\ min,max = %9.3f,%9.3f (%s)\n",
			short_name.c_str(),
			strings["standard_name"].c_str(), "", min, max,
			strings["glaciological_units"].c_str()); CHKERRQ(ierr);
    }

  } else {

    ierr = verbPrintf(2, com, 
		      " %-10s/ %-60s\n   %-16s\\ min,max = %9.3f,%9.3f (%s)\n",
		      short_name.c_str(),
		      strings["long_name"].c_str(), "", min, max,
		      strings["glaciological_units"].c_str()); CHKERRQ(ierr);
  }

  return 0;
}

//! Check if the range of a \b global Vec \c v is in the range specified by valid_min and valid_max attributes.
PetscErrorCode NCSpatialVariable::check_range(Vec v) {
  PetscScalar min, max;
  PetscErrorCode ierr;

  if (grid == NULL)
    SETERRQ(1, "NCVariable::check_range: grid is NULL.");

  ierr = VecMin(v, PETSC_NULL, &min); CHKERRQ(ierr);
  ierr = VecMax(v, PETSC_NULL, &max); CHKERRQ(ierr);

  string &units_string = strings["units"];

  if (has("valid_min") && has("valid_max")) {
    double valid_min = get("valid_min"),
      valid_max = get("valid_max");
    if ((min < valid_min) || (max > valid_max))
      ierr = verbPrintf(2, com,
			"PISM WARNING: some values of '%s' are outside the valid range [%f, %f] (%s)\n",
			short_name.c_str(), valid_min, valid_max, units_string.c_str()); CHKERRQ(ierr);
    
  } else if (has("valid_min")) {
    double valid_min = get("valid_min");
    if (min < valid_min) {
      ierr = verbPrintf(2, com,
			"PISM WARNING: some values of '%s' are less than the valid minimum %f (%s)\n",
			short_name.c_str(), valid_min, units_string.c_str()); CHKERRQ(ierr);
    }
    
  } else if (has("valid_max")) {
    double valid_max = get("valid_max");
    if (max > valid_max) {
      ierr = verbPrintf(2, com,
			"PISM WARNING: some values of '%s' are greater than the valid maximum %f (%s)\n",
			short_name.c_str(), valid_max, units_string.c_str()); CHKERRQ(ierr);
    }
  }
  return 0;
}

//! Define a NetCDF variable corresponding to a NCVariable object.
PetscErrorCode NCSpatialVariable::define(const NCTool &nc, int &varid, nc_type nctype,
                                         bool write_in_glaciological_units) const {
  int stat, i = 0, ndims, dimids[4],
    ncid = nc.get_ncid();

  bool exists;
  stat = nc.find_variable(short_name, &varid, exists); CHKERRQ(stat); 
  if (exists) return 0;

  if (rank != 0) {
    varid = 0;
    return 0;
  }

  stat = nc.define_mode(); CHKERRQ(stat); 

  if (!time_independent) {
    stat = nc_inq_dimid(ncid, "t", &dimids[i++]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  // Use t,x,y,z(zb) variable order: it is weird, but matches in-memory storage
  // order and is *a lot* faster.
#if (PISM_FAST_WRITE == 1)
  stat = nc_inq_dimid(ncid, "x", &dimids[i++]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  stat = nc_inq_dimid(ncid, "y", &dimids[i++]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
#endif

  // Use the t,y,x,z variable order: also weird, somewhat slower, but 2D fields
  // are stored in the "natural" order.
#if  (PISM_FAST_WRITE == 2)
  stat = nc_inq_dimid(ncid, "y", &dimids[i++]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  stat = nc_inq_dimid(ncid, "x", &dimids[i++]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
#endif

  switch (dims) {
  case GRID_3D:
    stat = nc_inq_dimid(ncid, "z", &dimids[i++]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    ndims = 4 - time_independent;
    break;
  case GRID_3D_BEDROCK:
    stat = nc_inq_dimid(ncid, "zb", &dimids[i++]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    ndims = 4 - time_independent;
    break;
  case GRID_2D:
  default:
    ndims = 3 - time_independent;
    break;
  }

  // Use the t,z(zb),y,x variables order: more natural for plotting and post-processing,
  // but requires transposing data while writing and is *a lot* slower.
#if (PISM_FAST_WRITE == 3)
  stat = nc_inq_dimid(ncid, "y", &dimids[i++]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  stat = nc_inq_dimid(ncid, "x", &dimids[i++]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
#endif
    
  stat = nc_def_var(ncid, short_name.c_str(), nctype, ndims, dimids, &varid);
  CHKERRQ(check_err(stat,__LINE__,__FILE__));

  // this is all we need to turn on compression:
  // stat = nc_def_var_deflate(ncid, varid, 0, 1, 5);
  // CHKERRQ(check_err(stat,__LINE__,__FILE__));

  stat = write_attributes(nc, varid, nctype, write_in_glaciological_units); CHKERRQ(stat);

  return 0;
}

//! Reset all the attributes.
PetscErrorCode NCVariable::reset() {

  strings.clear();
  short_name = "unnamed_variable";
  // long_name is unset
  // standard_name is unset
  // pism_intent is unset
  // coordinates is unset

  doubles.clear();
  // valid_min and valid_max are unset

  utClear(&units);
  utClear(&glaciological_units);

  return 0;
}

//! Checks if an attribute is present. Ignores empty strings.
bool NCVariable::has(string name) const {

  map<string,string>::const_iterator j = strings.find(name);
  if (j != strings.end()) {
    if ((j->second).empty())
      return false;
		   
    return true;
  }

  if (doubles.find(name) != doubles.end())
    return true;

  return false;
}

//! Set a scalar attribute to a single value.
void NCVariable::set(string name, double value) {
  doubles[name] = vector<double>(1, value);
}

//! Get a single-valued scalar attribute.
/*! Returns 0 if an attribute is not present, so that has() works correctly.
 */
double NCVariable::get(string name) const {
  map<string,vector<double> >::const_iterator j = doubles.find(name);
  if (j != doubles.end())
    return (j->second)[0];
  else
    return 0;
}

//! Set a string attribute.
void NCVariable::set_string(string name, string value) {
  if (name == "short_name") {
    short_name = name;
    return;
  }

  strings[name] = value;
}

//! Get a string attribute.
string NCVariable::get_string(string name) const {
  if (name == "short_name") return short_name;

  map<string,string>::const_iterator j = strings.find(name);
  if (j != strings.end())
    return j->second;
  else
    return string();
}

//! Check if a value \c a is in the valid range defined by \c valid_min and \c valid_min attributes.
bool NCVariable::is_valid(PetscScalar a) const {
  
  if (has("valid_min") && has("valid_max"))
    return (a >= get("valid_min")) && (a <= get("valid_max"));

  if (has("valid_min"))
    return a >= get("valid_min");

  if (has("valid_max"))       
    return a <= get("valid_max");

  return true;
}

//! Read boolean flags and double parameters from a NetCDF file.
/*!
  Erases all the present parameters before reading.
 */
PetscErrorCode NCConfigVariable::read(const char filename[]) {

  PetscErrorCode ierr;
  bool variable_exists;
  int varid, nattrs;
  NCTool nc(com, rank);

  strings.clear();
  doubles.clear();
  config_filename = filename;

  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);

  ierr = nc.find_variable(short_name, &varid, variable_exists); CHKERRQ(ierr);

  if (!variable_exists) {
    ierr = PetscPrintf(com,
		       "PISM ERROR: configuration variable %s was not found in %s.\n"
		       "            Exiting...\n",
		       short_name.c_str(), filename); CHKERRQ(ierr);
    PISMEnd();
  }

  ierr = nc.inq_nattrs(varid, nattrs); CHKERRQ(ierr);

  for (int j = 0; j < nattrs; ++j) {
    string attname;
    nc_type nctype;
    ierr = nc.inq_att_name(varid, j, attname); CHKERRQ(ierr);
    ierr = nc.inq_att_type(varid, attname.c_str(), nctype); CHKERRQ(ierr);

    if (nctype == NC_CHAR) {
      string value;
      ierr = nc.get_att_text(varid, attname.c_str(), value); CHKERRQ(ierr);

      strings[attname] = value;
    } else {
      vector<double> values;

      ierr = nc.get_att_double(varid, attname.c_str(), values); CHKERRQ(ierr);
      doubles[attname] = values;
    }
  } // end of for (int j = 0; j < nattrs; ++j)

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

//! Write a config variable to a file (with all its attributes).
PetscErrorCode NCConfigVariable::write(const char filename[]) const {
  PetscErrorCode ierr;
  int varid;
  bool variable_exists;
  NCTool nc(com, rank);

  ierr = nc.open_for_writing(filename); CHKERRQ(ierr);

  ierr = nc.find_variable(short_name, &varid, variable_exists); CHKERRQ(ierr);

  if (!variable_exists) {
    ierr = define(nc, varid, NC_BYTE, false); CHKERRQ(ierr);
  }

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

//! Define a configuration NetCDF variable.
PetscErrorCode NCConfigVariable::define(const NCTool &nc, int &varid, nc_type type,
                                        bool) const {
  int stat, var_id,
    ncid = nc.get_ncid();

  bool exists;
  stat = nc.find_variable(short_name, &varid, exists); CHKERRQ(stat); 
  if (exists) return 0;

  if (rank == 0) {
    stat = nc.define_mode(); CHKERRQ(stat); 

    stat = nc_def_var(ncid, short_name.c_str(), type, 0, NULL, &var_id);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  stat = MPI_Bcast(&var_id, 1, MPI_INT, 0, com); CHKERRQ(stat);

  varid = var_id;

  stat = write_attributes(nc, varid, NC_DOUBLE, false); CHKERRQ(stat); 

  return 0;
}


//! Returns a \c double parameter. Stops if it was not found.
double NCConfigVariable::get(string name) const {

  if (doubles.find(name) != doubles.end()) {
    return NCVariable::get(name);
  } else {
    PetscPrintf(com, "PISM ERROR: parameter '%s' is unset. (Parameters read from '%s'.)\n",
		name.c_str(), config_filename.c_str());
    PISMEnd();
  }

  return 0;			// can't happen
}

//! Returns a boolean flag by name. Unset flags are treated as if they are set to 'false'.
/*!
  Strings "false", "no", "off" are interpreted as 'false'; "true", "on", "yes" -- as 'true'.

  Any other string produces an error.
 */
bool NCConfigVariable::get_flag(string name) const {

  map<string,string>::const_iterator j = strings.find(name);
  if (j != strings.end()) {
    const string &value = j->second;

    if ((value == "false") ||
	(value == "no") ||
	(value == "off"))
      return false;

    if ((value == "true") ||
	(value == "yes") ||
	(value == "on"))
      return true;

    PetscPrintf(com,
		"PISM ERROR: Parameter '%s' (%s) cannot be interpreted as a boolean.\n"
		"            Please make sure that it is equal to one of 'true', 'yes', 'on', 'false', 'no', 'off'.\n",
		name.c_str(), value.c_str());
    PISMEnd();
  }

  PetscPrintf(com,
	      "PISM ERROR: Parameter '%s' was not set. (Read from '%s'.)\n",
	      name.c_str(), config_filename.c_str());
  PISMEnd();

  return true;			// will never happen
}

string NCConfigVariable::get_string(string name) const {

  map<string,string>::const_iterator j = strings.find(name);
  if (j != strings.end())
    return j->second;
  else {
    PetscPrintf(com,
		"PISM ERROR: Parameter '%s' was not set. (Read from '%s'.)\n",
		name.c_str(), config_filename.c_str());
    PISMEnd();
  }

  return string();		// will never happen
}

//! Set a value of a boolean flag.
void NCConfigVariable::set_flag(string name, bool value) {
  if (value) strings[name] = "true";
  else       strings[name] = "false";
}

//! Write attributes to a NetCDF variable. All attributes are equal here.
PetscErrorCode NCConfigVariable::write_attributes(const NCTool &nc, int varid, nc_type nctype,
						  bool /*write_in_glaciological_units*/) const {
  int ierr, ncid;

  if (rank != 0) return 0;

  ncid = nc.get_ncid();

  ierr = nc.define_mode(); CHKERRQ(ierr); 

  // Write text attributes:
  map<string, string>::const_iterator i;
  for (i = strings.begin(); i != strings.end(); ++i) {
    string name  = i->first;
    string value = i->second;

    if (value.empty()) continue;

    ierr = nc_put_att_text(ncid, varid, name.c_str(), value.size(), value.c_str()); 
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  // Write double attributes:
  map<string, vector<double> >::const_iterator j;
  for (j = doubles.begin(); j != doubles.end(); ++j) {
    string name  = j->first;
    vector<double> values = j->second;

    if (values.empty()) continue;

    ierr = nc_put_att_double(ncid, varid, name.c_str(), nctype, values.size(), &values[0]);
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  return 0;
}

//! Get a flag from a command-line option.
/*!
  If called as flag_from_option("foo", "foo"), checks both -foo and -no_foo.

  \li if -foo is set, calls set_flag("foo", true),

  \li if -no_foo is set, calls set_flag("foo", false),
  
  \li if both are set, prints an error message and stops,

  \li if none, does nothing.
  
 */
PetscErrorCode NCConfigVariable::flag_from_option(string name, string flag) {
  PetscErrorCode ierr;
  bool foo = false,
    no_foo = false;

  ierr = PISMOptionsIsSet("-" + name, get_string(flag + "_doc"), foo); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-no_" + name, no_foo); CHKERRQ(ierr);

  if (foo && no_foo) {
    PetscPrintf(com, "PISM ERROR: Inconsistent command-line options: both -%s and -no_%s are set.\n",
		name.c_str(), name.c_str());
    PISMEnd();
  }

  if (foo)
    set_flag(flag, true);

  if (no_foo)
    set_flag(flag, false);

  return 0;
}

//! Sets a configuration parameter from a command-line option.
/*!
  If called as scalar_from_option("foo", "foo"), checks -foo and calls set("foo", value).

  Does nothing if -foo was not set.

  Note that no unit conversion is performed; parameters should be stored in
  input units and converted as needed. (This allows saving parameters without
  converting again.)
 */
PetscErrorCode NCConfigVariable::scalar_from_option(string name, string parameter) {
  PetscErrorCode ierr;
  PetscReal value = get(parameter);
  bool flag;
  
  ierr = PISMOptionsReal("-" + name,
			 get_string(parameter + "_doc"),
			 value, flag); CHKERRQ(ierr);
  if (flag)
    this->set(parameter, value);
  
  return 0;
}

//! Print all the attributes of a configuration variable.
PetscErrorCode NCConfigVariable::print(PetscInt vt) const {
  PetscErrorCode ierr;

  ierr = verbPrintf(vt, com, "PISM parameters read from %s:\n",
		    config_filename.c_str());

  // Print text attributes:
  map<string, string>::const_iterator i;
  for (i = strings.begin(); i != strings.end(); ++i) {
    string name  = i->first;
    string value = i->second;

    if (value.empty()) continue;

    ierr = verbPrintf(vt, com, "  %s = \"%s\"\n",
		      name.c_str(), value.c_str()); CHKERRQ(ierr);
  }

  // Print double attributes:
  map<string, vector<double> >::const_iterator j;
  for (j = doubles.begin(); j != doubles.end(); ++j) {
    string name  = j->first;
    vector<double> values = j->second;

    if (values.empty()) continue;
    
    if ((fabs(values[0]) >= 1.0e7) || (fabs(values[0]) <= 1.0e-4)) {
      // use scientific notation if a number is big or small
      ierr = verbPrintf(vt, com, "  %s = %12.3e\n",
		        name.c_str(), values[0]); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(vt, com, "  %s = %12.5f\n",
		        name.c_str(), values[0]); CHKERRQ(ierr);
    }

  }

  return 0;
}


string NCConfigVariable::get_config_filename() const {
  return config_filename;
}


//! Imports values from the other config variable, silently overwriting present values.
void NCConfigVariable::import_from(const NCConfigVariable &other) {
  map<string, vector<double> >::const_iterator j;
  for (j = other.doubles.begin(); j != other.doubles.end(); ++j)
    doubles[j->first] = j->second;

  map<string,string>::const_iterator k;
  for (k = other.strings.begin(); k != other.strings.end(); ++k)
    strings[k->first] = k->second;
}

//! Update values from the other config variable, overwriting present values but avoiding adding new ones.
void NCConfigVariable::update_from(const NCConfigVariable &other) {
  map<string, vector<double> >::iterator j;
  map<string, vector<double> >::const_iterator i;
  for (j = doubles.begin(); j != doubles.end(); ++j) {
    i = other.doubles.find(j->first);
    if (i != other.doubles.end())
      j->second = i->second;
  }

  map<string,string>::iterator k;
  map<string,string>::const_iterator m;
  for (k = strings.begin(); k != strings.end(); ++k) {
    m = other.strings.find(k->first);
    if (m != other.strings.end())
      k->second = m->second;
  }
}

void NCTimeseries::init(string n, string dim_name, MPI_Comm c, PetscMPIInt r) {
  NCVariable::init(n, c, r);
  dimension_name = dim_name;
}


//! Read a time-series variable from a NetCDF file to a vector of doubles.
PetscErrorCode NCTimeseries::read(const char filename[], vector<double> &data) {

  PetscErrorCode ierr;
  NCTool nc(com, rank);
  int ncid, varid;
  bool variable_exists;
  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);

  ncid = nc.get_ncid();

  // Find the variable:
  ierr = nc.find_variable(short_name, strings["standard_name"],
			  &varid, variable_exists); CHKERRQ(ierr);

  if (!variable_exists) {
    ierr = PetscPrintf(com,
		      "PISM ERROR: Can't find '%s' (%s) in '%s'.\n",
		       short_name.c_str(),
		       strings["standard_name"].c_str(), filename);
    CHKERRQ(ierr);
    PISMEnd();
  }

  vector<int> dimids;
  ierr = nc.inq_dimids(varid, dimids); CHKERRQ(ierr);

  if (dimids.size() != 1) {
    ierr = PetscPrintf(com,
		       "PISM ERROR: Variable '%s' in '%s' depends on %d dimensions,\n"
		       "            but a time-series variable can only depend on 1 dimension.\n",
		       short_name.c_str(), filename, dimids.size()); CHKERRQ(ierr);
    PISMEnd();
  }

  ierr = nc.inq_dimname(dimids[0], dimension_name); CHKERRQ(ierr);

  int length;
  ierr = nc.get_dim_length(dimension_name.c_str(), &length); CHKERRQ(ierr);

  if (length <= 0) {
    ierr = PetscPrintf(com,
		       "PISM ERROR: Dimension %s has zero (or negative) length!\n",
		       dimension_name.c_str()); CHKERRQ(ierr);
    PISMEnd();
  }

  data.resize(length);		// memory allocation happens here

  ierr = nc.data_mode(); CHKERRQ(ierr);

  if (rank == 0) {
    ierr = nc_get_var_double(ncid, varid, &data[0]); 
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }
  ierr = MPI_Bcast(&data[0], length, MPI_DOUBLE, 0, com); CHKERRQ(ierr);

  bool input_has_units;
  utUnit input_units;

  ierr = nc.get_units(varid,
		      input_has_units, input_units); CHKERRQ(ierr);

  if ( has("units") && (!input_has_units) ) {
    string &units_string = strings["units"],
      &long_name = strings["long_name"];
    ierr = verbPrintf(2, com,
		      "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
		      "              Assuming that it is in '%s'.\n",
		      short_name.c_str(), long_name.c_str(),
		      units_string.c_str()); CHKERRQ(ierr);
    utCopy(&units, &input_units);
  }

  ierr = change_units(data, &input_units, &units); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode NCTimeseries::report_range(vector<double> &data) {
  double slope, intercept;
  PetscErrorCode ierr;
  PetscReal min, max;

  // Get the conversion coefficients:
  utConvert(&units, &glaciological_units, &slope, &intercept);

  min = *min_element(data.begin(), data.end());
  max = *max_element(data.begin(), data.end());

  min = min * slope + intercept;
  max = max * slope + intercept;

  ierr = verbPrintf(2, com, 
		    "  FOUND  %-10s/ %-60s\n   %-16s\\ min,max = %9.3f,%9.3f (%s)\n",
		    short_name.c_str(),
		    strings["long_name"].c_str(), "", min, max,
		    strings["glaciological_units"].c_str()); CHKERRQ(ierr);

  return 0;
}

//! Define a NetCDF variable corresponding to a time-series.
PetscErrorCode NCTimeseries::define(const NCTool &nc, int &varid, nc_type nctype,
                                    bool) const {
  PetscErrorCode ierr;
  int dimid, ncid = nc.get_ncid();

  bool exists;
  ierr = nc.find_variable(short_name, &varid, exists); CHKERRQ(ierr); 
  if (exists) return 0;

  if (rank == 0) {
    ierr = nc.define_mode(); CHKERRQ(ierr);

    ierr = nc_inq_dimid(ncid, dimension_name.c_str(), &dimid);
    if (ierr != NC_NOERR) {
      ierr = nc_def_dim(ncid, dimension_name.c_str(), NC_UNLIMITED, &dimid);
      CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    }

    ierr = nc_def_var(ncid, short_name.c_str(), nctype, 1, &dimid, &varid);
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }
  ierr = MPI_Bcast(&varid, 1, MPI_INT, 0, com); CHKERRQ(ierr);

  ierr = write_attributes(nc, varid, NC_FLOAT, true);

  return 0;
}

PetscErrorCode NCTimeseries::write(const char filename[], size_t start,
				   vector<double> &data, nc_type nctype) {

  PetscErrorCode ierr;
  NCTool nc(com, rank);
  bool variable_exists = false;
  int varid = -1;
 
  ierr = nc.open_for_writing(filename); CHKERRQ(ierr);

  ierr = nc.find_variable(short_name.c_str(), &varid, variable_exists); CHKERRQ(ierr);
  if (!variable_exists) {
    ierr = define(nc, varid, nctype, true); CHKERRQ(ierr);
  }

  // convert to glaciological units:
  ierr = change_units(data, &units, &glaciological_units); CHKERRQ(ierr);

  ierr = nc.data_mode(); CHKERRQ(ierr);

  size_t count = static_cast<size_t>(data.size());
  if (rank == 0) {
    ierr = nc_put_vara_double(nc.get_ncid(), varid, &start, &count, &data[0]);
  }


  ierr = nc.close(); CHKERRQ(ierr);

  // restore internal units:
  ierr = change_units(data, &glaciological_units, &units); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode NCTimeseries::write(const char filename[], size_t start,
				   double data, nc_type nctype) {
  vector<double> tmp(1);
  tmp[0] = data;
  return write(filename, start, tmp, nctype);
}

//! Convert \c data.
PetscErrorCode NCTimeseries::change_units(vector<double> &data, utUnit *from, utUnit *to) {
  PetscErrorCode ierr;
  double slope, intercept;
  string from_name, to_name;
  char *tmp;
  bool use_slope, use_intercept;

  // Get string representations of units:
  utPrint(from, &tmp);
  from_name = tmp;
  utPrint(to, &tmp);
  to_name = tmp;

  // Get the slope and the intercept of the linear transformation.
  ierr = utConvert(from, to, &slope, &intercept);

  if (ierr != 0) { 		// can't convert
    if (ierr == UT_ECONVERT) {	// because units are incompatible
      ierr = PetscPrintf(com,
			 "PISM ERROR: variable '%s': attempted to convert data from '%s' to '%s'.\n",
			 short_name.c_str(), from_name.c_str(), to_name.c_str());
      PISMEnd();
    } else {			// some other error
      return 2;
    }
  }

  use_slope     = PetscAbsReal(slope - 1.0) > 1e-16;
  use_intercept = PetscAbsReal(intercept)   > 1e-16;

  vector<double>::iterator j;
  if (use_slope) {
    for (j = data.begin(); j < data.end(); ++j)
      *j *= slope;
  }

  if (use_intercept)
    for (j = data.begin(); j < data.end(); ++j)
      *j += intercept;

  return 0;
}

PetscErrorCode NCGlobalAttributes::read(const char filename[]) {
  PetscErrorCode ierr;
  int nattrs;
  NCTool nc(com, rank);

  strings.clear();
  doubles.clear();
  config_filename = filename;

  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);

  ierr = nc.inq_nattrs(NC_GLOBAL, nattrs); CHKERRQ(ierr);

  for (int j = 0; j < nattrs; ++j) {
    string attname;
    nc_type nctype;
    ierr = nc.inq_att_name(NC_GLOBAL, j, attname); CHKERRQ(ierr);
    ierr = nc.inq_att_type(NC_GLOBAL, attname.c_str(), nctype); CHKERRQ(ierr);

    if (nctype == NC_CHAR) {
      string value;
      ierr = nc.get_att_text(NC_GLOBAL, attname.c_str(), value); CHKERRQ(ierr);

      strings[attname] = value;
    } else {
      vector<double> values;

      ierr = nc.get_att_double(NC_GLOBAL, attname.c_str(), values); CHKERRQ(ierr);
      doubles[attname] = values;
    }
  } // end of for (int j = 0; j < nattrs; ++j)

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

//! Writes global attributes to a file by calling write_attributes().
PetscErrorCode NCGlobalAttributes::write(const char filename[]) const {
  PetscErrorCode ierr;
  NCTool nc(com, rank);

  ierr = nc.open_for_writing(filename); CHKERRQ(ierr);

  ierr = write_attributes(nc, NC_GLOBAL, NC_DOUBLE, false); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

//! Writes global attributes to a file. Prepends the history string.
PetscErrorCode NCGlobalAttributes::write_attributes(const NCTool &nc, int, nc_type, bool) const {
  int ierr, ncid;
  string old_history;

  ierr = nc.get_att_text(NC_GLOBAL, "history", old_history);
  CHKERRQ(check_err(ierr,__LINE__,__FILE__));

  if (rank != 0) return 0;

  ncid = nc.get_ncid();

  ierr = nc.define_mode(); CHKERRQ(ierr); 

  // Write text attributes:
  map<string, string>::const_iterator i;
  for (i = strings.begin(); i != strings.end(); ++i) {
    string name  = i->first;
    string value = i->second;

    if (value.empty()) continue;

    if (name == "history") {
      // prepend:
      value = value + old_history;
    }      

    ierr = nc_put_att_text(ncid, NC_GLOBAL, name.c_str(), value.size(), value.c_str()); 
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  // Write double attributes:
  map<string, vector<double> >::const_iterator j;
  for (j = doubles.begin(); j != doubles.end(); ++j) {
    string name  = j->first;
    vector<double> values = j->second;

    if (values.empty()) continue;

    ierr = nc_put_att_double(ncid, NC_GLOBAL, name.c_str(), NC_DOUBLE, values.size(), &values[0]);
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  return 0;
}

//! Prepends \c message to the history string.
void NCGlobalAttributes::prepend_history(string message) {
  strings["history"] = message + strings["history"];

}
