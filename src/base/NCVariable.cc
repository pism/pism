// Copyright (C) 2009 Constantine Khroulev and Ed Bueler
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
#include "pism_const.hh"

NCVariable::NCVariable() {
  grid = NULL;
  dims = GRID_2D;
  reset();
}

//! Initialize a NCVariable instance.
void NCVariable::init(string name, IceGrid &g, GridType d) {
  short_name = name;
  grid = &g;
  dims = d;
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

//! Read a variable from a file into a \b global Vec v.
/*! This also converts the data from input units to internal units if needed.
 */
PetscErrorCode NCVariable::read(const char filename[], unsigned int time, Vec v) {
  PetscErrorCode ierr;
  bool variable_exists;
  int varid;
  NCTool nc(grid);

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
    ierr = PetscPrintf(grid->com,
		      "PISM ERROR: Can't find '%s' (%s) in '%s'.\n",
		       short_name.c_str(),
		       strings["standard_name"].c_str(), filename);
    CHKERRQ(ierr);
    PetscEnd();
  }

  ierr = nc.get_global_var(varid, v, dims, time); CHKERRQ(ierr);  

  bool input_has_units;
  utUnit input_units;

  ierr = nc.get_units(varid,
		      input_has_units, input_units); CHKERRQ(ierr);

  if ( has("units") && (!input_has_units) ) {
    string &units_string = strings["units"],
      &long_name = strings["long_name"];
    ierr = verbPrintf(2, grid->com,
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
PetscErrorCode NCVariable::write(const char filename[], nc_type nctype,
				 bool write_in_glaciological_units, Vec v) {
  PetscErrorCode ierr;
  bool exists;
  NCTool nc(grid);
  int varid;

  if (grid == NULL)
    SETERRQ(1, "NCVariable::write: grid is NULL.");

  ierr = nc.open_for_writing(filename, false); CHKERRQ(ierr); // replace = false, because
				// we want to *append* at this point

  // find or define the variable
  ierr = nc.find_variable(short_name, strings["standard_name"],
			  &varid, exists); CHKERRQ(ierr);

  if (!exists) {
    ierr = define(nc.ncid, nctype, varid); CHKERRQ(ierr);
  }

  if (write_in_glaciological_units) {
    ierr = change_units(v, &units, &glaciological_units); CHKERRQ(ierr);
  }

  // write the attributes
  write_attributes(nc.ncid, varid, nctype, write_in_glaciological_units);

  // Actually write data:
  ierr = nc.put_global_var(varid, v, dims); CHKERRQ(ierr);  
  
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
PetscErrorCode NCVariable::regrid(const char filename[], LocalInterpCtx &lic,
				  bool critical, bool set_default_value,
				  PetscScalar default_value,
				  MaskInterp * interpolation_mask,
				  Vec v) {
  int varid;
  bool exists;
  PetscErrorCode ierr;
  NCTool nc(grid);

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
      ierr = PetscPrintf(grid->com,
			"PISM ERROR: Can't find '%s' in the regridding file '%s'.\n",
			 short_name.c_str(), filename);
      CHKERRQ(ierr);
      PetscEnd();
    }

    if (set_default_value) {	// if it's not and we have a default value, set it
      double slope, intercept, tmp;
      utConvert(&units, &glaciological_units, &slope, &intercept);
      tmp = intercept + slope*default_value;
      
      ierr = verbPrintf(2, grid->com, 
			"  absent %-10s/ %-60s\n   %-16s\\ not found; using default constant %7.2f (%s)\n",
			short_name.c_str(),
			strings["long_name"].c_str(),
			"", tmp,
			strings["glaciological_units"].c_str());
      CHKERRQ(ierr);
      ierr = VecSet(v, default_value); CHKERRQ(ierr);
    } else {			// otherwise leave it alone
      ierr = verbPrintf(2, grid->com, 
			"  absent %-10s/ %-60s\n   %-16s\\ not found; continuing without setting it\n",
			short_name.c_str(),
			strings["long_name"].c_str(), "");
      CHKERRQ(ierr);
    }
  } else {			// the variable was found successfully
    // Check if it is discrete
    bool use_interpolation_mask = (interpolation_mask != NULL);
    if (use_interpolation_mask)
      nc.set_MaskInterp(interpolation_mask);

    ierr = nc.regrid_global_var(varid, dims, lic, v,
				use_interpolation_mask); CHKERRQ(ierr);

    // Now we need to get the units string from the file and convert the units,
    // because check_range and report_range expect the data to be in PISM (SI)
    // units.

    bool input_has_units;
    utUnit input_units;

    ierr = nc.get_units(varid, input_has_units, input_units); CHKERRQ(ierr);

    if ( has("units") && (!input_has_units) ) {
      ierr = verbPrintf(2, grid->com,
			"PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
			"              Assuming that it is in '%s'.\n",
			strings["short_name"].c_str(),
			strings["long_name"].c_str(),
			strings["units"].c_str()); CHKERRQ(ierr);
      utCopy(&units, &input_units);
    }

    // Convert data:
    ierr = change_units(v, &input_units, &units); CHKERRQ(ierr);

    // Read the valid range info:
    ierr = read_valid_range(nc.ncid, varid); CHKERRQ(ierr);

    // Check the range and warn the user if needed:
    ierr = check_range(v); CHKERRQ(ierr);

    // We can report the success, and the range now:
    ierr = verbPrintf(2, grid->com, "  FOUND ");
    ierr = report_range(v, found_by_standard_name); CHKERRQ(ierr);
  } // end of if(exists)

  ierr = nc.close(); CHKERRQ(ierr);
  return 0;
}

//! Read the valid range information from a file.
/*! Reads \c valid_min, \c valid_max and \c valid_range attributes; if \c
    valid_range is found, sets the pair \c valid_min and \c valid_max instead.
 */
PetscErrorCode NCVariable::read_valid_range(int ncid, int varid) {
  NCTool nc(grid);
  string input_units_string;
  utUnit input_units;
  vector<double> bounds;
  double slope, intercept;
  int stat;

  if (grid == NULL)
    SETERRQ(1, "NCVariable::read_valid_range: grid is NULL.");

  // Never reset valid_min/max if any of them was set internally.
  if (has("valid_min") || has("valid_max"))
    return 0;

  nc.ncid = ncid;

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

//! Converts \c v from the units corresponding to \c from to the ones corresponding to \c to.
/*!
  Does nothing if this transformation is trivial.
 */
PetscErrorCode NCVariable::change_units(Vec v, utUnit *from, utUnit *to) {
  PetscErrorCode ierr;
  double slope, intercept;
  char from_name[TEMPORARY_STRING_LENGTH], to_name[TEMPORARY_STRING_LENGTH], *tmp;
  bool use_slope, use_intercept;

  // Get string representations of units:
  utPrint(from, &tmp);
  strncpy(from_name, tmp, TEMPORARY_STRING_LENGTH);
  utPrint(to, &tmp);
  strncpy(to_name, tmp, TEMPORARY_STRING_LENGTH);

  // Get the slope and the intercept of the linear transformation.
  ierr = utConvert(from, to, &slope, &intercept);

  if (ierr != 0) { 		// can't convert
    if (ierr == UT_ECONVERT) {	// because units are incompatible
      ierr = verbPrintf(2, grid->com,
			"PISM ERROR: IceModelVec '%s': attempted to convert data from '%s' to '%s'.\n",
			short_name.c_str(), from_name, to_name);
      return 1;
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
PetscErrorCode NCVariable::write_attributes(int ncid, int varid, nc_type nctype,
					    bool write_in_glaciological_units) {
  int ierr;

  if (grid == NULL)
    SETERRQ(1, "NCVariable::write_attributes: grid is NULL.");

  if (grid->rank != 0) return 0;

  ierr = nc_redef(ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));

  // units, valid_min, valid_max and valid_range need special treatment:
  if (has("units")) {
    string output_units = strings["units"];
    if (write_in_glaciological_units)
      output_units = strings["glaciological_units"];
    ierr = nc_put_att_text(ncid, varid,
			   "units",
			   output_units.size(), output_units.c_str());
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  double bounds[2];
  // We need to save valid_min, valid_max and valid_range in the units
  // matching the ones in the output.
  if (write_in_glaciological_units) {
    double slope, intercept;

    ierr = utConvert(&units, &glaciological_units, &slope, &intercept); CHKERRQ(ierr);

    bounds[0] = intercept + slope*get("valid_min");
    bounds[1] = intercept + slope*get("valid_max");
  } else {
    bounds[0] = get("valid_min");
    bounds[1] = get("valid_max");
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
  map<string, string>::iterator i;
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
  map<string, vector<double> >::iterator j;
  for (j = doubles.begin(); j != doubles.end(); ++j) {
    string name  = j->first;
    vector<double> values = j->second;

    if (name == "valid_min" ||
	name == "valid_max" ||
	name == "valid_range" ||
	values.empty())
      continue;

    ierr = nc_put_att_double(ncid, varid, name.c_str(), nctype, values.size(), &values[0]);
    CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  ierr = nc_enddef(ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));

  return 0;
}


//! Report the range of a \b global Vec \c v.
PetscErrorCode NCVariable::report_range(Vec v, bool found_by_standard_name) {
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
      ierr = verbPrintf(2, grid->com, 
			" %-10s/ standard_name=%-60s\n   %-16s\\ min,max = %9.3f,%9.3f (%s)\n",
			short_name.c_str(),
			strings["standard_name"].c_str(), "", min, max,
			strings["glaciological_units"].c_str()); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid->com, 
			" %-10s/ WARNING! standard_name=%s is missing, found by short_name\n   %-16s\\ min,max = %9.3f,%9.3f (%s)\n",
			short_name.c_str(),
			strings["standard_name"].c_str(), "", min, max,
			strings["glaciological_units"].c_str()); CHKERRQ(ierr);
    }

  } else {

    ierr = verbPrintf(2, grid->com, 
		      " %-10s/ %-60s\n   %-16s\\ min,max = %9.3f,%9.3f (%s)\n",
		      short_name.c_str(),
		      strings["long_name"].c_str(), "", min, max,
		      strings["glaciological_units"].c_str()); CHKERRQ(ierr);
  }

  return 0;
}

//! Check if the range of a \b global Vec \c v is in the range specified by valid_min and valid_max attributes.
PetscErrorCode NCVariable::check_range(Vec v) {
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
      ierr = verbPrintf(2, grid->com,
			"PISM WARNING: some values of '%s' are outside the valid range [%f, %f] (%s)\n",
			short_name.c_str(), valid_min, valid_max, units_string.c_str()); CHKERRQ(ierr);
    
  } else if (has("valid_min")) {
    double valid_min = get("valid_min");
    if (min < valid_min) {
      ierr = verbPrintf(2, grid->com,
			"PISM WARNING: some values of '%s' are less than the valid minimum %f (%s)\n",
			short_name.c_str(), valid_min, units_string.c_str()); CHKERRQ(ierr);
    }
    
  } else if (has("valid_max")) {
    double valid_max = get("valid_max");
    if (max > valid_max) {
      ierr = verbPrintf(2, grid->com,
			"PISM WARNING: some values of '%s' are greater than the valid maximum %f (%s)\n",
			short_name.c_str(), valid_max, units_string.c_str()); CHKERRQ(ierr);
    }
  }
  return 0;
}

//! Define a NetCDF variable corresponding to a NCVariable object.
PetscErrorCode NCVariable::define(int ncid, nc_type nctype, int &varid) {
  int stat, dimids[4], var_id;

  if (grid == NULL)
    SETERRQ(1, "NCVariable::define: grid is NULL.");

  if (grid->rank == 0) {
    stat = nc_redef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "t", &dimids[0]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "y", &dimids[1]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "x", &dimids[2]); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    switch (dims) {
    case GRID_2D:
      stat = nc_def_var(ncid, short_name.c_str(), nctype, 3, dimids, &var_id);
      break;
    case GRID_3D:
      stat = nc_inq_dimid(ncid, "z", &dimids[3]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
      stat = nc_def_var(ncid, short_name.c_str(), nctype, 4, dimids, &var_id);
      break;
    case GRID_3D_BEDROCK:
      stat = nc_inq_dimid(ncid, "zb", &dimids[3]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
      stat = nc_def_var(ncid, short_name.c_str(), nctype, 4, dimids, &var_id);
    }
    
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_enddef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  stat = MPI_Bcast(&var_id, 1, MPI_INT, 0, grid->com); CHKERRQ(stat);

  varid = var_id;

  return 0;
}

//! Reset all the attributes.
PetscErrorCode NCVariable::reset() {

  strings.clear();
  short_name = "unnamed_variable";
  // long_name is unset
  // standard_name is unset
  // pism_intent is unset;
  strings["coordinates"] = "lat lon";

  doubles.clear();
  // valid_min and valid_max are unset

  utClear(&units);
  utClear(&glaciological_units);

  return 0;
}

//! Checks if an attribute is present. Ignores empty strings.
bool NCVariable::has(string name) {

  if (strings.find(name) != strings.end()) {
    if (strings[name].empty())
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
double NCVariable::get(string name) {
  if (doubles.find(name) != doubles.end())
    return doubles[name][0];
  else
    return 0;
}
