// Copyright (C) 2007-2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cstring>
#include "udunits.h"
#include "pism_const.hh"
#include "NCTool.hh"

int nc_check(int stat) {
  if (stat)
    SETERRQ1(PETSC_COMM_SELF, 1, "NC_ERR: %s\n", nc_strerror(stat));
  return 0;
}

int check_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    (void) fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    SETERRQ1(PETSC_COMM_SELF, 1, "NC_ERR: %s\n", nc_strerror(stat));
    //exit(1);
  }
  return 0;
}

NCTool::NCTool(MPI_Comm c, PetscMPIInt r)
  : com(c), rank(r) {
  ncid = -1;
  def_mode = false;

  // Initialize UDUNITS if needed
  if (utIsInit() == 0) {
    if (utInit(NULL) != 0) {
      PetscPrintf(com, "PISM ERROR: UDUNITS initialization failed.\n");
      PISMEnd();
    }
  }
}

NCTool::~NCTool() {
  if (ncid >= 0) {
    PetscPrintf(com, "NCTool::~NCTool() ncid >= 0.\n");
    PISMEnd();
  }
}

//! Returns ncid corresponding to the current NetCDF file.
int NCTool::get_ncid() const {
  return ncid;
}

//! Finds a NetCDF dimension by its name.
PetscErrorCode  NCTool::find_dimension(string short_name, int *dimid, bool &exists) const {
  PetscErrorCode ierr;
  int stat, found = 0, my_dimid;
  if (rank == 0) {
    stat = nc_inq_dimid(ncid, short_name.c_str(), &my_dimid);
    if (stat == NC_NOERR)
      found = 1;
  }
  ierr = MPI_Bcast(&found, 1, MPI_INT, 0, com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&my_dimid, 1, MPI_INT, 0, com); CHKERRQ(ierr);


  if (found) {
    exists = true;
    if (dimid != NULL)
      *dimid = my_dimid;
  } else {
    exists = false;
    // dimid is not modified
  }

  return 0;
}

//! Finds the variable by its standard_name attribute.
/*!
  Here's how it works:
  <ol>
  <li> Check if the current IceModelVec has a standard_name. If it does, go to
  step 2, otherwise go to step 4.

  <li> Find the variable with this standard_name. Bail out if two
  variables have the same standard_name, otherwise go to step 3.

  <li> If no variable was found, go to step 4, otherwise go to step 5.

  <li> Find the variable with the right variable name. Go to step 5.

  <li> Broadcast the existence flag and the variable ID.
  </ol>
 */
PetscErrorCode  NCTool::find_variable(string short_name, string standard_name,
				      int *varidp, bool &exists,
				      bool &found_by_standard_name) const {
  int ierr;
  int stat, found = 0, my_varid = -1, nvars;
  bool standard_name_match = false;

  if (standard_name != "") {
    if (rank == 0) {
      ierr = nc_inq_nvars(ncid, &nvars); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    }
    ierr = MPI_Bcast(&nvars, 1, MPI_INT, 0, com); CHKERRQ(ierr);

    string attribute;
    for (int j = 0; j < nvars; j++) {
      ierr = get_att_text(j, "standard_name", attribute); CHKERRQ(ierr);
      if (attribute == "")
	continue;

      if (attribute == standard_name) {
	if (!found) {
	  found = true;
	  standard_name_match = true;
	  my_varid = j;
	} else {
	  ierr = PetscPrintf(com,
			     "PISM ERROR: Inconsistency in the input file: "
			     "Variables #%d and #%d have the same standard_name ('%s').\n",
			     my_varid, j, attribute.c_str());
	  CHKERRQ(ierr);
	  PISMEnd();
	}
      }
    } // end of for (int j = 0; j < nvars; j++)
  } // end of if (standard_name != "")

  // Check the short name:
  if (!found) {
    if (rank == 0) {
      stat = nc_inq_varid(ncid, short_name.c_str(), &my_varid);
      if (stat == NC_NOERR)
	found = true;
    }
    ierr = MPI_Bcast(&found,    1, MPI_INT, 0, com); CHKERRQ(ierr);
    ierr = MPI_Bcast(&my_varid, 1, MPI_INT, 0, com); CHKERRQ(ierr);
  }

  if (found) {
    exists = true;
    found_by_standard_name = standard_name_match;
    if (varidp != NULL)
      *varidp = my_varid;
  } else {
    exists = false;
    // *varidp is not modified
  }

  return 0;
}

//! \brief Find a variable and discard the found_by_standard_name flag.
PetscErrorCode NCTool::find_variable(string short_name, string standard_name,
				     int *varidp, bool &exists) const {
  bool dummy;
  PetscErrorCode ierr = find_variable(short_name, standard_name, varidp, exists, dummy);
  CHKERRQ(ierr);
  return 0;
}
				     
//! \brief Find a variable without specifying its standard name.
PetscErrorCode NCTool::find_variable(string short_name, int *varid, bool &exists) const {
  return find_variable(short_name, "", varid, exists);
}

//! Read in the variables \c z and \c zb from the NetCDF file; <i>do not</i> assume they are equally-spaced.
/*!
  This function allocates arrays z_levels and zb_levels, and they have to be
  freed by the caller (using delete[]).
 */

//! Put the variable for a dimension in a NetCDF file.  Makes no assumption about spacing.
PetscErrorCode NCTool::put_dimension(int varid, const vector<double> &vals) const {
  PetscErrorCode ierr;
  int stat;

  if (rank != 0) return 0;

  ierr = data_mode(); CHKERRQ(ierr); 

  stat = nc_put_var_double(ncid, varid, &vals[0]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  return 0;
}

//! Writes \c history to the history global attribute of a NetCDF dataset.
/*!
  Appends if overwrite == false (default).
 */
PetscErrorCode NCTool::write_history(string history, bool overwrite) const {
  int stat;
  string old_history, new_history;

  // Produce the new history string:
  stat = get_att_text(NC_GLOBAL, "history", old_history); CHKERRQ(stat);
  
  if (overwrite) {
    new_history = history;
  } else {
    new_history = history + old_history;
  }

  // Write it:
  if (rank == 0) {
    stat = define_mode(); CHKERRQ(stat); 
    
    stat = nc_put_att_text(ncid, NC_GLOBAL, "history",
			   new_history.size(), new_history.c_str());
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  return 0;
}

//! Checks if the dimension dim in a NetCDF file is OK.
/*! A dimension is OK if is exists and its length is equal to len. If len < 0,
    then dimension length is ignored.

    On processor 0 returns true if OK, false otherwise. Always returns true on
    processors other than 0.
 */
bool NCTool::check_dimension(string name, const int len) const {
  int stat, dimid;
  size_t dimlen;

  if (rank == 0) {
    stat = nc_inq_dimid(ncid, name.c_str(), &dimid);
    if (stat == NC_NOERR) {
      if (len < 0)
	return true;

      stat = nc_inq_dimlen(ncid, dimid, &dimlen);
      if (stat != NC_NOERR)
	return false;

      if ((int)dimlen == len)
	return true;
    }

    return false;
  }

  return true;
}


//! Appends \c time to the t dimension.
/*!
 * Does nothing on processors other than 0.
 */
PetscErrorCode NCTool::append_time(string dimension_name, PetscReal time) const {
  int stat, t_id;

  if (rank != 0) return 0;

  size_t t_len;
  stat = nc_inq_dimid(ncid, dimension_name.c_str(), &t_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  stat = nc_inq_dimlen(ncid, t_id, &t_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));

  stat = nc_inq_varid(ncid, dimension_name.c_str(), &t_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));

  stat = data_mode(); CHKERRQ(stat); 

  stat = nc_put_var1_double(ncid, t_id, &t_len, &time); CHKERRQ(check_err(stat,__LINE__,__FILE__));

  return 0;
}

PetscErrorCode NCTool::append_time_bounds(string name, PetscReal t0, PetscReal t1) const {
  int stat, t_id;

  if (rank != 0) return 0;

  size_t t_len;
  stat = nc_inq_dimid(ncid, name.c_str(), &t_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  stat = nc_inq_dimlen(ncid, t_id, &t_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));

  stat = nc_inq_varid(ncid, name.c_str(), &t_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));

  stat = data_mode(); CHKERRQ(stat); 

  size_t start[] = {t_len, 0}, count[] = {1, 2};
  double data[] = {t0, t1};
  stat = nc_put_vara_double(ncid, t_id, start, count, data); CHKERRQ(check_err(stat,__LINE__,__FILE__));

  return 0;
}

//! Opens a NetCDF file for reading.
/*!
  Stops if a file does not exist or could not be opened.
 */
PetscErrorCode NCTool::open_for_reading(string filename) {
  PetscErrorCode ierr;
  int stat = 0;

  if (ncid >= 0) SETERRQ(com, 1, "NCTool::open_for_reading(): ncid >= 0 at the beginning of the call");

  if (rank == 0) {
    stat = nc_open(filename.c_str(), NC_NOWRITE, &ncid);
  }
  ierr = MPI_Bcast(&stat, 1, MPI_INT, 0, com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, com); CHKERRQ(ierr);

  if (stat != NC_NOERR) {
    ierr = PetscPrintf(com, "ERROR: Can't open file '%s'!\n",
		       filename.c_str()); CHKERRQ(ierr);
    PISMEnd();
  }
  
  return 0;
}

//! Closes a NetCDF file.
PetscErrorCode NCTool::close() {
  PetscErrorCode ierr;
  if (rank == 0) {
    ierr = nc_close(ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }
  ncid = -1;			// make it invalid
  def_mode = false;
  return 0;
}

//! Opens a file for writing if it exists, creates if it does not.
PetscErrorCode NCTool::open_for_writing(string filename) {
  int stat;

  if (ncid >= 0) SETERRQ(com, 1, "NCTool::open_for_writing(): ncid >= 0 at the beginning of the call");

  if (rank == 0) {
    bool file_exists = false;

    // Check if the file exists:
    if (FILE *f = fopen(filename.c_str(), "r")) {
      file_exists = true;
      fclose(f);
    } else {
      file_exists = false;
    }

    if (file_exists) {
      stat = nc_open(filename.c_str(), NC_WRITE, &ncid);
      if (stat != NC_NOERR) {
	stat = PetscPrintf(com, "PISM ERROR: Can't open file '%s'. NetCDF error: %s\n",
			   filename.c_str(), nc_strerror(stat)); CHKERRQ(stat);
	PISMEnd();
      }
      stat = nc_set_fill(ncid, NC_NOFILL, NULL); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    } else {

#if (PISM_WRITE_NETCDF4 == 1)
      stat = nc_create(filename.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid); 
#else
      stat = nc_create(filename.c_str(), NC_CLOBBER|NC_64BIT_OFFSET, &ncid); 
#endif

      def_mode = true;
      CHKERRQ(check_err(stat,__LINE__,__FILE__));
      stat = nc_set_fill(ncid, NC_NOFILL, NULL); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    }

  } // end of if (rank == 0)

  stat = MPI_Bcast(&ncid, 1, MPI_INT, 0, com); CHKERRQ(stat);

  return 0;
}

//! Finds the length of a dimension. Returns 0 if failed.
PetscErrorCode NCTool::get_dim_length(string name, unsigned int *len) const {
  int stat, dim_id;
  
  if (rank == 0) {
    size_t dim_len;
    stat = nc_inq_dimid(ncid, name.c_str(), &dim_id);
    if (stat == NC_NOERR) {
      stat = nc_inq_dimlen(ncid, dim_id, &dim_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    } else
      dim_len = 0;

    *len = static_cast<int>(dim_len);
  }

  stat = MPI_Bcast(len, 1, MPI_INT, 0, com); CHKERRQ(stat);

  return 0;
}

//! \brief Get the number of records in a file.
PetscErrorCode NCTool::get_nrecords(int &nrecords) const {
  PetscErrorCode stat;

  int dimid;

  if (rank == 0) {
    size_t dim_len;
    stat = nc_inq_unlimdim(ncid, &dimid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    if (dimid == -1) {
      nrecords = 1;
    } else {
      stat = nc_inq_dimlen(ncid, dimid, &dim_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));
      nrecords = static_cast<int>(dim_len);
    }
  }

  stat = MPI_Bcast(&nrecords, 1, MPI_INT, 0, com); CHKERRQ(stat);

  return 0;
}


//! \brief Gets dimension limits.
/*! Gets dimension limits (%i.e. min and max values of the coordinate variable).

  Sets min = 0 and max = 0 if the dimension \c name has length 0.

  Set \c min or \c max to NULL to omit the corresponding value.

  Converts time and distance units to 'seconds' and 'meters', correspondingly.
 */
PetscErrorCode NCTool::get_dim_limits(string name, double *min, double *max) const {
  PetscErrorCode ierr;
  unsigned int len = 0;
  int varid = -1;
  bool variable_exists = false;
  size_t start = 0, count;
  double range[2] = {0, 0};

  ierr = get_dim_length(name, &len); CHKERRQ(ierr);

  if (len != 0) {
    ierr = find_variable(name, &varid, variable_exists);
    if (!variable_exists) {
      ierr = PetscPrintf(com, "PISM ERROR: coordinate variable '%s' does not exist.\n",
			 name.c_str());
      CHKERRQ(ierr);
      PISMEnd();
    }

    ierr = data_mode(); CHKERRQ(ierr); 

    if (rank == 0) {
      double *data;
      data = new double[len];

      count = static_cast<size_t>(len);
      ierr = nc_get_vara_double(ncid, varid, &start, &count, data); CHKERRQ(check_err(ierr,__LINE__,__FILE__));

      range[0] = data[0];
      range[1] = data[0];
      for (unsigned int j = 1; j < len; j++) {
	range[0] = PetscMin(data[j], range[0]);
	range[1] = PetscMax(data[j], range[1]);
      }
      delete[] data;
    } // end of if(rank == 0)
    ierr = MPI_Bcast(range, 2, MPI_DOUBLE, 0, com); CHKERRQ(ierr);

    char internal_units[TEMPORARY_STRING_LENGTH];
    utUnit input, internal;
    bool input_has_units;

    AxisType dimtype;
    ierr = inq_dimtype(name, dimtype); CHKERRQ(ierr);

    switch (dimtype) {
    case T_AXIS:
      // Note that this units specification does *not* have a reference date.
      strcpy(internal_units, "seconds");
      break;
    case UNKNOWN_AXIS:
      {
        ierr = verbPrintf(3, com,
                          "PISM WARNING: Can't determine which direction dimension '%s' corresponds to.\n"
                          "     Assuming that it corresponds to a spatial (not temporal) coordinate.\n",
                          name.c_str()); CHKERRQ(ierr);
      }
    case X_AXIS:
    case Y_AXIS:
    case Z_AXIS:
      strcpy(internal_units, "meters");
    }

    if (utScan(internal_units, &internal) != 0) {
      SETERRQ(com, 1, "UDUNITS failed trying to scan internal units.");
    }

    // Get the units information:
    ierr = get_units(varid, input_has_units, input); CHKERRQ(ierr);
    if (!input_has_units) {
      ierr = verbPrintf(3, com,
			"PISM WARNING: dimensional variable '%s' does not have the units attribute.\n"
			"     Assuming that it is in '%s'.\n",
			name.c_str(), internal_units); CHKERRQ(ierr);
      utCopy(&internal, &input);
    }

    // Find the conversion coefficients:
    double slope, intercept;
    ierr = utConvert(&input, &internal, &slope, &intercept);
    if (ierr != 0) {
      if (ierr == UT_ECONVERT) {
	ierr = PetscPrintf(com,
			   "PISM ERROR: dimensional variable '%s' has the units that are not compatible with '%s'.\n",
			   name.c_str(), internal_units); CHKERRQ(ierr);
	PISMEnd();
      }
      SETERRQ1(com, 1, "UDUNITS failure: error code = %d", ierr);
    }


    // Change units and return limits:
    if (min != NULL)
      *min = intercept + range[0]*slope;
    if (max != NULL)
      *max = intercept + range[1]*slope;
    
    return 0;
  } // if (len != 0)

  if (min != NULL) *min = 0.0;
  if (max != NULL) *max = 0.0;

  return 0;
}

//! \brief Reads a coordinate variable.
PetscErrorCode NCTool::get_dimension(string name, vector<double> &result) const {
  PetscErrorCode ierr;
  unsigned int len = 0;
  int varid = -1;
  bool variable_exists = false;
  size_t start = 0, count;

  ierr = get_dim_length(name, &len); CHKERRQ(ierr);

  if (len != 0) {
    ierr = find_variable(name, &varid, variable_exists);
    if (!variable_exists) {
      ierr = PetscPrintf(com, "PISM ERROR: coordinate variable '%s' does not exist.\n",
			 name.c_str());
      CHKERRQ(ierr);
      PISMEnd();
    }

    ierr = data_mode(); CHKERRQ(ierr); 

    double *data = new double[len];

    if (rank == 0) {
      count = static_cast<size_t>(len);
      ierr = nc_get_vara_double(ncid, varid, &start, &count, data); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    }

    ierr = MPI_Bcast(data, len, MPI_DOUBLE, 0, com); CHKERRQ(ierr);

    char internal_units[TEMPORARY_STRING_LENGTH];
    utUnit input, internal;
    bool input_has_units;

    AxisType dimtype;
    ierr = inq_dimtype(name, dimtype); CHKERRQ(ierr);

    switch (dimtype) {
    case T_AXIS:
      // Note that this units specification does *not* have a reference date.
      strcpy(internal_units, "seconds");
      break;
    case UNKNOWN_AXIS:
      {
        ierr = verbPrintf(3, com,
                          "PISM WARNING: Can't determine which direction dimension '%s' corresponds to.\n"
                          "     Assuming that it corresponds to a spatial (not temporal) coordinate.\n",
                          name.c_str()); CHKERRQ(ierr);
      }
    case X_AXIS:
    case Y_AXIS:
    case Z_AXIS:
      strcpy(internal_units, "meters");
    }

    if (utScan(internal_units, &internal) != 0) {
      SETERRQ(com, 1, "UDUNITS failed trying to scan internal units.");
    }

    // Get the units information:
    ierr = get_units(varid, input_has_units, input); CHKERRQ(ierr);
    if (!input_has_units) {
      ierr = verbPrintf(3, com,
			"PISM WARNING: dimensional variable '%s' does not have the units attribute.\n"
			"     Assuming that it is in '%s'.\n",
			name.c_str(), internal_units); CHKERRQ(ierr);
      utCopy(&internal, &input);
    }

    // Find the conversion coefficients:
    double slope, intercept;
    ierr = utConvert(&input, &internal, &slope, &intercept);
    if (ierr != 0) {
      if (ierr == UT_ECONVERT) {
	ierr = PetscPrintf(com,
			   "PISM ERROR: dimensional variable '%s' has the units that are not compatible with '%s'.\n",
			   name.c_str(), internal_units); CHKERRQ(ierr);
	PISMEnd();
      }
      SETERRQ1(com, 1, "UDUNITS failure: error code = %d", ierr);
    }

    result.resize(len);
    for (unsigned int i = 0; i < len; ++i)
      result[i] = intercept + data[i]*slope;

    delete [] data;
    
    return 0;
  } // if (len != 0)

  result.resize(0);

  return 0;
}

//! Reads a text attribute from a NetCDF file.
/*!
  Missing and empty attributes are treated the same.
 */
PetscErrorCode NCTool::get_att_text(const int varid, string name, string &result) const {
  char *str = NULL;
  int ierr, stat, len;

  // Read and broadcast the attribute length:
  if (rank == 0) {
    size_t attlen;
    stat = nc_inq_attlen(ncid, varid, name.c_str(), &attlen);
    if (stat == NC_NOERR)
      len = static_cast<int>(attlen);
    else
      len = 0;
  }
  ierr = MPI_Bcast(&len, 1, MPI_INT, 0, com); CHKERRQ(ierr);

  // Allocate some memory or set result to NULL and return:
  if (len == 0) {
    result = "";
    return 0;
  }
  str = new char[len + 1];
  // Zealously clear the string, so that we don't risk moving unitialized bytes
  // over MPI (because Valgrind can't tell the difference between these
  // harmless bytes and potential memory errors)
  ierr = PetscMemzero(str, len+1);CHKERRQ(ierr);

  // Now read the string and broadcast stat to see if we succeeded:
  if (rank == 0) {
    stat = nc_get_att_text(ncid, varid, name.c_str(), str);
  }
  ierr = MPI_Bcast(&stat, 1, MPI_INT, 0, com); CHKERRQ(ierr);
  
  // On success, broadcast the string. On error, set str to "".
  if (stat == NC_NOERR) {
    stat = MPI_Bcast(str, len, MPI_CHAR, 0, com); CHKERRQ(stat);
  } else {
    strcpy(str, "");
  }

  result = str;

  delete[] str;
  return 0;
}

//! Reads a scalar attribute from a NetCDF file.
PetscErrorCode NCTool::get_att_double(const int varid, string name,
				      vector<double> &result) const {
  int ierr, stat, len;

  // Read and broadcast the attribute length:
  if (rank == 0) {
    size_t attlen;
    stat = nc_inq_attlen(ncid, varid, name.c_str(), &attlen);
    if (stat == NC_NOERR)
      len = static_cast<int>(attlen);
    else
      len = 0;
  }
  ierr = MPI_Bcast(&len, 1, MPI_INT, 0, com); CHKERRQ(ierr);

  if (len == 0) {
    result.clear();
    return 0;
  }

  result.resize(len);
  // Now read the data and broadcast stat to see if we succeeded:
  if (rank == 0) {
    stat = nc_get_att_double(ncid, varid, name.c_str(), &result[0]);
  }
  ierr = MPI_Bcast(&stat, 1, MPI_INT, 0, com); CHKERRQ(ierr);
  
  // On success, broadcast the data. On error, stop.
  if (stat == NC_NOERR) {
    ierr = MPI_Bcast(&result[0], len, MPI_DOUBLE, 0, com); CHKERRQ(ierr);
  } else {
    SETERRQ(com, 1, "Error reading an attribute.");
  }

  return 0;
}

//! Get variable's units information from a NetCDF file.
/*! Note that this function intentionally ignores the reference date ('years
  since 1-1-1', 'years since 2000-1-1' and 'years' produce the same result).
 */
PetscErrorCode NCTool::get_units(int varid, bool &has_units, utUnit &units) const {
  PetscErrorCode ierr;
  string units_string;

  // Get the string:
  ierr = get_att_text(varid, "units", units_string); CHKERRQ(ierr);

  // If a variables does not have the units attribute, set the flag and return:
  if (units_string.empty()) {
    has_units = false;
    utClear(&units);
    return 0;
  }

  /*!
    \note This method finds the string "since" in the units_string and
    terminates it on the first 's' of "since", if this sub-string was found.
    This is done to ignore the reference date in the time units string (the
    reference date specification always starts with this word).
  */
  int n = (int)units_string.find("since");
  if (n != -1) units_string.resize(n);

  ierr = utScan(units_string.c_str(), &units);
  if (ierr != 0) {
    ierr = PetscPrintf(com, "PISM ERROR: units specification '%s' is unknown or invalid.\n",
		       units_string.c_str());
    PISMEnd();
  }
  
  has_units = true;

  return 0;
}

//! Get the number of attributes of a variable.
PetscErrorCode NCTool::inq_nattrs(int varid, int &N) const {
  int stat;

  if (rank == 0) {
    stat = nc_inq_varnatts(ncid, varid, &N); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  stat = MPI_Bcast(&N, 1, MPI_INT, 0, com); CHKERRQ(stat);

  return 0;
}

//! Get the attribute type.
PetscErrorCode NCTool::inq_att_type(int varid, string name, nc_type &result) const {
  int stat, type;
  nc_type tmp;

  if (rank == 0) {
    stat = nc_inq_atttype(ncid, varid, name.c_str(), &tmp); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    type = static_cast<int>(tmp);
  } // end of if(rank == 0)
  stat = MPI_Bcast(&type, 1, MPI_INT, 0, com); CHKERRQ(stat);

  result = static_cast<nc_type>(type);

  return 0;
}

//! Gets the name of the n-th (counting from 0) attribute of a NetCDF variable.
PetscErrorCode NCTool::inq_att_name(int varid, int n, string &name) const {
  int stat;
  char tmp[NC_MAX_NAME];

  if (rank == 0) {
    stat = nc_inq_attname(ncid, varid, n, tmp); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  stat = MPI_Bcast(tmp, NC_MAX_NAME, MPI_CHAR, 0, com); CHKERRQ(stat);

  name = tmp;
  return 0;
}

//! Gets the list of dimension ids for dimensions a variable depends on.
/*!
  The length of the result (\c dimids) is the number of dimensions.
 */
PetscErrorCode NCTool::inq_dimids(int varid, vector<int> &dimids) const {
  int stat;
  int ndims;

  if (rank == 0) {
    stat = nc_inq_varndims(ncid, varid, &ndims); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  stat = MPI_Bcast(&ndims, 1, MPI_INT, 0, com); CHKERRQ(stat);

  if (ndims == 0) {
    dimids.clear();
    return 0;
  }

  dimids.resize(ndims);		// every processor allocates at least ndims
				// integers (if necessary)

  if (rank == 0) {
    stat = nc_inq_vardimid(ncid, varid, &dimids[0]); CHKERRQ(check_err(stat,__LINE__,__FILE__)); 
  }
  stat = MPI_Bcast(&dimids[0], ndims, MPI_INT, 0, com); CHKERRQ(stat);

  return 0;
}

//! Get a name of a dimension by a dimension ID.
PetscErrorCode NCTool::inq_dimname(int dimid, string &name) const {
  int stat;
  char tmp[NC_MAX_NAME];

  if (rank == 0) {
    stat = nc_inq_dimname(ncid, dimid, tmp); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  stat = MPI_Bcast(tmp, NC_MAX_NAME, MPI_CHAR, 0, com); CHKERRQ(stat);

  name = tmp;
  return 0;
}

//! Get the dimension ID of an unlimited dimension. Sets unlimdimid to -1 if there isn't one.
PetscErrorCode NCTool::inq_unlimdim(int &unlimdimid) const {
  int stat;

  if (rank == 0) {
    stat = nc_inq_unlimdim(ncid, &unlimdimid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }
  stat = MPI_Bcast(&unlimdimid, 1, MPI_INT, 0, com); CHKERRQ(stat);

  return 0;
}

//! Moves \c filename to \c filename~ if \c filename exists.
PetscErrorCode NCTool::move_if_exists(string filename) {
  PetscErrorCode ierr;

  if (rank != 0)
    return 0;

  // Check if the file exists:
  if (FILE *f = fopen(filename.c_str(), "r")) {
    fclose(f);
  } else {
    return 0;
  }
  
  string tmp = filename;
  tmp = tmp + "~";
      
  ierr = rename(filename.c_str(), tmp.c_str());
  if (ierr != 0) {
    ierr = verbPrintf(1, com, "PISM ERROR: can't move '%s' to '%s'.\n",
                      filename.c_str(), tmp.c_str());
    PISMEnd();
  }
  ierr = verbPrintf(2, com, 
                    "PISM WARNING: output file '%s' already exists. Moving it to '%s'.\n",
                    filename.c_str(), tmp.c_str());

  return 0;
}

//! \brief Puts a NetCDF file in define mode if it is not in define mode
//! already.
PetscErrorCode NCTool::define_mode() const {
  PetscErrorCode ierr;

  if (rank != 0) return 0;

  if (def_mode) return 0;
  
  ierr = nc_redef(ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));

  def_mode = true;

  return 0;
}

//! \brief Puts a NetCDF file in data mode if it is not in data mode already.
PetscErrorCode NCTool::data_mode() const {
  PetscErrorCode ierr;

  if (rank != 0) return 0;

  if (!def_mode) return 0;

  //! 50000 (below) means that we allocate 50Kb for metadata in NetCDF files
  //! created by PISM.
  ierr = nc__enddef(ncid,50000,4,0,4); CHKERRQ(check_err(ierr,__LINE__,__FILE__));

  def_mode = false;

  return 0;
}

PetscErrorCode NCTool::set_attrs(int varid, map<string,string> attrs) const {
  PetscErrorCode ierr;

  if (rank != 0) return 0;

  ierr = define_mode(); CHKERRQ(ierr);

  map<string,string>::iterator j = attrs.begin();
  while(j != attrs.end()) {
    string name = j->first,
      value = j->second;

    ierr = nc_put_att_text(ncid, varid, name.c_str(), value.size(), value.c_str());
    check_err(ierr,__LINE__,__FILE__);

    ++j;
  }

  return 0;
}

//! \brief Creates a dimension and the corresponding coordinate variable; sets
//! string sttributes.
PetscErrorCode NCTool::create_dimension(string name, int length, map<string,string> attrs,
                                        int &dimid, int &varid) const {
  PetscErrorCode ierr;

  if (rank != 0) return 0;

  ierr = define_mode(); CHKERRQ(ierr);

  ierr = nc_def_dim(ncid, name.c_str(), length, &dimid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  ierr = nc_def_var(ncid, name.c_str(), NC_DOUBLE, 1, &dimid, &varid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));

  ierr = set_attrs(varid, attrs); CHKERRQ(ierr);

  return 0;
}

//! \brief Determines what kind of dimension dimid refers to.
/*!
 * Checks units for compatibility with "seconds", then checks the
 * "standard_name" attribute, then the "axis" attribute, then the variable
 * name, which can be "x", "X" or starting with "x" or "X" (or "y", etc).
 */
PetscErrorCode NCTool::inq_dimtype(string name, AxisType &result) const {
  PetscErrorCode ierr;
  double slope, intercept;
  string axis, standard_name, units;
  utUnit tmp_units, ut_units;
  int varid;
  bool exists;

  ierr = find_variable(name, &varid, exists); CHKERRQ(ierr);

  if (!exists) {
    PetscPrintf(com, "ERROR: coordinate variable '%s' is not present!\n", name.c_str());
    PISMEnd();
  }

  ierr = get_att_text(varid, "axis", axis); CHKERRQ(ierr);
  ierr = get_att_text(varid, "standard_name", standard_name); CHKERRQ(ierr);
  ierr = get_att_text(varid, "units", units); CHKERRQ(ierr);

  // check the units attribute:
  bool has_time_units = false;

  // check if it has units compatible with "seconds":
  ierr = utScan(units.c_str(), &ut_units);
  if (ierr != 0) {
    ierr = PetscPrintf(com, "ERROR: units specification '%s' is unknown or invalid.\n",
		       units.c_str());
    PISMEnd();
  }

  utScan("seconds", &tmp_units);
  ierr = utConvert(&ut_units, &tmp_units, &slope, &intercept);
  if (ierr == 0) has_time_units = true;

  if (has_time_units) {
    result = T_AXIS;
    return 0;
  }

  // check the standard_name attribute:
  if (standard_name == "time") {
    result = T_AXIS;
    return 0;
  } else if (standard_name == "projection_x_coordinate") {
    result = X_AXIS;
    return 0;
  } else if (standard_name == "projection_y_coordinate") {
    result = Y_AXIS;
    return 0;
  }

  // check the axis attribute:
  if (axis == "T" || axis == "t") {
    result = T_AXIS;
    return 0;
  } else if (axis == "X" || axis == "x") {
    result = X_AXIS;
    return 0;
  } else if (axis == "Y" || axis == "y") {
    result = Y_AXIS;
    return 0;
  } else if (axis == "Z" || axis == "z") {
    result = Z_AXIS;
    return 0;
  }

  // check the variable name:
  if (name == "x" || name == "X" ||
      name.find("x") == 0 || name.find("X") == 0) {
    result = X_AXIS;
    return 0;
  } else if (name == "y" || name == "Y" ||
             name.find("y") == 0 || name.find("Y") == 0) {
    result = Y_AXIS;
    return 0;
  } else if (name == "z" || name == "Z" ||
             name.find("z") == 0 || name.find("Z") == 0) {
    result = Z_AXIS;
    return 0;
  } else if (name == "t" || name == "T" || name == "time" ||
             name.find("t") == 0 || name.find("T") == 0) {
    result = T_AXIS;
    return 0;
  }

  // we have no clue:
  result = UNKNOWN_AXIS;
  return 0;
}

//! \brief Get the number of records of variable \c varname.
PetscErrorCode NCTool::get_nrecords(string short_name, string standard_name,
                                    unsigned int &nrecords) const {
  PetscErrorCode ierr;

  bool exists;
  int varid;
  vector<int> dimids;

  ierr = find_variable(short_name, standard_name, &varid, exists); CHKERRQ(ierr);
  
  if (!exists) {
    nrecords = 0;
    return 0;
  }

  ierr = inq_dimids(varid, dimids); CHKERRQ(ierr);

  for (unsigned int i = 0; i < dimids.size(); ++i) {
    AxisType dimtype;
    string dimname;
    
    ierr = inq_dimname(dimids[i], dimname); CHKERRQ(ierr);
    ierr = inq_dimtype(dimname, dimtype); CHKERRQ(ierr);

    if (dimtype == T_AXIS) {
      ierr = get_dim_length(dimname, &nrecords); CHKERRQ(ierr);
      return 0;
    }
  }

  nrecords = 1;
  return 0;
}

