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
#include "nc_util.hh"

int nc_check(int stat) {
  if (stat)
    SETERRQ1(1, "NC_ERR: %s\n", nc_strerror(stat));
  return 0;
}

int check_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    (void) fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    SETERRQ1(1, "NC_ERR: %s\n", nc_strerror(stat));
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
PetscErrorCode  NCTool::find_dimension(const char short_name[], int *dimid, bool &exists) const {
  PetscErrorCode ierr;
  int stat, found = 0, my_dimid;
  if (rank == 0) {
      stat = nc_inq_dimid(ncid, short_name, &my_dimid);
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

PetscErrorCode NCTool::find_variable(string short_name, string standard_name,
				     int *varidp, bool &exists) const {
  bool dummy;
  PetscErrorCode ierr = find_variable(short_name, standard_name, varidp, exists, dummy);
  CHKERRQ(ierr);
  return 0;
}
				     

PetscErrorCode NCTool::find_variable(string short_name, int *varid, bool &exists) const {
  return find_variable(short_name, "", varid, exists);
}

//! Read the last value of the time variable t from a NetCDF file.
PetscErrorCode NCTool::get_last_time(double *time) const {
  PetscErrorCode ierr;
  ierr = get_dim_limits("t", NULL, time); CHKERRQ(ierr);
  return 0;
}


//! Read in the variables \c z and \c zb from the NetCDF file; <i>do not</i> assume they are equally-spaced.
/*!
  This function allocates arrays z_levels and zb_levels, and they have to be
  freed by the caller (using delete[]).
 */
PetscErrorCode NCTool::get_vertical_dims(double* &z_levels, double* &zb_levels) const {
  int stat;
  int z_id, zb_id, z_len, zb_len;
  size_t zero  = 0, nc_z_len, nc_zb_len;

  stat = get_dim_length("z", &z_len); CHKERRQ(stat);
  stat = get_dim_length("zb", &zb_len); CHKERRQ(stat);

  z_levels = new double[z_len];
  zb_levels = new double[zb_len];

  nc_z_len  = (size_t) z_len;
  nc_zb_len = (size_t) zb_len;

  if (rank == 0) {
    stat = nc_inq_varid(ncid, "z", &z_id);
    if (stat != NC_NOERR) {
      stat = PetscPrintf(com, "PISM ERROR: Can't find the 'z' coordinate variable.\n"); CHKERRQ(stat);
      PISMEnd();
    }
    
    stat = nc_inq_varid(ncid, "zb", &zb_id);
    if (stat != NC_NOERR) {
      stat = PetscPrintf(com, "PISM ERROR: Can't find the 'zb' coordinate variable.\n"); CHKERRQ(stat);
      PISMEnd();
    }

    stat = data_mode(); CHKERRQ(stat); 

    stat = nc_get_vara_double(ncid, z_id, &zero, &nc_z_len, z_levels);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_get_vara_double(ncid, zb_id, &zero, &nc_zb_len, zb_levels);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  MPI_Bcast(z_levels, z_len, MPI_DOUBLE, 0, com);
  MPI_Bcast(zb_levels, zb_len, MPI_DOUBLE, 0, com);
  return 0;
}


//! Put the variable for a dimension in a NetCDF file.
/*! Uses starting and ending values and a grid length for regularly-spaced
values.
 */
PetscErrorCode NCTool::put_dimension_regular(int varid, int len, double start, double end) const {
  PetscErrorCode ierr;
  int stat;
  double *v, delta;

  if (rank != 0) return 0;

  if (end <= start)
    SETERRQ2(1, "Can't write dimension: start = %f, end = %f",
	     start, end);

  delta = (end - start) / (len - 1);

  ierr = PetscMalloc(len * sizeof(double), &v); CHKERRQ(ierr);
  for (int i = 0; i < len; i++)
    v[i] = start + i * delta;

  // Sometimes v[len - 1] turns out to be greater than end (because of rounding
  // errors). If that happens, we need to fix v[len - 1] by setting it equal to
  // end.
  if (v[len - 1] > end) v[len - 1] = end;

  ierr = data_mode(); CHKERRQ(ierr);
  
  stat = nc_put_var_double(ncid, varid, v); CHKERRQ(check_err(stat,__LINE__,__FILE__));

  ierr = PetscFree(v); CHKERRQ(ierr);

  return 0;
}


//! Put the variable for a dimension in a NetCDF file.  Makes no assumption about spacing.
PetscErrorCode NCTool::put_dimension(int varid, int len, PetscScalar *vals) const {
  PetscErrorCode ierr;
  int stat;
  double *v;

  if (rank != 0) return 0;

  ierr = PetscMalloc(len * sizeof(double), &v); CHKERRQ(ierr);
  for (int i = 0; i < len; i++) {
    v[i] = (double)vals[i];
  }

  ierr = data_mode(); CHKERRQ(ierr); 

  stat = nc_put_var_double(ncid, varid, v); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  ierr = PetscFree(v); CHKERRQ(ierr);
  return 0;
}

//! Writes \c history to the history global attribute of a NetCDF dataset.
/*!
  Appends if overwrite == false (default).
 */
PetscErrorCode NCTool::write_history(const char history[], bool overwrite) const {
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
bool NCTool::check_dimension(const char name[], const int len) const {
  int stat, dimid;
  size_t dimlen;

  if (rank == 0) {
    stat = nc_inq_dimid(ncid, name, &dimid);
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
PetscErrorCode NCTool::append_time(PetscReal time) const {
  int stat, t_id;


  if (rank == 0) {
    size_t t_len;
    stat = nc_inq_dimid(ncid, "t", &t_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimlen(ncid, t_id, &t_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_inq_varid(ncid, "t", &t_id); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = data_mode(); CHKERRQ(stat); 

    stat = nc_put_var1_double(ncid, t_id, &t_len, &time); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  return 0;
}

//! Opens a NetCDF file for reading.
/*!
  Stops if a file does not exist or could not be opened.
 */
PetscErrorCode NCTool::open_for_reading(const char filename[]) {
  PetscErrorCode ierr;
  int stat = 0;

  if (ncid >= 0) SETERRQ(1, "NCTool::open_for_reading(): ncid >= 0 at the beginning of the call");

  if (rank == 0) {
    stat = nc_open(filename, NC_NOWRITE, &ncid);
  }
  ierr = MPI_Bcast(&stat, 1, MPI_INT, 0, com); CHKERRQ(ierr);
  ierr = MPI_Bcast(&ncid, 1, MPI_INT, 0, com); CHKERRQ(ierr);

  if (stat != NC_NOERR) {
    ierr = PetscPrintf(com, "ERROR: Can't open file '%s'!\n",
		       filename); CHKERRQ(ierr);
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
PetscErrorCode NCTool::open_for_writing(const char filename[]) {
  int stat;

  if (ncid >= 0) SETERRQ(1, "NCTool::open_for_writing(): ncid >= 0 at the beginning of the call");

  if (rank == 0) {
    bool file_exists = false;

    // Check if the file exists:
    if (FILE *f = fopen(filename, "r")) {
      file_exists = true;
      fclose(f);
    } else {
      file_exists = false;
    }

    if (file_exists) {
      stat = nc_open(filename, NC_WRITE, &ncid);
      if (stat != NC_NOERR) {
	stat = PetscPrintf(com, "PISM ERROR: Can't open file '%s'. NetCDF error: %s\n",
			   filename, nc_strerror(stat)); CHKERRQ(stat);
	PISMEnd();
      }
      stat = nc_set_fill(ncid, NC_NOFILL, NULL); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    } else {

#if (PISM_WRITE_NETCDF4 == 1)
      stat = nc_create(filename, NC_CLOBBER|NC_NETCDF4, &ncid); 
#else
      stat = nc_create(filename, NC_CLOBBER|NC_64BIT_OFFSET, &ncid); 
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
PetscErrorCode NCTool::get_dim_length(const char name[], int *len) const {
  int stat, dim_id;
  
  if (rank == 0) {
    size_t dim_len;
    stat = nc_inq_dimid(ncid, name, &dim_id);
    if (stat == NC_NOERR) {
      stat = nc_inq_dimlen(ncid, dim_id, &dim_len); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    } else
      dim_len = 0;

    *len = static_cast<int>(dim_len);
  }

  stat = MPI_Bcast(len, 1, MPI_INT, 0, com); CHKERRQ(stat);

  return 0;
}

//! Gets dimension limits.
/*! Gets dimension limits (%i.e. min and max values of the coordinate variable).

  Sets min = 0 and max = 0 if the dimension \c name has length 0.

  Set \c min or \c max to NULL to omit the corresponding value.
 */
PetscErrorCode NCTool::get_dim_limits(const char name[], double *min, double *max) const {
  PetscErrorCode ierr;
  int len = 0, varid = -1;
  bool variable_exists = false;
  size_t start = 0, count;
  double range[2] = {0, 0};

  ierr = get_dim_length(name, &len); CHKERRQ(ierr);

  if (len != 0) {
    ierr = find_variable(name, &varid, variable_exists);
    if (!variable_exists) {
      ierr = PetscPrintf(com, "PISM ERROR: coordinate variable '%s' does not exist.\n",
			 name);
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
      for (int j = 1; j < len; j++) {
	range[0] = PetscMin(data[j], range[0]);
	range[1] = PetscMax(data[j], range[1]);
      }
      delete[] data;
    } // end of if(rank == 0)
    ierr = MPI_Bcast(range, 2, MPI_DOUBLE, 0, com); CHKERRQ(ierr);

    char internal_units[TEMPORARY_STRING_LENGTH];
    utUnit input, internal;
    bool input_has_units;

    if (strcmp(name, "t") == 0) {
      // Note that this units specification does *not* have a reference date.
      strcpy(internal_units, "seconds");
    } else {
      strcpy(internal_units, "meters");
    }

    if (utScan(internal_units, &internal) != 0) {
      SETERRQ(1, "UDUNITS failed trying to scan internal units.");
    }

    // Get the units information:
    ierr = get_units(varid, input_has_units, input); CHKERRQ(ierr);
    if (!input_has_units) {
      ierr = verbPrintf(3, com,
			"PISM WARNING: dimensional variable '%s' does not have the units attribute.\n"
			"     Assuming that it is in '%s'.\n",
			name, internal_units); CHKERRQ(ierr);
      utCopy(&internal, &input);
    }

    // Find the conversion coefficients:
    double slope, intercept;
    ierr = utConvert(&input, &internal, &slope, &intercept);
    if (ierr != 0) {
      if (ierr == UT_ECONVERT) {
	ierr = PetscPrintf(com,
			   "PISM ERROR: dimensional variable '%s' has the units that are not compatible with '%s'.\n",
			   name, internal_units); CHKERRQ(ierr);
	PISMEnd();
      }
      SETERRQ1(1, "UDUNITS failure: error code = %d", ierr);
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

PetscErrorCode NCTool::get_dimension(const char name[], vector<double> &result) const {
  PetscErrorCode ierr;
  int len = 0, varid = -1;
  bool variable_exists = false;
  size_t start = 0, count;

  ierr = get_dim_length(name, &len); CHKERRQ(ierr);

  if (len != 0) {
    ierr = find_variable(name, &varid, variable_exists);
    if (!variable_exists) {
      ierr = PetscPrintf(com, "PISM ERROR: coordinate variable '%s' does not exist.\n",
			 name);
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

    if (strcmp(name, "t") == 0) {
      // Note that this units specification does *not* have a reference date.
      strcpy(internal_units, "seconds");
    } else {
      strcpy(internal_units, "meters");
    }

    if (utScan(internal_units, &internal) != 0) {
      SETERRQ(1, "UDUNITS failed trying to scan internal units.");
    }

    // Get the units information:
    ierr = get_units(varid, input_has_units, input); CHKERRQ(ierr);
    if (!input_has_units) {
      ierr = verbPrintf(3, com,
			"PISM WARNING: dimensional variable '%s' does not have the units attribute.\n"
			"     Assuming that it is in '%s'.\n",
			name, internal_units); CHKERRQ(ierr);
      utCopy(&internal, &input);
    }

    // Find the conversion coefficients:
    double slope, intercept;
    ierr = utConvert(&input, &internal, &slope, &intercept);
    if (ierr != 0) {
      if (ierr == UT_ECONVERT) {
	ierr = PetscPrintf(com,
			   "PISM ERROR: dimensional variable '%s' has the units that are not compatible with '%s'.\n",
			   name, internal_units); CHKERRQ(ierr);
	PISMEnd();
      }
      SETERRQ1(1, "UDUNITS failure: error code = %d", ierr);
    }

    result.resize(len);
    for (int i = 0; i < len; ++i)
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
PetscErrorCode NCTool::get_att_text(const int varid, const char name[], string &result) const {
  char *str = NULL;
  int ierr, stat, len;

  // Read and broadcast the attribute length:
  if (rank == 0) {
    size_t attlen;
    stat = nc_inq_attlen(ncid, varid, name, &attlen);
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
    stat = nc_get_att_text(ncid, varid, name, str);
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
PetscErrorCode NCTool::get_att_double(const int varid, const char name[],
				      vector<double> &result) const {
  int ierr, stat, len;

  // Read and broadcast the attribute length:
  if (rank == 0) {
    size_t attlen;
    stat = nc_inq_attlen(ncid, varid, name, &attlen);
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
    stat = nc_get_att_double(ncid, varid, name, &result[0]);
  }
  ierr = MPI_Bcast(&stat, 1, MPI_INT, 0, com); CHKERRQ(ierr);
  
  // On success, broadcast the data. On error, stop.
  if (stat == NC_NOERR) {
    ierr = MPI_Bcast(&result[0], len, MPI_DOUBLE, 0, com); CHKERRQ(ierr);
  } else {
    SETERRQ(1, "Error reading an attribute.");
  }

  return 0;
}

//! Get variable's units information from a NetCDF file.
/*!
  Note that this function intentionally ignores the reference date.
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

  // This finds the string "since" in the units_string and terminates
  // it on the first 's' of "since", if this sub-string was found. This
  // is done to ignore the reference date in the time units string (the
  // reference date specification always starts with this word).

  int n = units_string.find("since");
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
PetscErrorCode NCTool::inq_att_type(int varid, const char name[], nc_type &result) const {
  int stat, type;
  nc_type tmp;

  if (rank == 0) {
    stat = nc_inq_atttype(ncid, varid, name, &tmp); CHKERRQ(check_err(stat,__LINE__,__FILE__));
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

//! Gets the list of dimensions a variable depends on.
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
PetscErrorCode NCTool::move_if_exists(const char filename[]) {
  PetscErrorCode ierr;

  if (rank != 0)
    return 0;

  // Check if the file exists:
  if (FILE *f = fopen(filename, "r")) {
    fclose(f);
  } else {
    return 0;
  }
  
  string tmp = filename;
  tmp = tmp + "~";
      
  ierr = rename(filename, tmp.c_str());
  if (ierr != 0) {
    ierr = verbPrintf(1, com, "PISM ERROR: can't move '%s' to '%s'.\n",
                      filename, tmp.c_str());
    PISMEnd();
  }
  ierr = verbPrintf(2, com, 
                    "PISM WARNING: output file '%s' already exists. Moving it to '%s'.\n",
                    filename, tmp.c_str());

  return 0;
}

PetscErrorCode NCTool::define_mode() const {
  PetscErrorCode ierr;

  if (rank != 0) return 0;

  if (def_mode) return 0;
  
  ierr = nc_redef(ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));

  def_mode = true;

  return 0;
}

PetscErrorCode NCTool::data_mode() const {
  PetscErrorCode ierr;

  if (rank != 0) return 0;

  if (!def_mode) return 0;

  ierr = nc__enddef(ncid,50000,4,0,4); CHKERRQ(check_err(ierr,__LINE__,__FILE__))

  def_mode = false;

  return 0;
}

