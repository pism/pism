// Copyright (C) 2008--2009 Ed Bueler and Constantine Khroulev
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

#include "pism_const.hh"
#include "iceModelVec.hh"


IceModelVec::IceModelVec() {

  shallow_copy = false;
  v = PETSC_NULL;
  da = PETSC_NULL;
  grid = PETSC_NULL;
  array = PETSC_NULL;
  localp = true;
  use_interpolation_mask = false;

  strcpy(short_name,"*****UNKNOWN***** variable name");

  reset_attrs();
}


//! Creates a shallow copy of an \c IceModelVec.
/*!
  No data is copied to the new IceModelVec.

  The difference is that such a copy will not free the memory when deleted (or
  goes out of scope).
 */
IceModelVec::IceModelVec(const IceModelVec &o) {
  // This IceModelVec is a shallow copy!
  shallow_copy = true;

  v = o.v;
  da = o.da;
  dims = o.dims;
  grid = o.grid;
  array = o.array;
  localp = o.localp;

  use_interpolation_mask = o.use_interpolation_mask;
  interpolation_mask = o.interpolation_mask;

  strcpy(short_name, o.short_name);

  has_long_name = o.has_long_name;
  strcpy(long_name, o.long_name);

  has_units = o.has_units;
  strcpy(units_string, o.units_string);
  strcpy(glaciological_units_string, o.glaciological_units_string);

  has_standard_name = o.has_standard_name;
  strcpy(standard_name, o.standard_name);

  has_pism_intent = o.has_pism_intent;
  strcpy(pism_intent, o.pism_intent);

  has_coordinates = o.has_coordinates;
  strcpy(coordinates, o.coordinates);

  has_valid_min = o.has_valid_min;
  has_valid_max = o.has_valid_max;
  valid_min = o.valid_min;
  valid_max = o.valid_max;

  units = o.units;
  glaciological_units = o.glaciological_units;
}


IceModelVec::~IceModelVec() {
  // Only destroy the IceModelVec if it is not a shallow copy:
  if (!shallow_copy) destroy();
}

PetscErrorCode  IceModelVec::create(IceGrid &/*mygrid*/, const char /*my_short_name*/[], 
                                    bool /*local*/) {
  SETERRQ(1,"IceModelVec::create(...) is VIRTUAL ONLY: not implemented");
}


//! Returns true if create() was called and false otherwise.
bool IceModelVec::was_created() {
  return (v != PETSC_NULL);
}

//! Returns the grid type of an IceModelVec. (This is the way to figure out if an IceModelVec is 2D or 3D).
GridType IceModelVec::grid_type() {
  return dims;
}

PetscErrorCode  IceModelVec::destroy() {
  PetscErrorCode ierr;
  if (v != PETSC_NULL) {
    ierr = VecDestroy(v); CHKERRQ(ierr);
    v = PETSC_NULL;
  }
  if (da != PETSC_NULL) {
    ierr = DADestroy(da); CHKERRQ(ierr);
    da = PETSC_NULL;
  }
  return 0;
}


PetscErrorCode  IceModelVec::printInfo(const PetscInt verbosity) {
  PetscErrorCode ierr;
  
  if (grid == PETSC_NULL) {
    SETERRQ1(1,"ERROR: cannot print info for IceModelVec with short_name='%s'\n"
               "   because grid=PETSC_NULL.  ENDING.\n\n",short_name);
  }

  ierr = verbPrintf(verbosity,grid->com,
         "\nprinting info for IceModelVec with short_name='%s':\n",
         short_name); CHKERRQ(ierr);
  if (da == PETSC_NULL) {
    ierr = verbPrintf(verbosity,grid->com,
          "  WARNING:  da == PETSC_NULL for IceModelVec with short_name='%s'!\n",
          short_name); CHKERRQ(ierr);
  }
  if (v == PETSC_NULL) {
    ierr = verbPrintf(verbosity,grid->com,
          "  WARNING:  v == PETSC_NULL for IceModelVec with short_name='%s'!\n",
          short_name); CHKERRQ(ierr);
  }
  if (array == PETSC_NULL) {
    ierr = verbPrintf(verbosity,grid->com,
          "  WARNING:  array == PETSC_NULL for IceModelVec with short_name='%s'!\n",
          short_name); CHKERRQ(ierr);
  }
  
  ierr = verbPrintf(verbosity,grid->com,
           "  boolean flags:  localp = %d, has_standard_name = %d\n",
		    (int)localp, (int)has_standard_name);  CHKERRQ(ierr);

  ierr = verbPrintf(verbosity,grid->com,
           "                  long_name = '%s'\n", long_name);  CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,grid->com,
           "                  standard_name = '%s'\n", standard_name);  CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,grid->com,
           "                  units = '%s'\n", units_string);  CHKERRQ(ierr);
  ierr = verbPrintf(verbosity,grid->com,
           "                  pism_intent = '%s'\n\n", pism_intent);  CHKERRQ(ierr);
  return 0;
}

//! Result: min <- min(v[j]), max <- max(v[j]).
/*! 
PETSc manual correctly says "VecMin and VecMax are collective on Vec" but
GlobalMax,GlobalMin \e are needed, when localp==true, to get correct 
values because Vecs created with DACreateLocalVector() are of type 
VECSEQ and not VECMPI.  See src/trypetsc/localVecMax.c.
 */
PetscErrorCode IceModelVec::range(PetscReal &min, PetscReal &max) {
  PetscReal my_min, my_max, gmin, gmax;
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);

  ierr = VecMin(v, PETSC_NULL, &my_min); CHKERRQ(ierr);
  ierr = VecMax(v, PETSC_NULL, &my_max); CHKERRQ(ierr);

  if (localp) {
    // needs a reduce operation; use PetscGlobalMax;
    ierr = PetscGlobalMin(&my_min, &gmin, grid->com); CHKERRQ(ierr);
    ierr = PetscGlobalMax(&my_max, &gmax, grid->com); CHKERRQ(ierr);
    min = gmin;
    max = gmax;
  } else {
    min = my_min;
    max = my_max;
  }
  return 0;
}


//! Computes the norm of an IceModelVec by calling PETSc VecNorm.
/*! 
See comment for range(); because local Vecs are VECSEQ, needs a reduce operation.
See src/trypetsc/localVecMax.c.
 */
PetscErrorCode IceModelVec::norm(NormType n, PetscReal &out) {
  PetscErrorCode ierr;
  PetscReal my_norm, gnorm;
  ierr = checkAllocated(); CHKERRQ(ierr);

  ierr = VecNorm(v, n, &my_norm); CHKERRQ(ierr);

  if (localp) {
    // needs a reduce operation; use PetscGlobalMax if NORM_INFINITY,
    //   otherwise PetscGlobalSum; carefully in NORM_2 case
    if (n == NORM_1_AND_2) {
      SETERRQ1(1, 
         "IceModelVec::norm(...): NORM_1_AND_2 not implemented (called as %s.norm(...))\n",
         short_name);
    } else if (n == NORM_1) {
      ierr = PetscGlobalSum(&my_norm, &gnorm, grid->com); CHKERRQ(ierr);
    } else if (n == NORM_2) {
      my_norm = PetscSqr(my_norm);  // undo sqrt in VecNorm before sum
      ierr = PetscGlobalSum(&my_norm, &gnorm, grid->com); CHKERRQ(ierr);
      gnorm = sqrt(gnorm);
    } else if (n == NORM_INFINITY) {
      ierr = PetscGlobalMax(&my_norm, &gnorm, grid->com); CHKERRQ(ierr);
    } else {
      SETERRQ1(2, "IceModelVec::norm(...): unknown NormType (called as %s.norm(...))\n",
         short_name);
    }
    out = gnorm;
  } else {
    out = my_norm;
  }
  return 0;
}


//! Result: v <- sqrt(v), elementwise.  Calls VecSqrt(v).
/*!
Name avoids clash with sqrt() in math.h.
 */
PetscErrorCode IceModelVec::squareroot() {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);

  ierr = VecSqrt(v); CHKERRQ(ierr);
  return 0;
}


//! Result: v <- v + alpha * x. Calls VecAXPY.
PetscErrorCode IceModelVec::add(PetscScalar alpha, IceModelVec &x) {
  PetscErrorCode ierr;
  PetscInt X_size, Y_size;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = x.checkAllocated(); CHKERRQ(ierr);

  if (dims != x.dims) {
    SETERRQ(1, "IceModelVec::add(...): operands have different numbers of dimensions");
  }

  ierr = VecGetSize(v, &X_size); CHKERRQ(ierr);
  ierr = VecGetSize(x.v, &Y_size); CHKERRQ(ierr);
  if (X_size != Y_size)
    SETERRQ1(1, "IceModelVec::add(...): incompatible Vec sizes (called as %s.add(...))\n", short_name);

  ierr = VecAXPY(v, alpha, x.v); CHKERRQ(ierr);
  return 0;
}

//! Result: result <- v + alpha * x. Calls VecWAXPY.
PetscErrorCode IceModelVec::add(PetscScalar alpha, IceModelVec &x, IceModelVec &result) {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = x.checkAllocated(); CHKERRQ(ierr);
  ierr = result.checkAllocated(); CHKERRQ(ierr);

  if ((dims != x.dims) || (dims != result.dims)) {
    SETERRQ(1, "IceModelVec::add(...): operands have different numbers of dimensions");
  }

  ierr = VecWAXPY(result.v, alpha, x.v, v); CHKERRQ(ierr);
  return 0;
}

//! Result: v[j] <- v[j] + alpha for all j. Calls VecShift.
PetscErrorCode IceModelVec::shift(PetscScalar alpha) {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);

  ierr = VecShift(v, alpha); CHKERRQ(ierr);
  return 0;
}

//! Result: v <- v * alpha. Calls VecScale.
PetscErrorCode IceModelVec::scale(PetscScalar alpha) {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);

  ierr = VecScale(v, alpha); CHKERRQ(ierr);
  return 0;
}

//! Result: result <- v .* x.  Calls VecPointwiseMult.
PetscErrorCode  IceModelVec::multiply_by(IceModelVec &x, IceModelVec &result) {
  PetscErrorCode ierr;
  PetscInt X_size, Y_size;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = x.checkAllocated(); CHKERRQ(ierr);

  if ((dims != x.dims) ||
      (dims != result.dims)) {
    SETERRQ(1, "IceModelVec::multiply_by(...): operands have different numbers of dimensions");
  }

  ierr = VecGetSize(v, &X_size); CHKERRQ(ierr);
  ierr = VecGetSize(x.v, &Y_size); CHKERRQ(ierr);

  if (X_size != Y_size)
    SETERRQ1(1, "IceModelVec::multiply_by(...): incompatible Vec sizes (called as %s.multiply_by(...))\n", short_name);

  ierr = VecPointwiseMult(result.v, v, x.v); CHKERRQ(ierr);
  return 0;
}

//! Result: v <- v .* x.  Calls VecPointwiseMult.
PetscErrorCode  IceModelVec::multiply_by(IceModelVec &x) {
  PetscErrorCode ierr;
  PetscInt X_size, Y_size;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = x.checkAllocated(); CHKERRQ(ierr);

  if ((dims != x.dims)) {
    SETERRQ(1, "IceModelVec::multiply_by(...): operands have different numbers of dimensions");
  }

  ierr = VecGetSize(v, &X_size); CHKERRQ(ierr);
  ierr = VecGetSize(x.v, &Y_size); CHKERRQ(ierr);

  if (X_size != Y_size)
    SETERRQ1(1, "IceModelVec::multiply_by(...): incompatible Vec sizes (called as %s.multiply_by(...))\n", short_name);

  ierr = VecPointwiseMult(v, v, x.v); CHKERRQ(ierr);
  return 0;
}

//! Copies v to a global vector 'destination'. Ghost points are discarded.
/*! This is potentially dangerous: make sure that \c destination has the same
    dimensions as the current IceModelVec.
 */
PetscErrorCode  IceModelVec::copy_to_global(Vec destination) {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);

  if (!localp)
    SETERRQ2(1, "Use copy_to(...). (Called %s.copy_to_global(...) and %s is global)", 
    short_name, short_name);

  ierr = DALocalToGlobal(da, v, INSERT_VALUES, destination); CHKERRQ(ierr);
  return 0;
}

//! Result: destination <- v.  Leaves metadata alone but copies values in Vec.  Uses VecCopy.
PetscErrorCode  IceModelVec::copy_to(IceModelVec &destination) {
  PetscErrorCode ierr;
  PetscInt X_size, Y_size;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = destination.checkAllocated(); CHKERRQ(ierr);

  if ((dims != destination.dims)) {
    SETERRQ(1, "IceModelVec::copy_to(...): operands have different numbers of dimensions");
  }

  ierr = VecGetSize(v, &X_size); CHKERRQ(ierr);
  ierr = VecGetSize(destination.v, &Y_size); CHKERRQ(ierr);

  if (X_size != Y_size)
    SETERRQ1(1, "IceModelVec::copy_to(...): incompatible Vec sizes (called as %s.copy_to(...))\n", short_name);

  ierr = VecCopy(v, destination.v); CHKERRQ(ierr);
  return 0;
}

//! Result: v <- source.  Leaves metadata alone but copies values in Vec.  Uses VecCopy.
PetscErrorCode  IceModelVec::copy_from(IceModelVec &source) {
  PetscErrorCode ierr;
  PetscInt X_size, Y_size;
  ierr = checkAllocated(); CHKERRQ(ierr);

  if ((dims != source.dims)) {
    SETERRQ(1, "IceModelVec::copy_from(...): operands have different numbers of dimensions");
  }
  
  ierr = VecGetSize(source.v, &X_size); CHKERRQ(ierr);
  ierr = VecGetSize(v, &Y_size); CHKERRQ(ierr);

  if (X_size != Y_size)
    SETERRQ1(1, "IceModelVec::copy_from(...): incompatible Vec sizes (called as %s.copy_from(...))\n", short_name);

  ierr = VecCopy(source.v, v); CHKERRQ(ierr);
  return 0;
}

//! Sets the variable name to \c name.
PetscErrorCode  IceModelVec::set_name(const char name[]) {
  reset_attrs();
  strcpy(short_name, name);
  return 0;
}

//! Sets the glaciological units of an IceModelVec.
/*!
This affects IceModelVec::report_range() and IceModelVec::write().  In write(),
if IceModelVec::write_in_glaciological_units == true, then that variable is written
with a conversion to the glaciological units set here.
 */
PetscErrorCode  IceModelVec::set_glaciological_units(const char my_units[]) {
  double a, b;			// dummy variables
  if (utScan(my_units, &glaciological_units) != 0) {
    SETERRQ2(1, "PISM ERROR: IceModelVec '%s': unknown or invalid units specification '%s'.", short_name, my_units);
  }
  
  if (utConvert(&units, &glaciological_units, &a, &b) == UT_ECONVERT) {
    SETERRQ3(1, "PISM ERROR: IceModelVec '%s': attempted to set glaciological units to '%s', which is not compatible with '%s'.\n",
	     short_name, my_units, units_string);
  }

  // Save the human-friendly version of the string; this is to avoid getting
  // things like '3.16887646408185e-08 meter second-1' instead of 'm year-1'
  // (and thus violating the CF conventions).
  strncpy(glaciological_units_string, my_units, PETSC_MAX_PATH_LEN);
  return 0;
}

//! Resets most IceModelVec attributes.
PetscErrorCode IceModelVec::reset_attrs() {

  write_in_glaciological_units = false;

  strcpy(long_name,                  "unknown long name");
  has_long_name = false;

  strcpy(units_string,               "unitless");
  strcpy(glaciological_units_string, "unitless");
  has_units = false;

  strcpy(pism_intent,                "unknown pism_intent");
  has_pism_intent = false;

  strcpy(standard_name,"unknown CF standard_name");
  has_standard_name = false;

  strcpy(coordinates, "lat lon");
  has_coordinates = true;	// true by default

  has_valid_min = false;
  has_valid_max = false;
  valid_min = 0;
  valid_max = 0;

  utClear(&units);
  utClear(&glaciological_units);

  return 0;
}

//! Sets NetCDF attributes of an IceModelVec object.
/*! Call set_attrs("new long name", "new units", "new pism_intent", NULL) if a
  variable does not have a standard name. Similarly, by putting NULL in an
  appropriate spot, it is possible tp leave long_name, units or pism_intent
  unmodified.

  If my_units != NULL, this also resets glaciological_units, so that they match
  internal units.
 */
PetscErrorCode  IceModelVec::set_attrs(const char my_pism_intent[],
				       const char my_long_name[], const char my_units[],
				       const char my_standard_name[]) {

  if (my_long_name != NULL) {
    has_long_name = true;
    strcpy(long_name,my_long_name);
  } else {
    has_long_name = false;
  }

  if (my_units != NULL) {
    has_units = true;
    strcpy(units_string, my_units);

    if (utScan(my_units, &units) != 0) {
      SETERRQ2(1, "PISM ERROR: IceModelVec '%s': unknown or invalid units specification '%s'.", short_name, my_units);
    }

    // Set the glaciological units too:
    utCopy(&units, &glaciological_units);
    strncpy(glaciological_units_string, my_units, PETSC_MAX_PATH_LEN);
  } else {
    has_units = false;
  }

  if (my_pism_intent != NULL) {
    has_pism_intent = true;
    strcpy(pism_intent,my_pism_intent);
  } else {
    has_pism_intent = false;
  }

  if (my_standard_name != NULL) {
    strcpy(standard_name,my_standard_name);
    has_standard_name = true;
  } else {
    has_standard_name = false;
  }
  return 0;
}


//! Defines a netcdf variable corresponding to an IceModelVec object.
PetscErrorCode IceModelVec::define_netcdf_variable(int ncid, nc_type nctype, int *varidp) {
  int stat, dimids[4], var_id;

  if (grid->rank == 0) {
    stat = nc_redef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "t", &dimids[0]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "y", &dimids[1]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
    stat = nc_inq_dimid(ncid, "x", &dimids[2]); CHKERRQ(check_err(stat,__LINE__,__FILE__));

    switch (dims) {
    case GRID_2D:
      stat = nc_def_var(ncid, short_name, nctype, 3, dimids, &var_id);
      break;
    case GRID_3D:
      stat = nc_inq_dimid(ncid, "z", &dimids[3]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
      stat = nc_def_var(ncid, short_name, nctype, 4, dimids, &var_id);
      break;
    case GRID_3D_BEDROCK:
      stat = nc_inq_dimid(ncid, "zb", &dimids[3]); CHKERRQ(check_err(stat,__LINE__,__FILE__));
      stat = nc_def_var(ncid, short_name, nctype, 4, dimids, &var_id);
    }
    
    CHKERRQ(check_err(stat,__LINE__,__FILE__));

    stat = nc_enddef(ncid); CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  stat = MPI_Bcast(&var_id, 1, MPI_INT, 0, grid->com); CHKERRQ(stat);

  *varidp = var_id;

  return 0;
}

//! Writes NetCDF attributes to a dataset.
/*! Call this <b>after</b> making sure that the NetCDF variable is defined.
 */
PetscErrorCode IceModelVec::write_attrs(const int ncid, nc_type nctype) {
  bool exists;
  int varid, ierr;
  NCTool nc(grid);
  nc.ncid = ncid;

  ierr = nc.find_variable(short_name, standard_name, &varid, exists); CHKERRQ(ierr);
  if (!exists) {
    SETERRQ(1, "Can't write attributes of an undefined variable.");
  }

  if (grid->rank == 0) {
    ierr = nc_redef(ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));

    if (has_pism_intent) {
      ierr = nc_put_att_text(ncid, varid,"pism_intent", strlen(pism_intent), pism_intent);
      CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    }

    if (has_units) {
      char *output_units = units_string;
      if (write_in_glaciological_units)
	output_units = glaciological_units_string;
      ierr = nc_put_att_text(ncid, varid,"units", strlen(output_units), output_units);
      CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    }

    if (has_long_name) {
      ierr = nc_put_att_text(ncid, varid,"long_name", strlen(long_name), long_name);
      CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    }

    if (has_standard_name) {
      ierr = nc_put_att_text(ncid, varid,"standard_name", strlen(standard_name), standard_name);
      CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    }

    if (has_coordinates) {
      ierr = nc_put_att_text(ncid, varid,"coordinates", strlen(coordinates), coordinates);
      CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    }

    double bounds[2];
    // We need to save valid_min, valid_max and valid_range in the units
    // matching the ones in the output.
    if (write_in_glaciological_units) {
      double slope, intercept;

      ierr = utConvert(&units, &glaciological_units, &slope, &intercept); CHKERRQ(ierr);

      bounds[0] = intercept + slope*valid_min;
      bounds[1] = intercept + slope*valid_max;
    } else {
      bounds[0] = valid_min;
      bounds[1] = valid_max;
    }

    if (has_valid_min && has_valid_max) {
      ierr = nc_put_att_double(ncid, varid, "valid_range", nctype, 2, bounds);
      CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    } else if (has_valid_min) {
      ierr = nc_put_att_double(ncid, varid, "valid_min", nctype, 1, &bounds[0]);
      CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    } else if (has_valid_max) {
      ierr = nc_put_att_double(ncid, varid, "valid_max", nctype, 1, &bounds[1]);
      CHKERRQ(check_err(ierr,__LINE__,__FILE__));
    }

    ierr = nc_enddef(ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  return 0;
}

//! Gets an IceModelVec from a file \c filename, interpolating onto the current grid.
/*! Stops if the variable was not found and \c critical == true.
 */
PetscErrorCode IceModelVec::regrid(const char filename[], LocalInterpCtx &lic, 
                                   bool critical) {
  PetscErrorCode ierr;
  // Signature:
  // regrid_from_netcdf(filename, lic, critical, set_default_value, default_value)
  ierr = regrid_from_netcdf(filename, lic, critical, false, 0.0); CHKERRQ(ierr);
  return 0;
}

//! Gets an IceModelVec from a file \c filename, interpolating onto the current grid.
/*! Sets all the values to \c default_value if the variable was not found..
 */
PetscErrorCode IceModelVec::regrid(const char filename[], 
                                   LocalInterpCtx &lic, PetscScalar default_value) {
  PetscErrorCode ierr;
  // Signature:
  // regrid_from_netcdf(filename, lic, critical, set_default_value, default_value)
  ierr = regrid_from_netcdf(filename, lic, false, true, default_value); CHKERRQ(ierr);
  return 0;
}

//! Calls the appropriate NCTool method to read a NetCDF variable into the IceModelVec.
/*!
  <ol>
  <li> Opens the file by calling NCTool::open_for_reading(...)
  <li> Finds the variable by calling NCTool::find_variable(...)
  <li> Reads data by calling NCTool::get_global_var(...) or NCTool::get_local_var(...)
  </ol>
 */
PetscErrorCode IceModelVec::read(const char filename[], const unsigned int time) {           
  PetscErrorCode ierr;
  bool variable_exists;
  int varid;
  NCTool nc(grid);

  ierr = checkAllocated(); CHKERRQ(ierr);
  if (grid->da2 == PETSC_NULL)
    SETERRQ(1, "IceModelVec::read_from_netcdf: grid.da2 is NULL.");

  // Open the file:
  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
  
  // Find the variable:
  ierr = nc.find_variable(short_name, standard_name, &varid, variable_exists); CHKERRQ(ierr);
  if (!variable_exists) {
    ierr = PetscPrintf(grid->com,
		      "PISM ERROR: Can't find '%s' (%s) in '%s'.\n",
		       short_name, standard_name, filename);
    CHKERRQ(ierr);
    PetscEnd();
  }

  if (localp) {
    ierr = nc.get_local_var(varid, da, v, dims, time); CHKERRQ(ierr);  
  } else {
    ierr = nc.get_global_var(varid, v, dims, time); CHKERRQ(ierr);  
  }

  bool input_has_units;
  utUnit input_units;

  ierr = nc.get_units(short_name, standard_name,
		      input_has_units, input_units); CHKERRQ(ierr);

  if ( has_units && (!input_has_units) ) {
    ierr = verbPrintf(2, grid->com,
		      "PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
		      "              Assuming that it is in '%s'.\n",
		      short_name, long_name, units_string); CHKERRQ(ierr);
    utCopy(&units, &input_units);
  }

  // Convert data:
  ierr = change_units(&input_units, &units); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);
  return 0;
}

//! Writes an IceModelVec to a NetCDF file.
/*!
  1) Get the last time.
  2) Find the variable in the file. Call define_variable if not found.
  3) Call put_global_var or put_local_var.
 */
PetscErrorCode IceModelVec::write(const char filename[], nc_type nctype) {
  PetscErrorCode ierr;
  bool exists;
  NCTool nc(grid);
  int varid;

  ierr = checkAllocated(); CHKERRQ(ierr);
  
  ierr = nc.open_for_writing(filename, false); CHKERRQ(ierr); // replace = false, because
				// we want to *append* at this point

  // find or define the variable
  ierr = nc.find_variable(short_name, standard_name, &varid, exists); CHKERRQ(ierr);
  if (!exists) {
    ierr = define_netcdf_variable(nc.ncid, nctype, &varid); CHKERRQ(ierr);
  }

  if (write_in_glaciological_units) {
    ierr = change_units(&units, &glaciological_units); CHKERRQ(ierr);
  }

  // write the attributes
  write_attrs(nc.ncid, nctype);

  // Actually write data:
  if (localp) {
    ierr = nc.put_local_var(varid, da, v, dims); CHKERRQ(ierr);  
  } else {
    ierr = nc.put_global_var(varid, v, dims); CHKERRQ(ierr);  
  }
  
  if (write_in_glaciological_units) {
    ierr = change_units(&glaciological_units, &units); CHKERRQ(ierr); // restore the units
  }

  ierr = nc.close(); CHKERRQ(ierr);
  return 0;
}

//! Reports the range of an IceModelVec, with the appropriate units.
/*! Always uses the glaciological units, since they match the internal units if
    not set.
 */
PetscErrorCode IceModelVec::report_range() {
  double slope, intercept;
  PetscErrorCode ierr;
  PetscReal min, max;

  // Get the conversion coefficients:
  utConvert(&units, &glaciological_units, &slope, &intercept);

  ierr = range(min, max);
  // Note that in some cases the following conversion does nothing
  min = min * slope + intercept;
  max = max * slope + intercept;

  ierr = verbPrintf(2, grid->com, 
		    " %-10s/ %-60s\n   %-16s\\ min,max = %9.3f,%9.3f (%s)\n",
		    short_name, long_name, "", min, max,
		    glaciological_units_string); CHKERRQ(ierr);

  return 0;
}

//! Reads data from a NetCDF file, with regridding.
/*!
  <ol>
  <li> Open the file
  <li> Find the variable
  <li> If \c critical == true, break if it was not found.
  <li> If \c critical == false, regrid the variable or set the default value if asked to.
  </ol>
 */
PetscErrorCode IceModelVec::regrid_from_netcdf(const char filename[],
					       LocalInterpCtx &lic, bool critical,
					       bool set_default_value,
					       PetscScalar default_value) {
  int varid;
  bool exists;
  PetscErrorCode ierr;
  NCTool nc(grid);

  ierr = checkAllocated(); CHKERRQ(ierr);
  if (grid->da2 == PETSC_NULL)
    SETERRQ(1, "IceModelVec::regrid_from_netcdf: grid.da2 is NULL.");

  // Open the file
  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);

  // Find the variable
  ierr = nc.find_variable(short_name, standard_name, &varid, exists); CHKERRQ(ierr);

  if (!exists) {		// couldn't find the variable
    if (critical) {		// if it's critical, print an error message and stop
      // SETERRQ1(1, "Variable '%s' was not found.\n", short_name);
      ierr = PetscPrintf(grid->com,
			"PISM ERROR: Can't find '%s' in the regridding file '%s'.\n",
			short_name, filename);
      CHKERRQ(ierr);
      PetscEnd();
    }

    if (set_default_value) {	// if it's not and we have a default value, set it
      double slope, intercept, tmp;
      utConvert(&units, &glaciological_units, &slope, &intercept);
      tmp = intercept + slope*default_value;
      
      ierr = verbPrintf(2, grid->com, 
			"  absent %-10s/ %-60s\n   %-16s\\ not found; using default constant %7.2f (%s)\n",
			short_name, long_name, "", tmp, glaciological_units_string);
      CHKERRQ(ierr);
      ierr = set(default_value); CHKERRQ(ierr);
    } else {			// otherwise leave it alone
      ierr = verbPrintf(2, grid->com, 
			"  absent %-10s/ %-60s\n   %-16s\\ not found; continuing without setting it\n",
			short_name, long_name, "");
      CHKERRQ(ierr);
    }
  } else {			// the variable was found successfully
    // Check if it is discrete
    if (use_interpolation_mask)
      nc.set_MaskInterp(&interpolation_mask);

    if (localp) {
      ierr = nc.regrid_local_var(varid, dims, lic, da, v,
				 use_interpolation_mask); CHKERRQ(ierr);
    } else {
      ierr = nc.regrid_global_var(varid, dims, lic, v,
				  use_interpolation_mask); CHKERRQ(ierr);
    }
    // We are done reading, and the data is in the units specified in the
    // bootstrapping file now. We need to check the range before changing the
    // units to avoid converting all the NetCDF attributes (i.e. valid_range
    // and others).
    ierr = check_range(nc.ncid, varid); CHKERRQ(ierr);
    // Now we need to get the units string from the file and convert the units,
    // because report_range expects the data to be in PISM (SI) units.
    
    bool input_has_units;
    utUnit input_units;

    ierr = nc.get_units(short_name, standard_name,
			input_has_units, input_units); CHKERRQ(ierr);

    if ( has_units && (!input_has_units) ) {
      ierr = verbPrintf(2, grid->com,
			"PISM WARNING: Variable '%s' ('%s') does not have the units attribute.\n"
			"              Assuming that it is in '%s'.\n",
			short_name, long_name, units_string); CHKERRQ(ierr);
      utCopy(&units, &input_units);
    }

    // Convert data:
    ierr = change_units(&input_units, &units); CHKERRQ(ierr);

    // Read the valid range info:
    ierr = read_valid_range(nc.ncid, varid); CHKERRQ(ierr);

    // We can report the success, and the range now:
    ierr = verbPrintf(2, grid->com, "  FOUND ");
    ierr = report_range(); CHKERRQ(ierr);
  } // end of if(exists)

  ierr = nc.close(); CHKERRQ(ierr);
  return 0;
}

//! Checks if an IceModelVec is allocated.  Terminates if not.
PetscErrorCode  IceModelVec::checkAllocated() {
  if (v == PETSC_NULL) {
    SETERRQ1(1,"IceModelVec ERROR: IceModelVec with short_name='%s' WAS NOT allocated\n",
             short_name);
  }
  return 0;
}

//! Checks if the access to the array is available.
PetscErrorCode  IceModelVec::checkHaveArray() {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  if (array == PETSC_NULL) {
    SETERRQ1(1,"array for IceModelVec with short_name='%s' not available\n"
               "  (REMEMBER TO RUN begin_access() before access and end_access() after access)\n",
               short_name);
  }
  return 0;
}


/*
PetscErrorCode  IceModelVec::checkSelfOwnsIt(const PetscInt i, const PetscInt j) {
  if (allocated == PETSC_FALSE) {
    SETERRQ3(1,"IceModelVec ERROR: (i,j)=(%d,%d) not in ownership range of processor %d\n",
             i,j,grid->rank);
  }
  return 0;
}
*/

//! Checks if an IceModelVec is allocated and calls DAVecGetArray.
PetscErrorCode  IceModelVec::begin_access() {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  if (array == NULL) {
    ierr = DAVecGetArray(da, v, &array); CHKERRQ(ierr);
  }
  return 0;
}

//! Checks if an IceModelVec is allocated and calls DAVecRestoreArray.
PetscErrorCode  IceModelVec::end_access() {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  if (array != NULL) {
    ierr = DAVecRestoreArray(da, v, &array); CHKERRQ(ierr);
    array = PETSC_NULL;
  }
  return 0;
}

//! Starts the communication of ghost points.
PetscErrorCode  IceModelVec::beginGhostComm() {
  PetscErrorCode ierr;
  if (!localp) {
    SETERRQ1(1,"makes no sense to communicate ghosts for GLOBAL IceModelVec! (has short_name='%s')\n",
               short_name);
  }
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(da, v, INSERT_VALUES, v);  CHKERRQ(ierr);
  return 0;
}

//! Ends the communication of ghost points.
PetscErrorCode  IceModelVec::endGhostComm() {
  PetscErrorCode ierr;
  if (!localp) {
    SETERRQ1(1,"makes no sense to communicate ghosts for GLOBAL IceModelVec! (has short_name='%s')\n",
               short_name);
  }
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(da, v, INSERT_VALUES, v); CHKERRQ(ierr);
  return 0;
}

//! Starts the communication of ghost points to IceModelVec destination.
PetscErrorCode  IceModelVec::beginGhostComm(IceModelVec &destination) {
  PetscErrorCode ierr;
  if (!localp) {
    SETERRQ1(1,"makes no sense to communicate ghosts for GLOBAL IceModelVec! (has short_name='%s')",
               short_name);
  }

  if (!destination.localp) {
    SETERRQ(1, "IceModelVec::beginGhostComm(): destination has to be local.");
  }

  if (dims != destination.dims) {
    SETERRQ(1, "IceModelVec::beginGhostComm(): operands have different numbers of dimensions.");
  }

  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = DALocalToLocalBegin(da, v, INSERT_VALUES, destination.v);  CHKERRQ(ierr);
  return 0;
}

//! Ends the communication of ghost points to IceModelVec destination.
PetscErrorCode  IceModelVec::endGhostComm(IceModelVec &destination) {
  PetscErrorCode ierr;
  if (!localp) {
    SETERRQ1(1,"makes no sense to communicate ghosts for GLOBAL IceModelVec! (has short_name='%s')\n",
               short_name);
  }

  if (!destination.localp) {
    SETERRQ(1, "IceModelVec::beginGhostComm(): destination has to be local.");
  }

  if (dims != destination.dims) {
    SETERRQ(1, "IceModelVec::beginGhostComm(): operands have different numbers of dimensions.");
  }

  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = DALocalToLocalEnd(da, v, INSERT_VALUES, destination.v); CHKERRQ(ierr);
  return 0;
}

//! Result: v[j] <- c for all j.
PetscErrorCode  IceModelVec::set(const PetscScalar c) {
  PetscErrorCode ierr;
  ierr = checkAllocated(); CHKERRQ(ierr);
  ierr = VecSet(v,c); CHKERRQ(ierr);
  return 0;
}

//! Checks the range of an IceModelVec read from a netCDF file.
/*! Uses valid_range, valid_min, valid_max, missing_value and _FillValue
    attributes, if present.

    This should be done after reading but before units conversion.
 */
PetscErrorCode IceModelVec::check_range(const int ncid, const int varid) {
  PetscErrorCode ierr;
  double bounds[2] = {0, 0};
  int use_min = 0, use_max = 0;
  NCTool nc(grid);
  nc.ncid = ncid;

  int stat;
  bool missing_value_exists = false,
    fillvalue_exists = false,
    valid_min_exists = false,
    valid_max_exists = false,
    valid_range_exists = false;
  double missing_value, fillvalue, valid_min, valid_max, valid_range[2];

  // Try to read all the attributes:

  stat = nc.get_att_double(varid, "missing_value", 1, &missing_value);
  if (stat == 0)
    missing_value_exists = true;

  stat = nc.get_att_double(varid, "_FillValue", 1, &fillvalue);
  if (stat == 0)
    fillvalue_exists = true;

  stat = nc.get_att_double(varid, "valid_min", 1, &valid_min);
  if (stat == 0)
    valid_min_exists = true;

  stat = nc.get_att_double(varid, "valid_max", 1, &valid_max);
  if (stat == 0)
    valid_max_exists = true;

  stat = nc.get_att_double(varid, "valid_range", 2, valid_range);
  if (stat == NC_NOERR)
    valid_range_exists = true;

  // Process attributes:
  if (valid_range_exists) {
    // consistency check
    if (valid_min_exists || valid_max_exists || missing_value_exists || fillvalue_exists) {
      ierr = verbPrintf(2, grid->com,
	 "PISM WARNING (variable '%s'): both valid_range and at least one of valid_min, valid_max\n"
	 "  missing_value and _FillValue are set; using valid_range\n",
	 short_name);
      CHKERRQ(ierr);
    }	// end of consistency check

    use_min = use_max = 1;
    bounds[0] = valid_range[0];
    bounds[1] = valid_range[1];
    // end of if(valid_range_exists)
  } else if (valid_min_exists && valid_max_exists) {
    // consistency check
    if (missing_value_exists || fillvalue_exists) {
      ierr = verbPrintf(2, grid->com,
	 "PISM WARNING (variable '%s'): both valid_min, valid_max and at least one of\n"
	 "  missing_value and _FillValue are set; using the valid_min, valid_max pair\n",
	 short_name);
         CHKERRQ(ierr);
    }	// end of consistency check

    use_min = use_max = 1;
    bounds[0] = valid_min;
    bounds[1] = valid_max;
  } else if (valid_min_exists) {
    // consistency check
    if (missing_value_exists || fillvalue_exists) {
      ierr = verbPrintf(2, grid->com,
	 "PISM WARNING (variable '%s'): both valid_min and at least one of\n"
	 "  missing_value and _FillValue are set; using valid_min\n",short_name);
         CHKERRQ(ierr);
    }	// end of consistency check

    use_min = 1; use_max = 0;
    bounds[0] = valid_min;
    // end of if (valid_min_exists)
  } else if (valid_max_exists) {
    // consistency check
    if (missing_value_exists || fillvalue_exists) {
      ierr = verbPrintf(2, grid->com,
	 "PISM WARNING (variable '%s'): both valid_max and at least one of\n"
	 "  missing_value and _FillValue are set; using valid_max\n",short_name);
         CHKERRQ(ierr);
    }	// end of consistency check

    use_min = 0; use_max = 1;
    bounds[1] = valid_max;
    // end of if (valid_max_exists)
  } else if (fillvalue_exists) {
    // consistency check
    if (missing_value_exists) {
      ierr = verbPrintf(2, grid->com,
	 "PISM WARNING (variable '%s'): both _FillValue and missing_value are set;\n"
	 "  using _FillValue\n",short_name);
	 CHKERRQ(ierr);
    }
      
    if (fillvalue < 0) {
      use_min = 1; use_max = 0;
      bounds[0] = fillvalue + 3e-6; // a trick to exclude this value
    } else {
      use_min = 0; use_max = 1;
      bounds[1] = fillvalue - 3e-6; // a trick to exclude this value
    }
    // end of if (fillvalue_exists)
  } else if (missing_value_exists) {
    ierr = verbPrintf(2, grid->com,
       "PISM WARNING (variable '%s'): missing_value attribute is deprecated by the NetCDF User's Guide;\n"
       "  ignoring it; please use valid_min, valid_max, valid_range or _FillValue attributes\n",
       short_name);
       CHKERRQ(ierr);
    use_min = use_max = 0;
    // end of if(missing_value_exists)
  }

  const double eps = 2e-16;
  double min, max;
  if (use_min || use_max) {
    ierr = range(min, max); CHKERRQ(ierr);
  }

  if (use_min && use_max) {
    if ((min < bounds[0] - eps) || (max > bounds[1] + eps)) {
      ierr = verbPrintf(2, grid->com,
	 "PISM WARNING: some values of '%s' are outside the valid range [%f, %f] %s\n",
	 short_name, bounds[0], bounds[1], units_string); CHKERRQ(ierr);
    }
  } else if (use_min) {
    if (min < bounds[0] - eps) {
      ierr = verbPrintf(2, grid->com,
	 "PISM WARNING: some values of '%s' are less than the valid minimum (%f)\n",
	 short_name, bounds[0]); CHKERRQ(ierr);
    }
  } else if (use_max) {
    if (max > bounds[1] + eps) {
      ierr = verbPrintf(2, grid->com,
	 "PISM WARNING: some values of '%s' are greater than the valid maximum (%f)\n",
	 short_name, bounds[0]); CHKERRQ(ierr);
    }
  }

  return 0;
}

//! Unit conversion.
/*! Convert the data stored in an IceModelVec from the units corresponding the
  the \from argument to ones corresponding to \c to.

  Note that this method *does not* change the units_string, so it has to be
  updated manually to reflect the change.

  Returns 0 on success, 1 if given units are incompatible, 2 on all the other
  UDUNITS-related errors.
*/
PetscErrorCode IceModelVec::change_units(utUnit *from, utUnit *to) {
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
			short_name, from_name, to_name);
      return 1;
    } else {			// some other error
      return 2;
    }
  }

  use_slope     = PetscAbsReal(slope - 1.0) > 1e-16;
  use_intercept = PetscAbsReal(intercept)   > 1e-16;

  if (use_slope && use_intercept) {
    ierr = scale(slope); CHKERRQ(ierr);
    ierr = shift(intercept); CHKERRQ(ierr);
  } else if (use_slope && !use_intercept) {
    ierr = scale(slope); CHKERRQ(ierr);
  } else if (!use_slope && use_intercept) {
    ierr = shift(intercept); CHKERRQ(ierr);
  }

  return 0;
}

//! Write an extra text attribute (not stored in IceModelVec).
/*! Call this after making sure that the variable exists in the output.

  \c tp has to point to a null-terminated string.
 */
PetscErrorCode IceModelVec::write_text_attr(const char filename[], const char name[], 
                                            const char *tp) {
  bool exists;
  int varid, ierr;
  NCTool nc(grid);

  ierr = nc.open_for_writing(filename, false); CHKERRQ(ierr); // do not replace the file

  ierr = nc.find_variable(short_name, standard_name, &varid, exists); CHKERRQ(ierr);
  if (!exists)
    SETERRQ(1, "Can't write attributes of an undefined variable.");

  if (grid->rank == 0) {
    ierr = nc_redef(nc.ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));

    ierr = nc_put_att_text(nc.ncid, varid, name, strlen(tp), tp);
       CHKERRQ(check_err(ierr,__LINE__,__FILE__));

    ierr = nc_enddef(nc.ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }

  return 0;
}

//! Write an extra scalar attribute (not stored in IceModelVec).
/*! Call this after making sure that the variable exists in the output 
 */
PetscErrorCode IceModelVec::write_scalar_attr(const char filename[], const char name[], 
                                              nc_type nctype, size_t len, const double *dp) {
  bool exists;
  int varid, ierr;
  NCTool nc(grid);

  ierr = nc.open_for_writing(filename, false); CHKERRQ(ierr); // do not replace the file

  ierr = nc.find_variable(short_name, standard_name, &varid, exists); CHKERRQ(ierr);
  if (!exists)
    SETERRQ(1, "Can't write attributes of an undefined variable.");

  if (grid->rank == 0) {
    ierr = nc_redef(nc.ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));

    ierr = nc_put_att_double(nc.ncid, varid, name, nctype, len, dp);
       CHKERRQ(check_err(ierr,__LINE__,__FILE__));

    ierr = nc_enddef(nc.ncid); CHKERRQ(check_err(ierr,__LINE__,__FILE__));
  }
  return 0;
}

//! Sets the valid_min NetCDF attribute.
PetscErrorCode IceModelVec::set_valid_min(PetscReal new_min) {
  valid_min = new_min;
  has_valid_min = true;
  return 0;
}

//! Sets the valid_max NetCDF attribute.
PetscErrorCode IceModelVec::set_valid_max(PetscReal new_max) {
  valid_max = new_max;
  has_valid_max = true;
  return 0;
}

//! Sets the valid_range NetCDF attribute.
PetscErrorCode IceModelVec::set_valid_range(PetscReal new_min, PetscReal new_max) {
  valid_min = new_min;
  has_valid_min = true;
  valid_max = new_max;
  has_valid_max = true;
  return 0;
}

//! Sets the coordinates NetCDF attribute.
/*! Coordinates are present by default; this method is here only to allow
    removing the attribute if it is not applicable.
 */
PetscErrorCode IceModelVec::set_coordinates(const char name[]) {
  if (name != NULL) {
    strncpy(coordinates, name, PETSC_MAX_PATH_LEN);
    has_coordinates = true;
  } else {
    has_coordinates = false;
  }

  return 0;
}

//! Reads the valid_range, valid_min valid_max attributes from an input file.
/*! Reads these attributes and and sets IceModelVec::valid_min and
  IceModelVec::valid_max if they were not set already.

  Note that if valid_range is present, then valid_min and valid_max are
  ignored.
 */
PetscErrorCode IceModelVec::read_valid_range(const int ncid, const int varid) {
  NCTool nc(grid);
  char  *input_units_string;
  utUnit input_units;
  double bounds[2], slope, intercept;
  int stat;

  // Never reset valid_min/max if any of them was set internally.
  if (has_valid_min || has_valid_max)
    return 0;

  nc.ncid = ncid;

  // Read the units: The following code ignores the units in the input file if
  // a) they are absent :-) b) they are invalid c) they are not compatible with
  // internal units.
  stat = nc.get_att_text(varid, "units", NULL, &input_units_string); CHKERRQ(stat);
  if (input_units_string != NULL) {
    stat = utScan(input_units_string, &input_units);
    if (stat != 0)
      utCopy(&units, &input_units);

    delete[] input_units_string;
  }

  if (utConvert(&input_units, &units, &slope, &intercept) != 0) {
    slope = 1;
    intercept = 0;
  }

  stat = nc.get_att_double(varid, "valid_range", 2, bounds);
  if (stat == 0) {		// valid_range is present
      has_valid_min = true;
      valid_min = intercept + slope*bounds[0];
      has_valid_max = true;
      valid_max = intercept + slope*bounds[1];
  } else {			// valid_range is absent or an error occured
    stat = nc.get_att_double(varid, "valid_min", 1, bounds);
    if (stat == 0) {		// valid_min is present
	has_valid_min = true;
	valid_min = intercept + slope*bounds[0];
    }

    stat = nc.get_att_double(varid, "valid_max", 1, bounds);
    if (stat == 0) {		// valid_max is present
	has_valid_max = true;
	valid_max = intercept + slope*bounds[0];
    }
  }

  return 0;
}

//! Checks if a value \c a in in the range of valid values of an IceModelVec.
bool IceModelVec::is_valid(PetscScalar a) {
  if (has_valid_min && has_valid_max)
    return (a >= valid_min) && (a <= valid_max);

  if (has_valid_min)
    return a >= valid_min;

  if (has_valid_max)       
    return a <= valid_max;

  return true;
}

/********* IceModelVec3 and IceModelVec3Bedrock: SEE SEPARATE FILE  iceModelVec3.cc    **********/

/********* IceModelVec2: SEE SEPARATE FILE  iceModelVec2.cc    **********/
