// Copyright (C) 2004-2010 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <sstream>
#include <cstring>
#include "iceModel.hh"
#include "pism_signal.h"
#include <petscvec.h>

//! Virtual.  Does nothing in \c IceModel.  Derived classes can do more computation in each time step.
PetscErrorCode IceModel::additionalAtStartTimestep() {
  return 0;
}


//! Virtual.  Does nothing in \c IceModel.  Derived classes can do more computation in each time step.
PetscErrorCode IceModel::additionalAtEndTimestep() {
  return 0;
}

//! Catch signals -USR1 and -TERM; in the former case save and continue; in the latter, save and stop.
/*!
Signal \c SIGTERM makes PISM end, saving state under original \c -o name 
(or default name).  We also add an indication to the history attribute 
of the output NetCDF file.

Signal \c SIGUSR1 makes PISM save state under a filename based on the
the name of the executable (e.g. \c pismr or \c pismv) and the current 
model year.  In addition the time series (\c -ts_file, etc.) is flushed out
There is no indication of these actions in the history attribute of the output (\c -o)
NetCDF file because there is no effect on it, but there is an indication at \c stdout.

Signal \c SIGUSR2 makes PISM flush time-series, without saving model state.
 */
int IceModel::endOfTimeStepHook() {
  
  if (pism_signal == SIGTERM) {
    verbPrintf(1, grid.com, 
       "\ncaught signal SIGTERM:  EXITING EARLY and saving with original filename.\n");
    char str[TEMPORARY_STRING_LENGTH];
    snprintf(str, sizeof(str), 
       "EARLY EXIT caused by signal SIGTERM.  Completed timestep at year=%.3f.",
       grid.year);
    stampHistory(str);
    return 1;
  }
  
  if (pism_signal == SIGUSR1) {
    char file_name[PETSC_MAX_PATH_LEN];
    snprintf(file_name, PETSC_MAX_PATH_LEN, "%s-%5.3f.nc",
             executable_short_name.c_str(), grid.year);
    verbPrintf(1, grid.com, 
       "\ncaught signal SIGUSR1:  Writing intermediate file `%s' and flushing time series.\n\n",
       file_name);
    pism_signal = 0;
    dumpToFile(file_name);

    // flush all the time-series buffers:
    vector<DiagnosticTimeseries*>::iterator i;
    for (i = timeseries.begin(); i < timeseries.end(); ++i) {
      (*i)->flush();
    }
  }

  if (pism_signal == SIGUSR2) {
    verbPrintf(1, grid.com, 
       "\ncaught signal SIGUSR2:  Flushing time series.\n\n");
    pism_signal = 0;

    // flush all the time-series buffers:
    vector<DiagnosticTimeseries*>::iterator i;
    for (i = timeseries.begin(); i < timeseries.end(); ++i) {
      (*i)->flush();
    }
  }

  return 0;
}


//! Build a history string from the command which invoked PISM.
PetscErrorCode  IceModel::stampHistoryCommand() {
  PetscErrorCode ierr;
  PetscInt argc;
  char **argv;
  
  ierr = PetscGetArgs(&argc, &argv); CHKERRQ(ierr);
  
  char startstr[TEMPORARY_STRING_LENGTH];

  snprintf(startstr, sizeof(startstr), 
           "PISM (%s) started on %d procs.", PISM_Revision, (int)grid.size);
  ierr = stampHistory(string(startstr)); CHKERRQ(ierr);

  // Create a string with space-separated command-line arguments:
  string cmdstr;
  for (int j = 0; j < argc; j++)
    cmdstr += string(" ") + argv[j];
  cmdstr += "\n";

  global_attributes.prepend_history(cmdstr);

  return 0;
}


//! Build the particular history string associated to the end of a PISM run.
PetscErrorCode  IceModel::stampHistoryEnd() {
  PetscErrorCode ierr;

  // timing stats
  PetscLogDouble current_time;
  PetscReal wall_clock_hours, proc_hours, mypph;
  ierr = PetscGetTime(&current_time); CHKERRQ(ierr);
  wall_clock_hours = (current_time - start_time) / 3600.0;
  proc_hours = grid.size * wall_clock_hours;
  mypph = (grid.year - grid.start_year) / proc_hours;

  // get PETSc's reported number of floating point ops (*not* per time) on this
  //   process, then sum over all processes
  PetscLogDouble flops, my_flops;
  MPI_Datatype mpi_type;
  ierr = PetscGetFlops(&my_flops); CHKERRQ(ierr);
  ierr = PetscDataTypeToMPIDataType(PETSC_DOUBLE, &mpi_type); CHKERRQ(ierr);
  MPI_Reduce(&my_flops, &flops, 1, mpi_type, MPI_SUM, 0, grid.com);

  // build and put string into global attribute "history"
  char str[TEMPORARY_STRING_LENGTH];
  snprintf(str, sizeof(str), 
    "PISM done.  Performance stats: %.4f wall clock hours, %.4f proc.-hours, %.4f model years per proc.-hour, PETSc MFlops = %.2f.",
    wall_clock_hours, proc_hours, mypph, flops * 1.0e-6);
  ierr = stampHistory(str); CHKERRQ(ierr);
  
  return 0;
}


//! Get time and user/host name and add it to the given string.
PetscErrorCode  IceModel::stampHistory(string str) {

  global_attributes.prepend_history(pism_username_prefix() + (str + "\n"));
  
  return 0;
}

//! Check if the thickness of the ice is too large and extend the grid if necessary.
/*!
  Extends the grid such that the new one has 2 (two) levels above the ice.
 */
PetscErrorCode IceModel::check_maximum_thickness() {
  PetscErrorCode  ierr;
  PetscReal H_min, H_max, dz_top;
  double *new_zlevels, *new_zblevels;
  const int old_Mz = grid.Mz;
  int N = 0; 			// the number of new levels

  ierr = vH.range(H_min, H_max); CHKERRQ(ierr);
  if (grid.Lz >= H_max) return 0;

  if (grid.initial_Mz == 0)
    grid.initial_Mz = grid.Mz;
  else if (grid.Mz > grid.initial_Mz * 2) {
    ierr = PetscPrintf(grid.com,
		       "\n"
		       "PISM ERROR: Max ice thickness (%7.4f m) is greater than the height of the computational box (%7.4f m)"
		       " AND the grid has twice the initial number of vertical levels (%d) already. Exiting...\n",
		       H_max, grid.Lz, grid.initial_Mz); CHKERRQ(ierr);
    PetscEnd();
  }

  // So, we need to extend the grid. We find dz at the top of the grid,
  // create new zlevels and zblevels, then extend all the IceModelVec3s.

  // Find the vertical grid spacing at the top of the grid:
  dz_top = grid.Lz - grid.zlevels[old_Mz - 2];

  // Find the number of new levels:
  while (grid.Lz + N * dz_top <= H_max) N++;
  // This makes sure that we always have *exactly* two levels strictly above the ice:
  if (grid.Lz + N * dz_top > H_max) N += 1;
  else N += 2;

  ierr = verbPrintf(2, grid.com,
		    "\n"
		    "PISM WARNING: max ice thickness (%7.4f m) is greater than the height of the computational box (%7.4f m)...\n"
		    "              Adding %d new grid levels %7.4f m apart...\n",
		    H_max, grid.Lz, N, dz_top); CHKERRQ(ierr);

  // Create new zlevels and zblevels:
  new_zblevels = new double[grid.Mbz];
  new_zlevels = new double[old_Mz + N];

  for (int j = 0; j < grid.Mbz; j++)
    new_zblevels[j] = grid.zblevels[j];
  for (int j = 0; j < old_Mz; j++)
    new_zlevels[j] = grid.zlevels[j];

  // Fill the new levels:
  for (int j = 0; j < N; j++)
    new_zlevels[old_Mz + j] = grid.Lz + dz_top * (j + 1);

  ierr = grid.set_vertical_levels(old_Mz + N, grid.Mbz,
				  new_zlevels, new_zblevels); CHKERRQ(ierr);
  delete[] new_zlevels;
  delete[] new_zblevels;

  // Done with the grid. Now we need to extend IceModelVec3s.

  // We use surface temperatures to extend T3 and Tnew3. We get them from the
  // PISMSurfaceModel.

  if (surface != PETSC_NULL) {
    ierr = surface->ice_surface_temperature(grid.year, 0.0, artm); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: surface == PETSC_NULL");
  }

  // for extending the variables Enth3 and Enthnew3 vertically, put into
  //   vWork2d[0] the enthalpy of the air
  ierr = vWork2d[0].begin_access(); CHKERRQ(ierr);
  ierr = artm.begin_access(); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      ierr = EC->getEnthPermissive(
         artm(i,j),0.0,EC->getPressureFromDepth(0.0),vWork2d[0](i,j));
         CHKERRQ(ierr);
    }
  }
  ierr = vWork2d[0].end_access(); CHKERRQ(ierr);
  ierr = artm.end_access(); CHKERRQ(ierr);

  // Model state 3D vectors:
  ierr =     u3.extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  ierr =     v3.extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  ierr =     w3.extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  ierr = Sigma3.extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  ierr =  Enth3.extend_vertically(old_Mz, vWork2d[0]); CHKERRQ(ierr);

  // Work 3D vectors:
  ierr =      Enthnew3.extend_vertically(old_Mz, vWork2d[0]); CHKERRQ(ierr);
  ierr = Sigmastag3[0].extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  ierr = Sigmastag3[1].extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  ierr =     Istag3[0].extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  ierr =     Istag3[1].extend_vertically(old_Mz, 0); CHKERRQ(ierr);

  if (config.get_flag("do_cold_ice_methods")) {
    ierr =    T3.extend_vertically(old_Mz, artm); CHKERRQ(ierr);
    ierr = Tnew3.extend_vertically(old_Mz, artm); CHKERRQ(ierr);
  }

  // deal with 3D age conditionally
  if (config.get_flag("do_age")) {
    ierr = tau3.extend_vertically(old_Mz, 0); CHKERRQ(ierr);
    ierr = taunew3.extend_vertically(old_Mz, 0); CHKERRQ(ierr);
  }
  
  ierr = check_maximum_thickness_hook(old_Mz); CHKERRQ(ierr);

  if (save_snapshots && (!split_snapshots) ) {
    char tmp[20];
    snprintf(tmp, 20, "%d", grid.Mz);
    
    snapshots_filename = pism_filename_add_suffix(snapshots_filename, "-Mz", tmp);
    snapshots_file_is_ready = false;

    ierr = verbPrintf(2, grid.com,
		      "NOTE: Further snapshots will be saved to '%s'...\n",
		      snapshots_filename.c_str()); CHKERRQ(ierr);
  }

  return 0;
}


//! Allows derived classes to extend their own IceModelVec3's in vertical.
/*! Base class version does absolutely nothing. */
PetscErrorCode IceModel::check_maximum_thickness_hook(const int /*old_Mz*/) {
  return 0;
}


PetscErrorCode IceModel::report_grid_parameters() {
  PetscErrorCode ierr;

  ierr = verbPrintf(2,grid.com, "computational domain and grid:\n"); CHKERRQ(ierr);
  // report on computational box
  ierr = verbPrintf(2,grid.com, 
           "           spatial domain   %.2f km x %.2f km",
           2*grid.Lx/1000.0,2*grid.Ly/1000.0); CHKERRQ(ierr);
  if (grid.Mbz > 1) {
    ierr = verbPrintf(2,grid.com," x (%.2f m + %.2f m bedrock)\n"
         ,grid.Lz,grid.Lbz); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com," x %.2f m\n",grid.Lz); CHKERRQ(ierr);
  }

  ierr = verbPrintf(2, grid.com,
           "            time interval   [ %.2f a, %.2f a ]; run length = %.4f a\n",
		    grid.start_year, grid.end_year, grid.end_year - grid.start_year);
  
  // report on grid cell dims
  ierr = verbPrintf(2,grid.com, 
           "     horizontal grid cell   %.2f km x %.2f km\n",
                    grid.dx/1000.0,grid.dy/1000.0); CHKERRQ(ierr);
  if (grid.ice_vertical_spacing == EQUAL) {
    ierr = verbPrintf(2,grid.com, 
           "  vertical spacing in ice   dz = %.3f m (equal spacing)\n",
                    grid.dzMIN); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2,grid.com, 
           "  vertical spacing in ice   uneven, %d levels, %.3f m < dz < %.3f m\n",
		    grid.Mz, grid.dzMIN, grid.dzMAX); CHKERRQ(ierr);
  }

  if (grid.Mbz > 1) {
    if (grid.bed_vertical_spacing == EQUAL) {
      ierr = verbPrintf(2,grid.com, 
           "  vert spacing in bedrock   dz = %.3f m (equal spacing)\n",
			grid.zblevels[1]-grid.zblevels[0]); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2,grid.com, 
			"  vert spacing in bedrock   uneven, %d levels, %.3f m < dz < %.3f m\n",
			grid.Mbz, grid.dzbMIN, grid.dzbMAX); CHKERRQ(ierr);
    }
    ierr = verbPrintf(3,grid.com, 
           "  fine spacing in conservation of energy and age:\n"
           "                            fMz = %d, fdz = %.3f m, fMbz = %d m\n",
           grid.Mz_fine, grid.dz_fine, grid.Mbz_fine); CHKERRQ(ierr);
  } else { // no bedrock case
    ierr = verbPrintf(3,grid.com, 
           "   fine spacing used in energy/age   fMz = %d, fdz = %.3f m\n",
           grid.Mz_fine, grid.dz_fine); CHKERRQ(ierr);
  }
  if (grid.Mz_fine > 1000) {
    ierr = verbPrintf(2,grid.com,
      "\n\nWARNING: Using more than 1000 ice vertical levels internally in energy/age computation!\n\n");
      CHKERRQ(ierr);
  }

  // if -verbose (=-verbose 3) then (somewhat redundantly) list parameters of grid
  ierr = grid.printInfo(3); CHKERRQ(ierr);

  // if -verbose 5 then more stuff
  ierr = verbPrintf(5,grid.com,
       "  REALLY verbose output on IceGrid:\n"); CHKERRQ(ierr);
  ierr = grid.printVertLevels(5); CHKERRQ(ierr);
  
  return 0;
}


bool IceModel::issounding(const PetscInt i, const PetscInt j){ 
  return ((i == id) && (j == jd));
}

void IceModel::attach_surface_model(PISMSurfaceModel *my_surface) {
  surface = my_surface;
}

void IceModel::attach_ocean_model(PISMOceanModel *my_ocean) {
  ocean = my_ocean;
}

//! Computes the area of a triangle using vector cross product.
static PetscReal triangle_area(PetscReal *A, PetscReal *B, PetscReal *C) {
  PetscReal V1[3], V2[3];
  for (int j = 0; j < 3; ++j) {
    V1[j] = B[j] - A[j];
    V2[j] = C[j] - A[j];
  }

  return 0.5*sqrt(PetscSqr(V1[1]*V2[2] - V2[1]*V1[2]) +
		  PetscSqr(V1[0]*V2[2] - V2[0]*V1[2]) +
		  PetscSqr(V1[0]*V2[1] - V2[0]*V1[1]));
}

static inline PetscReal sin_degrees(PetscReal deg) {
  return sin(deg * M_PI / 180.0);
}

static inline PetscReal cos_degrees(PetscReal deg) {
  return cos(deg * M_PI / 180.0);
}

//! Computes geocentric x-coordinate of a point using reference ellipsoid parameters.
//! (see http://en.wikipedia.org/wiki/Reference_ellipsoid)
static PetscReal geo_x(PetscReal a, PetscReal b,
		       PetscReal lon, PetscReal lat) {
  const PetscReal oe = acos(b/a),
    N = a/sqrt( 1 - PetscSqr( sin(oe)*sin_degrees(lat) ) );

  return N * cos_degrees(lat) * cos_degrees(lon);
}

//! Computes geocentric y-coordinate of a point using reference ellipsoid parameters.
//! (see http://en.wikipedia.org/wiki/Reference_ellipsoid)
static PetscReal geo_y(PetscReal a, PetscReal b,
		       PetscReal lon, PetscReal lat) {
  const PetscReal oe = acos(b/a),
    N = a/sqrt( 1 - PetscSqr( sin(oe)*sin_degrees(lat) ) );

  return N * cos_degrees(lat) * sin_degrees(lon);
}

//! Computes geocentric z-coordinate of a point using reference ellipsoid parameters.
//! (see http://en.wikipedia.org/wiki/Reference_ellipsoid)
static PetscReal geo_z(PetscReal a, PetscReal b,
		       PetscReal lon, PetscReal lat) {
  const PetscReal oe = acos(b/a),
    N = a/sqrt( 1 - PetscSqr( sin(oe)*sin_degrees(lat) ) );

  return N * sin_degrees(lat) * PetscSqr(b/a);
}

/*!
    Allocate and compute corrected cell areas. Uses linear interpolation to
    find latitudes and longitudes of grid corners, WGS84 parameters to compute
    cartesian coordinates of grid corners and vector products to compute areas
    of resulting triangles.
   
    \note Latitude and longitude fields are \b not periodic, so computing
    corrected areas for cells at the grid boundary is not feasible. This should
    not matter, since these cells should be (and in most cases are) ice-free.
    
    \note Using linear interpolation introduces errors in lon/lat coordinates
    of cell corners, but the corresponding ice volume error (present day
    Greenland, 5km grid) is about 3.11e-04 %.
*/
PetscErrorCode IceModel::correct_cell_areas() {
  PetscErrorCode ierr;

  if (!mapping.has("grid_mapping_name"))
    return 0;

  if (!config.get_flag("correct_cell_areas"))
    return 0;

  ierr = verbPrintf(2,grid.com,
		    "* Computing corrected cell areas using WGS84 datum...\n"); CHKERRQ(ierr);

  // allocate cell_area
  ierr = cell_area.create(grid, "cell_area", false); CHKERRQ(ierr);
  ierr = cell_area.set_attrs("diagnostic", "corrected cell areas", "m2", ""); CHKERRQ(ierr);
  cell_area.time_independent = true;
  ierr = cell_area.set_glaciological_units("km2"); CHKERRQ(ierr);
  cell_area.write_in_glaciological_units = true;
  ierr = variables.add(cell_area); CHKERRQ(ierr);

  PetscReal a = config.get("WGS84_semimajor_axis"),
    b = config.get("WGS84_semiminor_axis");

  // value to use near the grid boundary
  ierr = cell_area.set(grid.dx * grid.dy);

// Cell layout:
// (nw)        (ne)
// +-----------+
// |           |
// |           |
// |     o     |
// |           |
// |           |
// +-----------+
// (sw)        (se)

  PetscScalar **lat, **lon;
  ierr =  cell_area.begin_access(); CHKERRQ(ierr);
  ierr =  vLatitude.get_array(lat); CHKERRQ(ierr);
  ierr = vLongitude.get_array(lon); CHKERRQ(ierr);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    if ((i == 0) || (i == grid.Mx-1))
      continue;

    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if ((j == 0) || (j == grid.My-1))
	continue;

      PetscReal 
	// latitudes
        lat_nw = 0.5 * ( lat[i][j] + lat[i-1][j+1] ),
        lat_sw = 0.5 * ( lat[i][j] + lat[i-1][j-1] ),
        lat_se = 0.5 * ( lat[i][j] + lat[i+1][j-1] ),
        lat_ne = 0.5 * ( lat[i][j] + lat[i+1][j+1] ),
	// longitudes
        lon_nw = 0.5 * ( lon[i][j] + lon[i-1][j+1] ),
        lon_sw = 0.5 * ( lon[i][j] + lon[i-1][j-1] ),
        lon_se = 0.5 * ( lon[i][j] + lon[i+1][j-1] ),
        lon_ne = 0.5 * ( lon[i][j] + lon[i+1][j+1] );

      // geocentric coordinates:
      PetscReal sw[3] = {geo_x(a, b, lon_sw, lat_sw),
			 geo_y(a, b, lon_sw, lat_sw),
			 geo_z(a, b, lon_sw, lat_sw)};

      PetscReal se[3] = {geo_x(a, b, lon_se, lat_se),
			 geo_y(a, b, lon_se, lat_se),
			 geo_z(a, b, lon_se, lat_se)};

      PetscReal ne[3] = {geo_x(a, b, lon_ne, lat_ne),
			 geo_y(a, b, lon_ne, lat_ne),
			 geo_z(a, b, lon_ne, lat_ne)};

      PetscReal nw[3] = {geo_x(a, b, lon_nw, lat_nw),
			 geo_y(a, b, lon_nw, lat_nw),
			 geo_z(a, b, lon_nw, lat_nw)};

      cell_area(i,j) = triangle_area(sw, se, ne) + triangle_area(ne, nw, sw);
    }
  }
  ierr =  cell_area.end_access(); CHKERRQ(ierr);
  ierr =  vLatitude.end_access(); CHKERRQ(ierr);
  ierr = vLongitude.end_access(); CHKERRQ(ierr);

  return 0;
}

