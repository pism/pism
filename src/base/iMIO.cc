// Copyright (C) 2004-2009 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <cstdio>
#include <petscda.h>
#include "iceModel.hh"

bool IceModel::hasSuffix(const char* fname, const char *suffix) const {
  int flen = strlen(fname);
  int slen = strlen(suffix);
  if (strcmp(fname + flen - slen, suffix) == 0) {
    return true;
  } else {
    return false;
  }
}

PetscErrorCode  IceModel::setStartRunEndYearsFromOptions(const PetscTruth grid_p_year_VALID) {
  PetscErrorCode ierr;

  // read options about year of start, year of end, number of run years;
  // note grid.year has already been set from input file
  PetscScalar usrStartYear, usrEndYear, usrRunYears;
  PetscTruth ysSet, yeSet, ySet;
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-ys", &usrStartYear, &ysSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-ye", &usrEndYear, &yeSet); CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL, "-y", &usrRunYears, &ySet); CHKERRQ(ierr);
  if (ysSet == PETSC_TRUE) {
    // user option overwrites data
    ierr = setStartYear(usrStartYear); CHKERRQ(ierr);
    grid.year = usrStartYear;
  } else if (grid_p_year_VALID == PETSC_TRUE) {
    ierr = setStartYear(grid.year); CHKERRQ(ierr);
  } // else do nothing; defaults are set
  if (yeSet == PETSC_TRUE) {
    if (usrEndYear < startYear) {
      ierr = PetscPrintf(grid.com,
			"PISM ERROR: -ye (%3.3f) is less than -ys (%3.3f) (or input file year or default).\n"
			"PISM cannot run backward in time.\n",
			usrEndYear, startYear); CHKERRQ(ierr);
      PetscEnd();
    }
    if (ySet == PETSC_TRUE) {
      ierr = verbPrintf(1,grid.com,"WARNING: -y option ignored.  -ye used instead.\n"); CHKERRQ(ierr);
    }
    ierr = setEndYear(usrEndYear); CHKERRQ(ierr);
  } else if (ySet == PETSC_TRUE) {
    ierr = setEndYear(usrRunYears + startYear); CHKERRQ(ierr);
  } else {
    ierr = setEndYear(run_year_default + startYear); CHKERRQ(ierr);
  }
  
  yearsStartRunEndDetermined = PETSC_TRUE;
  return 0;
}
  
//! Save model state in NetCDF format (and save variables in Matlab format if desired).
/*! 
Optionally allows saving of full velocity field.

Calls dumpToFile() and writeMatlabVars() to do the actual work.
 */
PetscErrorCode  IceModel::writeFiles(const char* default_filename,
                                     const PetscTruth forceFullDiagnostics) {
  PetscErrorCode ierr;
  char filename[PETSC_MAX_PATH_LEN];

  if (doPDD == PETSC_TRUE) { // want to save snow accumulation map, not net accumulation
    ierr = putBackSnowAccumPDD(); CHKERRQ(ierr);
  }
  
  ierr = stampHistoryEnd(); CHKERRQ(ierr);

  PetscTruth o_set;
  ierr = PetscOptionsGetString(PETSC_NULL, "-o", filename, PETSC_MAX_PATH_LEN, &o_set); CHKERRQ(ierr);

  // Use the default if the output file name was not given:
  if (!o_set)
    strncpy(filename, default_filename, PETSC_MAX_PATH_LEN);

  if (!hasSuffix(filename, ".nc")) {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: output file name does not have the '.nc' suffix!\n");
    CHKERRQ(ierr);
  }

  PetscTruth userWantsFull;
  ierr = PetscOptionsHasName(PETSC_NULL, "-f3d", &userWantsFull); CHKERRQ(ierr);

  if ((forceFullDiagnostics == PETSC_TRUE) || (userWantsFull == PETSC_TRUE)) {
    ierr = verbPrintf(2, grid.com, 
		      "Writing model state, with full 3D velocities, to file `%s'",
		      filename); CHKERRQ(ierr);

    ierr = dumpToFile(filename); CHKERRQ(ierr);
    ierr = write3DPlusToFile(filename); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, "Writing model state to file `%s'",
		      filename); CHKERRQ(ierr);
    ierr = dumpToFile(filename); CHKERRQ(ierr);
  }

  // write out individual variables out to Matlab file
  char       matf[PETSC_MAX_PATH_LEN];
  PetscTruth matoSet, matvSet;
  ierr = PetscOptionsGetString(PETSC_NULL, "-mato", matf, PETSC_MAX_PATH_LEN, &matoSet); 
           CHKERRQ(ierr);
  if (matoSet == PETSC_FALSE) {// put default name in matf; perhaps user set "-matv" only
    strcpy(matf, "pism_views");
  }
  strcpy(matlabOutVars, "\0");
  ierr = PetscOptionsGetString(PETSC_NULL, "-matv", matlabOutVars, PETSC_MAX_PATH_LEN, &matvSet); 
            CHKERRQ(ierr);
  if (matvSet == PETSC_TRUE) {
    strcat(matf, ".m");
    ierr = verbPrintf(1, grid.com, 
       "\n ... writing variables %s to Matlab file `%s'", matlabOutVars, matf); CHKERRQ(ierr);
    ierr = writeMatlabVars(matf); CHKERRQ(ierr); // see iMmatlab.cc
  }

  return 0;
}

PetscErrorCode IceModel::dumpToFile(const char *filename) {
  PetscErrorCode ierr;
  PetscTruth append = PETSC_FALSE;
  NCTool nc(&grid);

  ierr = PetscOptionsHasName(PETSC_NULL, "-a", &append); CHKERRQ(ierr);
  if (append) {
    ierr = verbPrintf(2, grid.com, "\nWill append to '%s' if possible.\n", filename); CHKERRQ(ierr);
  }

  // Prepare the file
  ierr = nc.open_for_writing(filename, append == PETSC_FALSE); CHKERRQ(ierr);
  ierr = nc.append_time(grid.year * secpera); CHKERRQ(ierr);
  ierr = nc.write_history(history); CHKERRQ(ierr); // append the history
  ierr = nc.write_polar_stereographic(psParams.svlfp, psParams.lopo, psParams.sp); CHKERRQ(ierr);
  ierr = nc.write_global_attrs(useSSAVelocity, "CF-1.0"); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = write_model_state(filename);  CHKERRQ(ierr);
  ierr = write_extra_fields(filename); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::write_model_state(const char filename[]) {
  PetscErrorCode ierr;

  // 2D model quantities
  ierr =         vH.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  ierr = vLongitude.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  ierr =  vLatitude.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  ierr =      vMask.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  ierr =     vHmelt.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  ierr =       vbed.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  ierr =    vuplift.write(filename, NC_DOUBLE); CHKERRQ(ierr);

  if (useSSAVelocity) {
    ierr = vubarSSA.write(filename, NC_DOUBLE); CHKERRQ(ierr);
    ierr = vvbarSSA.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  }

  // 3D model quantities
  ierr =   T3.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  ierr =  Tb3.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  ierr = tau3.write(filename, NC_DOUBLE); CHKERRQ(ierr);

  // 2D climate quantities
  ierr =    vTs.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  ierr =   vGhf.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  ierr = vAccum.write(filename, NC_DOUBLE); CHKERRQ(ierr);

  // write tillphi = till friction angle in degrees
  ierr = vtillphi.write(filename, NC_DOUBLE); CHKERRQ(ierr);

  // 2D diagnostic quantities
  // note h is diagnostic because it is recomputed by h=H+b at each time step
  // these are not written in MKS units because they are intended to be viewed,
  // not read by programs; IS THIS THE RIGHT CHOICE?
  ierr = vh.write(filename, NC_FLOAT); CHKERRQ(ierr);

  ierr = vdHdt.write(filename, NC_FLOAT); CHKERRQ(ierr);

  // Create the mask of zeros and ones
  // 1 - ice thickness > 0
  // 0 - ice-free location
  PetscScalar **M, **H;
  ierr = vWork2d[5].get_array(M);
  ierr = vH.get_array(H);
  for (PetscInt i=grid.xs; i<grid.xs+grid.xm; ++i) {
    for (PetscInt j=grid.ys; j<grid.ys+grid.ym; ++j) {
      if (H[i][j] > 0.0) {
	M[i][j] = 1.0;
      } else {
	M[i][j] = 0.0;
      }
    }
  }
  ierr = vWork2d[5].end_access(); CHKERRQ(ierr);
  ierr = vH.end_access(); CHKERRQ(ierr);

  // compute cbar = sqrt(ubar^2 + vbar^2) and save it
  ierr = getMagnitudeOf2dVectorField(vubar, vvbar, vWork2d[0]); CHKERRQ(ierr);
  ierr = vWork2d[0].multiply_by(vWork2d[5]); CHKERRQ(ierr); // mask out ice-free areas

  ierr = vWork2d[0].set_name("cbar"); CHKERRQ(ierr);
  ierr = vWork2d[0].set_attrs("diagnostic", 
            "magnitude of vertically-integrated horizontal velocity of ice",
	    "m s-1", NULL); CHKERRQ(ierr);
  ierr = vWork2d[0].set_glaciological_units("m year-1", secpera); CHKERRQ(ierr);
  ierr = vWork2d[0].write(filename, NC_FLOAT); CHKERRQ(ierr);

  // compute cflx = cbar .* thk and save it
  ierr = vWork2d[0].multiply_by(vH, vWork2d[1]); CHKERRQ(ierr);
  ierr = vWork2d[1].set_name("cflx"); CHKERRQ(ierr);
  ierr = vWork2d[1].set_attrs("diagnostic", 
             "magnitude of vertically-integrated horizontal flux of ice",
	     "m2 s-1", NULL); CHKERRQ(ierr);
  ierr = vWork2d[1].set_glaciological_units("m2 year-1", secpera); CHKERRQ(ierr);
  ierr = vWork2d[1].write(filename, NC_FLOAT); CHKERRQ(ierr);

  // compute cbase  = sqrt(u|_{z=0}^2 + v|_{z=0}^2) and save it
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = u3.getHorSlice(vWork2d[0], 0.0); CHKERRQ(ierr); // vWork2d[0] = u_{z=0}
  ierr = v3.getHorSlice(vWork2d[1], 0.0); CHKERRQ(ierr); // vWork2d[1] = v_{z=0}
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);

  ierr = getMagnitudeOf2dVectorField(vWork2d[0],vWork2d[1],vWork2d[2]); CHKERRQ(ierr);
  ierr = vWork2d[2].multiply_by(vWork2d[5]); CHKERRQ(ierr); // mask out ice-free areas

  ierr = vWork2d[2].set_name("cbase"); CHKERRQ(ierr);
  ierr = vWork2d[2].set_attrs("diagnostic", 
             "magnitude of horizontal velocity of ice at base of ice",
	     "m s-1", NULL); CHKERRQ(ierr);
  ierr = vWork2d[2].set_glaciological_units("m year-1", secpera); CHKERRQ(ierr);
  ierr = vWork2d[2].write(filename, NC_FLOAT); CHKERRQ(ierr);

  // compute csurf = sqrt(u|_surface^2 + v|_surface^2) and save it
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = u3.getSurfaceValues(vWork2d[0], vH); CHKERRQ(ierr);
  ierr = v3.getSurfaceValues(vWork2d[1], vH); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);

  ierr = getMagnitudeOf2dVectorField(vWork2d[0],vWork2d[1],vWork2d[0]); CHKERRQ(ierr);
  ierr = vWork2d[0].multiply_by(vWork2d[5]); CHKERRQ(ierr); // mask out ice-free areas

  ierr = vWork2d[0].set_name("csurf"); CHKERRQ(ierr);
  ierr = vWork2d[0].set_attrs("diagnostic", 
             "magnitude of horizontal velocity of ice at ice surface",
	     "m s-1", NULL); CHKERRQ(ierr);
  ierr = vWork2d[0].set_glaciological_units("m year-1", secpera); CHKERRQ(ierr);
  ierr = vWork2d[0].write(filename, NC_FLOAT); CHKERRQ(ierr);

  // compute wvelsurf, the surface values of vertical velocity
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = w3.getSurfaceValues(vWork2d[0], vH); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);

  ierr = vWork2d[0].set_name("wvelsurf"); CHKERRQ(ierr);
  ierr = vWork2d[0].set_attrs("diagnostic", "vertical velocity of ice at ice surface",
			      "m s-1", NULL); CHKERRQ(ierr);
  ierr = vWork2d[0].set_glaciological_units("m year-1", secpera); CHKERRQ(ierr);
  ierr = vWork2d[0].write(filename, NC_FLOAT); CHKERRQ(ierr);

  // compute magnitude of basal shear stress = rho g H |grad h|
  ierr = computeDrivingStress(vWork2d[0],vWork2d[1]); CHKERRQ(ierr);
  ierr = getMagnitudeOf2dVectorField(vWork2d[0],vWork2d[1],vWork2d[2]); CHKERRQ(ierr);

  ierr = vWork2d[2].set_name("taud"); CHKERRQ(ierr);
  ierr = vWork2d[2].set_attrs("diagnostic",
             "magnitude of driving shear stress at base of ice",
	     "Pa", NULL); CHKERRQ(ierr);
  ierr = vWork2d[2].set_glaciological_units("Pa", 1.0); CHKERRQ(ierr);
  ierr = vWork2d[2].write(filename, NC_FLOAT); CHKERRQ(ierr);

  // write out yield stress
  ierr = vtauc.write(filename, NC_FLOAT); CHKERRQ(ierr);

  return 0;
}

// Writes extra fields to the output file \c filename. Does nothing in the base
// class.
PetscErrorCode IceModel::write_extra_fields(const char filename[]) {
  // Do nothing.
  return 0;
}

PetscErrorCode IceModel::write3DPlusToFile(const char filename[]) {
  PetscErrorCode ierr;

  // 3D velocity fields
  ierr = u3.write(filename, NC_FLOAT); CHKERRQ(ierr);
  ierr = v3.write(filename, NC_FLOAT); CHKERRQ(ierr);
  ierr = w3.write(filename, NC_FLOAT); CHKERRQ(ierr);

  // horizontal components of surface values of velocity; note wvelsurf already written
  ierr = u3.begin_access(); CHKERRQ(ierr);
  ierr = u3.getSurfaceValues(vWork2d[0], vH); CHKERRQ(ierr);
  ierr = u3.end_access(); CHKERRQ(ierr);

  ierr = vWork2d[0].set_name("uvelsurf"); CHKERRQ(ierr);
  ierr = vWork2d[0].set_attrs(
              "diagnostic", "x component of velocity of ice at ice surface",
	      "m s-1", NULL); CHKERRQ(ierr);
  ierr = vWork2d[0].set_glaciological_units("m year-1", secpera); CHKERRQ(ierr);
  ierr = vWork2d[0].write(filename, NC_FLOAT); CHKERRQ(ierr);

  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = v3.getSurfaceValues(vWork2d[0], vH); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);

  ierr = vWork2d[0].set_name("vvelsurf"); CHKERRQ(ierr);
  ierr = vWork2d[0].set_attrs(
              "diagnostic", "y component of velocity of ice at ice surface",
	      "m s-1", NULL); CHKERRQ(ierr);
  ierr = vWork2d[0].set_glaciological_units("m year-1", secpera); CHKERRQ(ierr);
  ierr = vWork2d[0].write(filename, NC_FLOAT); CHKERRQ(ierr);

  return 0;
}


//! When reading a saved PISM model state, warn the user if options <tt>-Mx,-My,-Mz,-Mbz</tt> have been ignored.
PetscErrorCode IceModel::warnUserOptionsIgnored(const char *fname) {
  PetscErrorCode ierr;
  PetscInt       ignor;
  PetscTruth     M_Set;

  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mx", &ignor, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,
             "WARNING: user option -Mx ignored; value read from file %s\n", fname);
             CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetInt(PETSC_NULL, "-My", &ignor, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,
             "WARNING: user option -My ignored; value read from file %s\n", fname);
             CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mz", &ignor, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,
             "WARNING: user option -Mz ignored; value read from file %s\n", fname);
             CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mbz", &ignor, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,
              "WARNING: user option -Mbz ignored; value read from file %s\n", fname);
              CHKERRQ(ierr);
  }
  return 0;
}



//! Read a saved PISM model state in NetCDF format, for complete initialization of an evolution or diagnostic run.
/*! 
When initializing from a NetCDF input file, the input file determines 
the number of grid points (Mx,My,Mz,Mbz) and the dimensions (Lx,Ly,Lz) of the computational box.   
The user is warned when their command line options "-Mx", "-My", "-Mz", "-Mbz" are overridden.  
 */
PetscErrorCode IceModel::initFromFile(const char *fname) {
  PetscErrorCode  ierr;
  int         stat;
  NCTool nc(&grid);

  ierr = verbPrintf(2,grid.com,"initializing from NetCDF file '%s'...\n",
                     fname); CHKERRQ(ierr);

  bool file_exists = false;
  ierr = nc.open_for_reading(fname, file_exists); CHKERRQ(ierr);
  if (!file_exists) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: Can't open file '%s'.\n", fname); CHKERRQ(ierr);
    PetscEnd();
  }

  // note user option setting of -Lx,-Ly,-Lz will overwrite the corresponding settings from 
  // this file but that the file's settings of Mx,My,Mz,Mbz will overwrite the user options 
  grid_info g;
  ierr = nc.get_grid_info(g);
  grid.year = g.time / secpera;
  grid.Mx = g.x_len;
  grid.My = g.y_len;
  grid.Mz = g.z_len;
  grid.Mbz = g.zb_len;
  // grid.Lx, grid.Ly set from g.x_max, g.y_max below in call to grid.rescale_using_zlevels()

  double *zlevs, *zblevs;
  zlevs = new double[grid.Mz];
  zblevs = new double[grid.Mbz];
  ierr = nc.get_vertical_dims(grid.Mz, grid.Mbz, zlevs, zblevs); CHKERRQ(ierr);

  // re-allocate and fill grid.zlevels & zblevels with read values
  delete [] grid.zlevels;
  delete [] grid.zblevels;
  grid.zlevels = new PetscScalar[grid.Mz];
  grid.zblevels = new PetscScalar[grid.Mbz];
  for (PetscInt k = 0; k < grid.Mz; k++) {
    grid.zlevels[k] = (PetscScalar) zlevs[k];
  }
  for (PetscInt k = 0; k < grid.Mbz; k++) {
    grid.zblevels[k] = (PetscScalar) zblevs[k];
  }
  delete [] zlevs;  delete [] zblevs;  // done with these
  
  // set DA and create vecs; we have enough already to do so
  ierr = warnUserOptionsIgnored(fname); CHKERRQ(ierr);  
  ierr = grid.createDA(); CHKERRQ(ierr);
  // FIXME: note we *can* determine from the input file whether the hor. dims are truely periodic,
  // but this has not been done; here we simply require it is not periodic
  ierr = grid.rescale_using_zlevels(-g.x_min, -g.y_min); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);
  
  // set IceModel::startYear, IceModel::endYear, grid.year, but respecting grid.year
  // which came from -if file, _unless_ -ys set by user
  ierr = setStartRunEndYearsFromOptions(PETSC_TRUE);  CHKERRQ(ierr);

  void *a_mpi;
  int a_len, max_a_len;
  max_a_len = a_len = grid.xm * grid.ym * grid.Mz;
  MPI_Reduce(&a_len, &max_a_len, 1, MPI_INT, MPI_MAX, 0, grid.com);
  ierr = PetscMalloc(max_a_len * sizeof(double), &a_mpi); CHKERRQ(ierr);

  // 2-D mapping
  ierr = vLongitude.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr =  vLatitude.read(fname, g.t_len - 1); CHKERRQ(ierr);

  // 2-D model quantities: discrete
  ierr = vMask.read(fname, g.t_len - 1); CHKERRQ(ierr);

  // 2-D model quantities: double
  ierr =  vHmelt.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr =      vH.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr =    vbed.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr = vuplift.read(fname, g.t_len - 1); CHKERRQ(ierr);

  // Read vubarSSA and vvbarSSA if SSA is on, if not asked to ignore them and
  // if they are present in the input file.
  PetscTruth dontreadSSAvels = PETSC_FALSE;
  ierr = PetscOptionsHasName(PETSC_NULL, "-dontreadSSAvels", &dontreadSSAvels); CHKERRQ(ierr);

  if ((grid.rank == 0) && useSSAVelocity && (!dontreadSSAvels)) {

    int flag;
    stat = nc_get_att_int(nc.ncid, NC_GLOBAL, "ssa_velocities_are_valid", &flag);
    if (stat == NC_NOERR)
      have_ssa_velocities = flag;
  }
  ierr = MPI_Bcast(&have_ssa_velocities, 1, MPI_INT, 0, grid.com); CHKERRQ(ierr);
  
  if (have_ssa_velocities == 1) {
    ierr = verbPrintf(2,grid.com,"Reading vubarSSA and vvbarSSA...\n"); CHKERRQ(ierr);

    ierr = vubarSSA.read(fname, g.t_len - 1);
    ierr = vvbarSSA.read(fname, g.t_len - 1);
  }

  // 2-D climate/bdry quantities
  ierr =      vTs.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr =     vGhf.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr =   vAccum.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr = vtillphi.read(fname, g.t_len - 1); CHKERRQ(ierr);

  // 3-D model quantities
  ierr =   T3.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr =  Tb3.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr = tau3.read(fname, g.t_len - 1); CHKERRQ(ierr);

  ierr = PetscFree(a_mpi); CHKERRQ(ierr);

  // read the polar_stereographic if present
  ierr = nc.read_polar_stereographic(psParams.svlfp,
				     psParams.lopo,
				     psParams.sp); CHKERRQ(ierr);

  // Get the current history length
  unsigned int history_len;	// used for communication
  if (grid.rank == 0) {
    size_t H;
    stat = nc_inq_attlen(nc.ncid, NC_GLOBAL, "history", &H);
    history_len = (int)H;
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  // Broadcast the history length
  MPI_Bcast(&history_len, 1, MPI_INT, 0, grid.com);

  // Allocate some memory (if necessary)
  if (history_len > history_size - 1) {
    history_size = history_len + 1;
    delete[] history;
    history = new char[history_size];
  }

  // Zero out the allocated memory so that we don't need to worry about the
  // trailing zero in the history string (which might be absent).
  memset(history, 0, history_size);

  // Read the history string and close the file
  if (grid.rank == 0) {
    stat = nc_get_att_text(nc.ncid, NC_GLOBAL, "history", history);
    CHKERRQ(check_err(stat,__LINE__,__FILE__));
  }

  ierr = nc.close(); CHKERRQ(ierr);

  // Broadcast the string
  MPI_Bcast(history, history_size, MPI_CHAR, 0, grid.com);

  initialized_p = PETSC_TRUE;
  return 0;
}


//! Manage regridding based on user options.  Call IceModelVec::regrid() to do each selected variable.
/*!
For each variable selected by option <tt>-regrid_vars</tt>, we regrid it onto the current grid from 
the NetCDF file specified by <tt>-regrid</tt>.

The default, if <tt>-regrid_vars</tt> is not given, is to regrid the 3 dimensional 
quantities \c tau3, \c T3, \c Tb3.  This is consistent with one standard purpose of 
regridding, which is to stick with current geometry through the downscaling procedure.  
Most of the time the user should carefully specify which variables to regrid.
 */
PetscErrorCode IceModel::regrid(const char *filename) {
  PetscErrorCode ierr;
  PetscTruth regridVarsSet;
  char regridVars[PETSC_MAX_PATH_LEN];
  NCTool nc(&grid);

  const int  npossible = 7;
  const char possible[20] = "bBehHLT";
  
  ierr = PetscOptionsGetString(PETSC_NULL, "-regrid_vars", regridVars,
                               PETSC_MAX_PATH_LEN, &regridVarsSet); CHKERRQ(ierr);
  if (regridVarsSet == PETSC_FALSE) {
    strcpy(regridVars, "");
  }
  ierr = verbPrintf(2,grid.com, 
           "regridding variables with single character flags `%s' from NetCDF file `%s':", 
           regridVars,filename); CHKERRQ(ierr);

  // following are dimensions, limits and lengths, and id for *source* NetCDF file (regridFile)

  // create "local interpolation context" from dimensions, limits, and lengths extracted from regridFile,
  //   and from information about the part of the grid owned by this processor

  bool file_exists = false;
  ierr = nc.open_for_reading(filename, file_exists);
  if (!file_exists) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: Can't open file '%s'.\n", filename); CHKERRQ(ierr);
    PetscEnd();
  }
  
  grid_info g;
  ierr = nc.get_grid_info(g); CHKERRQ(ierr);

  double *zlevs, *zblevs;
  zlevs = new double[g.z_len];
  zblevs = new double[g.zb_len];
  ierr = nc.get_vertical_dims(g.z_len, g.zb_len, zlevs, zblevs); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  { // explicit scoping means destructor will be called for lic
    LocalInterpCtx lic(g, zlevs, zblevs, grid);
    // ierr = lic.printGrid(grid.com); CHKERRQ(ierr);
    
    for (PetscInt k = 0; k < npossible; k++) {
      if (strchr(regridVars, possible[k])) {
       switch (possible[k]) {
         case 'b':
           ierr = verbPrintf(2, grid.com, "\n   b: regridding 'topg' ... \n"); CHKERRQ(ierr);
	   ierr = vbed.regrid(filename, lic, true); CHKERRQ(ierr);
           break;
         case 'B':
           ierr = verbPrintf(2, grid.com, "\n   B: regridding 'litho_temp' ... \n"); CHKERRQ(ierr);
           ierr = Tb3.regrid(filename, lic, true); CHKERRQ(ierr);
           break;
         case 'e':
           ierr = verbPrintf(2, grid.com, "\n   e: regridding 'age' ... \n"); CHKERRQ(ierr);
           ierr = tau3.regrid(filename, lic, true); CHKERRQ(ierr);
           break;
         case 'h':
           ierr = verbPrintf(2, grid.com, "\n   h: regridding 'usurf' ... \n"); CHKERRQ(ierr);
	   ierr = vh.regrid(filename, lic, true); CHKERRQ(ierr);
           break;
         case 'H':
           ierr = verbPrintf(2, grid.com, "\n   H: regridding 'thk' ... \n"); CHKERRQ(ierr);
	   ierr = vH.regrid(filename, lic, true); CHKERRQ(ierr);
           break;
         case 'L':
           ierr = verbPrintf(2, grid.com, "\n   L: regridding 'bwat' ... \n"); CHKERRQ(ierr);
	   ierr = vHmelt.regrid(filename, lic, true); CHKERRQ(ierr);
           break;
         case 'T':
           ierr = verbPrintf(2, grid.com, "\n   T: regridding 'temp' ... \n"); CHKERRQ(ierr);
           ierr = T3.regrid(filename, lic, true); CHKERRQ(ierr);
           break;
       }
      }
    }

  }
  
  delete [] zlevs;  delete [] zblevs;
  ierr = verbPrintf(2,grid.com, "\n"); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode IceModel::init_snapshots_from_options() {
  PetscErrorCode ierr;
  PetscTruth save_at_set = PETSC_FALSE, save_to_set = PETSC_FALSE;
  char tmp[TEMPORARY_STRING_LENGTH];
  first_snapshot = next_snapshot = last_snapshot = snapshot_dt = 0;

  ierr = PetscOptionsGetString(PETSC_NULL, "-save_to", snapshots_filename,
			       PETSC_MAX_PATH_LEN, &save_to_set); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-save_at", tmp,
			       TEMPORARY_STRING_LENGTH, &save_at_set); CHKERRQ(ierr);

  if (save_to_set && !save_at_set) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: -save_to is set, but -save_at is not.\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

  if (save_at_set && !save_to_set) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: -save_at is set, but -save_to is not.\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

  // check if the user specified a MATLAB-style range (a:dt:b). We assume that
  // any string containing a colon defines a range like this
  char *endptr1, *endptr2;
  bool parsing_failed = false;
  if (strchr(tmp, ':')) {	// if a string contains a colon...
      
    // try to read the first number
    first_snapshot = strtod(tmp, &endptr1);
    if (tmp == endptr1)
      parsing_failed = true;

    endptr1++;			// skip the colon separating a and dt (a:dt)
    // try to read the second number
    snapshot_dt = strtod(endptr1, &endptr2);
    if (endptr1 == endptr2)
      parsing_failed = true;

    endptr2++;			// skip the colon separating dt and b (dt:b)
    // try to read the third number
    last_snapshot = strtod(endptr2, &endptr1);
    if (endptr1 == endptr2)
      parsing_failed = true;

    if (parsing_failed) {
      ierr = PetscPrintf(grid.com, "PISM ERROR: Parsing the -save_at argument failed.\n");
         CHKERRQ(ierr);
      PetscEnd();
    }

    if (first_snapshot >= last_snapshot) {
      ierr = PetscPrintf(grid.com,
         "PISM ERROR: Error in the -save_at argument: a >= b in the range specification '%s'.\n",
	 tmp); CHKERRQ(ierr);
      PetscEnd();
    }

    if (snapshot_dt <= 0) {
      ierr = PetscPrintf(grid.com,
	 "PISM ERROR: Error in the -save_at argument: dt <= 0 in the range specification '%s'.\n",
	 tmp); CHKERRQ(ierr);
      PetscEnd();
    }

    save_at_equal_intervals = PETSC_TRUE;
    next_snapshot = first_snapshot;
    last_snapshot += snapshot_dt; // make it save after the last time in the range, too
  } else {			// no colon; it must be a list of numbers
    n_snapshots = max_n_snapshots;
    ierr = PetscOptionsGetRealArray(PETSC_NULL, "-save_at", save_at, &n_snapshots, PETSC_NULL);
    qsort(save_at, n_snapshots, sizeof(double), compare_doubles); // sort the values
    save_at_equal_intervals = PETSC_FALSE;
    current_snapshot = 0;
  }

  if (save_to_set && save_at_set) {
    save_snapshots = true;
    file_is_ready = false;
    split_snapshots = false;

    PetscTruth split;
    ierr = PetscOptionsHasName(PETSC_NULL, "-split_snapshots", &split); CHKERRQ(ierr);
    if (split) {
      split_snapshots = true;
    } else if (!hasSuffix(snapshots_filename, ".nc")) {
      ierr = verbPrintf(2, grid.com,
			"PISM WARNING: snapshots file name does not have the '.nc' suffix!\n");
      CHKERRQ(ierr);
    }

    if (split) {
      ierr = verbPrintf(2, grid.com, "saving snapshots to '%s'+year.nc; ", snapshots_filename); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid.com, "saving snapshots to '%s'; ", snapshots_filename); CHKERRQ(ierr);
    }

    if (save_at_equal_intervals) {
      ierr = verbPrintf(2, grid.com, "times requested: %3.3f:%3.3f:%3.3f\n", first_snapshot, snapshot_dt, last_snapshot); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid.com, "times requested: "); CHKERRQ(ierr);
      for (int j = 0; j < n_snapshots; j++) {
	ierr = verbPrintf(2, grid.com, "%3.3f, ", save_at[j]); CHKERRQ(ierr);
      }
      ierr = verbPrintf(2, grid.com, "\b\b\n"); CHKERRQ(ierr);
    }
  } else {
    save_snapshots = false;
  }

  return 0;
}

PetscErrorCode IceModel::write_snapshot() {
  PetscErrorCode ierr;
  NCTool nc(&grid);
  bool save_now = false;
  double saving_after;
  char filename[PETSC_MAX_PATH_LEN];

  // determine if the user set the -save_at and -save_to options
  if (!save_snapshots)
    return 0;

  // do we need to save *now*?
  if (save_at_equal_intervals) {
    if (grid.year >= next_snapshot && grid.year < last_snapshot) {
      save_now = true;
      saving_after = next_snapshot;

      while (next_snapshot <= grid.year)
	next_snapshot += snapshot_dt;
    }
  } else if ( (grid.year >= save_at[current_snapshot]) && (current_snapshot < n_snapshots) ) {
    save_now = true;
    saving_after = save_at[current_snapshot];

    while (save_at[current_snapshot] <= grid.year)
      current_snapshot++;
  }

  if (save_now) {
    if (split_snapshots) {
      file_is_ready = false;	// each snapshot is written to a separate file
      snprintf(filename, PETSC_MAX_PATH_LEN, "%s-%06.0f.nc", snapshots_filename, grid.year);
    } else {
      strncpy(filename, snapshots_filename, PETSC_MAX_PATH_LEN);
    }

    ierr = verbPrintf(2, grid.com, "Saving a model state snapshot at t = %3.5f.\n", grid.year);
    CHKERRQ(ierr);

    char tmp[TEMPORARY_STRING_LENGTH];
    snprintf(tmp, TEMPORARY_STRING_LENGTH,
	     " %s: saving a snapshot at %10.5f (years), following time-step at %10.5f (years).\n",
	     executable_short_name, grid.year, saving_after);

    if (!file_is_ready) {
      // Prepare the snapshots file:
      ierr = nc.open_for_writing(filename, true); CHKERRQ(ierr);
      ierr = nc.write_history(history); CHKERRQ(ierr); // append the history
      ierr = nc.write_polar_stereographic(psParams.svlfp, psParams.lopo, psParams.sp); CHKERRQ(ierr);
      ierr = nc.write_global_attrs(useSSAVelocity, "CF-1.0"); CHKERRQ(ierr);
      ierr = nc.close(); CHKERRQ(ierr);
      file_is_ready = true;
    }
    
    ierr = nc.open_for_writing(filename, false); CHKERRQ(ierr); // replace == false
    ierr = nc.append_time(grid.year * secpera); CHKERRQ(ierr);
    ierr = nc.write_history(tmp); CHKERRQ(ierr); // append the history
    ierr = nc.close(); CHKERRQ(ierr);

    ierr = write_model_state(filename);  CHKERRQ(ierr);
    ierr = write_extra_fields(filename); CHKERRQ(ierr);
  }

  return 0;
}
