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
  size_t flen = strlen(fname);
  size_t slen = strlen(suffix);
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
		      "Writing model state, with full 3D velocities, to file `%s'\n",
		      filename); CHKERRQ(ierr);

    ierr = dumpToFile(filename); CHKERRQ(ierr);
    ierr = write3DPlusToFile(filename); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com, "Writing model state to file `%s'\n",
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
  ierr = nc.write_global_attrs(useSSAVelocity, "CF-1.3"); CHKERRQ(ierr);
  ierr = nc.close(); CHKERRQ(ierr);

  ierr = write_model_state(filename);  CHKERRQ(ierr);

  if (atmosPCC != PETSC_NULL) {
    ierr = atmosPCC->writeCouplingFieldsToFile(filename); CHKERRQ(ierr);
  } else {
    SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL");
  }

  ierr = write_extra_fields(filename); CHKERRQ(ierr); // chance for derived classes to do more

  return 0;
}


PetscErrorCode IceModel::write_model_state(const char filename[]) {
  PetscErrorCode ierr;

  // 2D model quantities
  ierr =         vH.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  ierr = vLongitude.write(filename, NC_DOUBLE); CHKERRQ(ierr);
  ierr =  vLatitude.write(filename, NC_DOUBLE); CHKERRQ(ierr);

  ierr = vMask.write(filename, NC_BYTE); CHKERRQ(ierr);
  double mask_values[4] = {MASK_SHEET, MASK_DRAGGING, MASK_FLOATING, MASK_FLOATING_OCEAN0};
  ierr = vMask.write_scalar_attr(filename, "flag_values", NC_BYTE, 4, mask_values); CHKERRQ(ierr);
  ierr = vMask.write_text_attr(filename, "flag_meanings",
			       "sheet dragging floating floating_at_time_0"); CHKERRQ(ierr);

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

  // 2D earth quantity; like climate but always steady and always owned by IceModel
  ierr =   vGhf.write(filename, NC_DOUBLE); CHKERRQ(ierr);

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
  ierr = vWork2d[0].set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vWork2d[0].write_in_glaciological_units = true;
  vWork2d[0].set_valid_min(0.0);
  ierr = vWork2d[0].write(filename, NC_FLOAT); CHKERRQ(ierr);

  // compute cflx = cbar .* thk and save it
  ierr = vWork2d[0].multiply_by(vH, vWork2d[1]); CHKERRQ(ierr);
  ierr = vWork2d[1].set_name("cflx"); CHKERRQ(ierr);
  ierr = vWork2d[1].set_attrs("diagnostic", 
             "magnitude of vertically-integrated horizontal flux of ice",
	     "m2 s-1", NULL); CHKERRQ(ierr);
  ierr = vWork2d[1].set_glaciological_units("m2 year-1"); CHKERRQ(ierr);
  vWork2d[1].write_in_glaciological_units = true;
  vWork2d[1].set_valid_min(0.0);
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
  ierr = vWork2d[2].set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vWork2d[2].write_in_glaciological_units = true;
  vWork2d[2].set_valid_min(0.0);
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
  ierr = vWork2d[0].set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vWork2d[0].write_in_glaciological_units = true;
  vWork2d[0].set_valid_min(0.0);
  ierr = vWork2d[0].write(filename, NC_FLOAT); CHKERRQ(ierr);

  // compute wvelsurf, the surface values of vertical velocity
  ierr = w3.begin_access(); CHKERRQ(ierr);
  ierr = w3.getSurfaceValues(vWork2d[0], vH); CHKERRQ(ierr);
  ierr = w3.end_access(); CHKERRQ(ierr);

  ierr = vWork2d[0].set_name("wvelsurf"); CHKERRQ(ierr);
  ierr = vWork2d[0].set_attrs("diagnostic", "vertical velocity of ice at ice surface",
			      "m s-1", NULL); CHKERRQ(ierr);
  ierr = vWork2d[0].set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vWork2d[0].write_in_glaciological_units = true;
  ierr = vWork2d[0].write(filename, NC_FLOAT); CHKERRQ(ierr);

  // compute magnitude of basal shear stress = rho g H |grad h|
  ierr = computeDrivingStress(vWork2d[0],vWork2d[1]); CHKERRQ(ierr);
  ierr = getMagnitudeOf2dVectorField(vWork2d[0],vWork2d[1],vWork2d[2]); CHKERRQ(ierr);

  ierr = vWork2d[2].set_name("taud"); CHKERRQ(ierr);
  ierr = vWork2d[2].set_attrs("diagnostic",
             "magnitude of driving shear stress at base of ice",
	     "Pa", NULL); CHKERRQ(ierr);
  vWork2d[2].set_valid_min(0.0);
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
  ierr = vWork2d[0].set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vWork2d[0].write_in_glaciological_units = true;
  ierr = vWork2d[0].write(filename, NC_FLOAT); CHKERRQ(ierr);

  ierr = v3.begin_access(); CHKERRQ(ierr);
  ierr = v3.getSurfaceValues(vWork2d[0], vH); CHKERRQ(ierr);
  ierr = v3.end_access(); CHKERRQ(ierr);

  ierr = vWork2d[0].set_name("vvelsurf"); CHKERRQ(ierr);
  ierr = vWork2d[0].set_attrs(
              "diagnostic", "y component of velocity of ice at ice surface",
	      "m s-1", NULL); CHKERRQ(ierr);
  ierr = vWork2d[0].set_glaciological_units("m year-1"); CHKERRQ(ierr);
  vWork2d[0].write_in_glaciological_units = true;
  ierr = vWork2d[0].write(filename, NC_FLOAT); CHKERRQ(ierr);

  return 0;
}


//! When reading a saved PISM model state, warn the user if options <tt>-Mx,-My,-Mz,-Mbz</tt> have been ignored.
PetscErrorCode IceModel::warnUserOptionsIgnored(const char *fname) {
  PetscErrorCode ierr;
  PetscInt       ignore;
  PetscTruth     M_Set;

  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mx", &ignore, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,
             "WARNING: user option -Mx ignored; value read from file %s\n", fname);
             CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetInt(PETSC_NULL, "-My", &ignore, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,
             "WARNING: user option -My ignored; value read from file %s\n", fname);
             CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mz", &ignore, &M_Set); CHKERRQ(ierr);
  if (M_Set == PETSC_TRUE) {
    ierr = verbPrintf(1,grid.com,
             "WARNING: user option -Mz ignored; value read from file %s\n", fname);
             CHKERRQ(ierr);
  }
  ierr = PetscOptionsGetInt(PETSC_NULL, "-Mbz", &ignore, &M_Set); CHKERRQ(ierr);
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
  grid.Mx   = g.x_len;
  grid.My   = g.y_len;
  grid.Mz   = g.z_len;
  grid.Mbz  = g.zb_len;
  grid.x0   = g.x0;
  grid.y0   = g.y0;
  // grid.Lx, grid.Ly are set from g.Lx, g.Ly below in call to grid.rescale_using_zlevels()

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
  ierr = grid.rescale_using_zlevels(g.Lx, g.Ly); CHKERRQ(ierr);
  ierr = createVecs(); CHKERRQ(ierr);

  // set IceModel::startYear, IceModel::endYear, grid.year, but respecting grid.year
  // which came from -i file, _unless_ -ys set by user
  ierr = setStartRunEndYearsFromOptions(PETSC_TRUE);  CHKERRQ(ierr);
 
  // 2-D mapping
  ierr = vLongitude.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr =  vLatitude.read(fname, g.t_len - 1); CHKERRQ(ierr);

  // 2-D model quantities: discrete
  ierr = vMask.read(fname, g.t_len - 1); CHKERRQ(ierr);

  // 2-D model quantities: double
  ierr =   vHmelt.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr =       vH.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr =     vbed.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr =  vuplift.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr = vtillphi.read(fname, g.t_len - 1); CHKERRQ(ierr);


  if (useSSAVelocity) {
    double flag;
    stat = nc.get_att_double(NC_GLOBAL, "ssa_velocities_are_valid", 1, &flag);
    if (stat == 0)
      have_ssa_velocities = (int)flag;
  }

  // Read vubarSSA and vvbarSSA if SSA is on, if not asked to ignore them and
  // if they are present in the input file.
  PetscTruth dontreadSSAvels = PETSC_FALSE;
  ierr = PetscOptionsHasName(PETSC_NULL, "-dontreadSSAvels", &dontreadSSAvels); CHKERRQ(ierr);
  
  if ((have_ssa_velocities == 1)  && (!dontreadSSAvels)) {
    ierr = verbPrintf(3,grid.com,"Reading vubarSSA and vvbarSSA...\n"); CHKERRQ(ierr);

    ierr = vubarSSA.read(fname, g.t_len - 1);
    ierr = vvbarSSA.read(fname, g.t_len - 1);
  }

  // 2-D earth quantity; like climate but owned by IceModel
  ierr =     vGhf.read(fname, g.t_len - 1); CHKERRQ(ierr);

  // 3-D model quantities
  ierr =   T3.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr =  Tb3.read(fname, g.t_len - 1); CHKERRQ(ierr);
  ierr = tau3.read(fname, g.t_len - 1); CHKERRQ(ierr);

  // read the polar_stereographic if present
  ierr = nc.read_polar_stereographic(psParams.svlfp,
				     psParams.lopo,
				     psParams.sp); CHKERRQ(ierr);

  int hist_len;
  char *hist;
  stat = nc.get_att_text(NC_GLOBAL, "history", &hist_len, &hist);
  if (hist != NULL) {
    delete[] history;
    history = hist;
    history_size = hist_len;
  }

  ierr = nc.close(); CHKERRQ(ierr);

  initialized_p = PETSC_TRUE;
  return 0;
}


//! Manage regridding based on user options.  Call IceModelVec::regrid() to do each selected variable.
/*!
For each variable selected by option <tt>-regrid_vars</tt>, we regrid it onto the current grid from 
the NetCDF file specified by <tt>-regrid_from</tt>.

The default, if <tt>-regrid_vars</tt> is not given, is to regrid the 3 dimensional 
quantities \c tau3, \c T3, \c Tb3.  This is consistent with one standard purpose of 
regridding, which is to stick with current geometry through the downscaling procedure.  
Most of the time the user should carefully specify which variables to regrid.
 */
PetscErrorCode IceModel::regrid(const char *filename) {
  PetscErrorCode ierr;
  PetscTruth regridVarsSet;
  bool regrid_2d_only = false;
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
           "regridding variables with single character flags `%s' from NetCDF file `%s':\n", 
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
  // Note that after this call g.z_len and g.zb_len are zero if the
  // corresponding dimension does not exist.
  ierr = nc.get_grid_info(g); CHKERRQ(ierr);

  double *zlevs = NULL, *zblevs = NULL; // NULLs correspond to 2D-only regridding
  if ((g.z_len != 0) && (g.zb_len != 0)) {
    zlevs  = new double[g.z_len];
    zblevs = new double[g.zb_len];
    ierr = nc.get_vertical_dims(g.z_len, g.zb_len, zlevs, zblevs); CHKERRQ(ierr);
  } else {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: at least one of 'z' and 'zb' is absent in '%s'.\n"
		      "              3D regridding is disabled.\n",
		      filename);
    regrid_2d_only = true;
    CHKERRQ(ierr);
  }
  ierr = nc.close(); CHKERRQ(ierr);

  { // explicit scoping means destructor will be called for lic
    LocalInterpCtx lic(g, zlevs, zblevs, grid);
    // ierr = lic.printGrid(grid.com); CHKERRQ(ierr);
    
    for (PetscInt k = 0; k < npossible; k++) {
      if (strchr(regridVars, possible[k])) {
       switch (possible[k]) {
         case 'b':
           ierr = verbPrintf(2, grid.com, "  b: regridding 'topg' ... \n"); CHKERRQ(ierr);
	   ierr = vbed.regrid(filename, lic, true); CHKERRQ(ierr);
           break;
         case 'B':
	   if (regrid_2d_only) {
	     ierr = verbPrintf(2, grid.com, "  B: skipping 'litho_temp'...\n"); CHKERRQ(ierr);
	   } else {
	     ierr = verbPrintf(2, grid.com, "  B: regridding 'litho_temp' ... \n"); CHKERRQ(ierr);
	     ierr = Tb3.regrid(filename, lic, true); CHKERRQ(ierr);
	   }
           break;
         case 'e':
	   if (regrid_2d_only) {
	     ierr = verbPrintf(2, grid.com, "  e: skipping 'age'...\n"); CHKERRQ(ierr);
	   } else {
	     ierr = verbPrintf(2, grid.com, "  e: regridding 'age' ... \n"); CHKERRQ(ierr);
	     ierr = tau3.regrid(filename, lic, true); CHKERRQ(ierr);
	   }
           break;
         case 'h':
           ierr = verbPrintf(2, grid.com, "  h: regridding 'usurf' ... \n"); CHKERRQ(ierr);
	   ierr = vh.regrid(filename, lic, true); CHKERRQ(ierr);
           break;
         case 'H':
           ierr = verbPrintf(2, grid.com, "  H: regridding 'thk' ... \n"); CHKERRQ(ierr);
	   ierr = vH.regrid(filename, lic, true); CHKERRQ(ierr);
           break;
         case 'L':
           ierr = verbPrintf(2, grid.com, "  L: regridding 'bwat' ... \n"); CHKERRQ(ierr);
	   ierr = vHmelt.regrid(filename, lic, true); CHKERRQ(ierr);
           break;
         case 'T':
	   if (regrid_2d_only) {
	     ierr = verbPrintf(2, grid.com, "  T: skipping 'temp'...\n"); CHKERRQ(ierr);
	   } else {
	     ierr = verbPrintf(2, grid.com, "  T: regridding 'temp' ... \n"); CHKERRQ(ierr);
	     ierr = T3.regrid(filename, lic, true); CHKERRQ(ierr);
	   }
           break;
       }
      }
    }

  }

  // Note that deleting a null pointer is safe.
  delete [] zlevs;  delete [] zblevs;
  return 0;
}


PetscErrorCode IceModel::init_snapshots_from_options() {
  PetscErrorCode ierr;
  PetscTruth save_at_set = PETSC_FALSE, save_to_set = PETSC_FALSE;
  char tmp[TEMPORARY_STRING_LENGTH] = "\0";
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
    last_snapshot += snapshot_dt; // make it saves after the last time in the range, too
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
      ierr = verbPrintf(2, grid.com, "saving snapshots to '%s+year.nc'; ", snapshots_filename); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid.com, "saving snapshots to '%s'; ", snapshots_filename); CHKERRQ(ierr);
    }

    if (save_at_equal_intervals) {
      ierr = verbPrintf(2, grid.com, "times requested: %3.3f:%3.3f:%3.3f\n",
			first_snapshot, snapshot_dt, last_snapshot - snapshot_dt); CHKERRQ(ierr);
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

    ierr = verbPrintf(2, grid.com, 
       "Saving a model state snapshot at %3.5f (years), the time-step after goal %3.5f.\n",
       grid.year,saving_after);
    CHKERRQ(ierr);

    char tmp[TEMPORARY_STRING_LENGTH];
    snprintf(tmp, TEMPORARY_STRING_LENGTH,
       " %s: saving model state snapshot at %10.5f (years), the time-step after goal %10.5f.\n",
       executable_short_name, grid.year, saving_after);

    if (!file_is_ready) {
      // Prepare the snapshots file:
      ierr = nc.open_for_writing(filename, true); CHKERRQ(ierr);
      ierr = nc.write_history(history); CHKERRQ(ierr); // append the history
      ierr = nc.write_polar_stereographic(psParams.svlfp, psParams.lopo, psParams.sp); CHKERRQ(ierr);
      ierr = nc.write_global_attrs(useSSAVelocity, "CF-1.3"); CHKERRQ(ierr);
      ierr = nc.close(); CHKERRQ(ierr);
      file_is_ready = true;
    }
    
    ierr = nc.open_for_writing(filename, false); CHKERRQ(ierr); // replace == false
    ierr = nc.append_time(grid.year * secpera); CHKERRQ(ierr);
    ierr = nc.write_history(tmp); CHKERRQ(ierr); // append the history
    ierr = nc.close(); CHKERRQ(ierr);

    ierr = write_model_state(filename);  CHKERRQ(ierr);

    if (atmosPCC != PETSC_NULL) {
      ierr = atmosPCC->writeCouplingFieldsToFile(filename); CHKERRQ(ierr);
    } else {
      SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL");
    }
    
    ierr = write_extra_fields(filename); CHKERRQ(ierr);
  }

  return 0;
}
