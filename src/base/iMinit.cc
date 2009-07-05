// Copyright (C) 2009 Ed Bueler and Constantine Khroulev
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

//This file contains various initialization routines. See the IceModel::init()
//documentation comment in iceModel.cc for the order in which they are called.

#include <petscda.h>
#include "iceModel.hh"

//! Set default values of grid parameters.
/*!
  Derived classes (IceCompModel, for example) reimplement this to change the
  grid initialization when no -i option is set.
 */
PetscErrorCode IceModel::set_grid_defaults() {
  PetscErrorCode ierr;
  PetscTruth Mx_set, My_set, Mz_set, Lz_set, boot_from_set;
  char filename[PETSC_MAX_PATH_LEN];
  grid_info gi;

  // Get the bootstrapping file name:
  ierr = check_old_option_and_stop(grid.com, "-bif", "-boot_from"); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-boot_from",
			       filename, PETSC_MAX_PATH_LEN, &boot_from_set); CHKERRQ(ierr);

  if (!boot_from_set) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: Please specify an input file using -i or -boot_from.\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

  // Use a bootstrapping file to set some grid parameters (they can be
  // overridden later, in IceModel::set_grid_from_options()).

  // Determine the grid extent from a bootstrapping file:
  NCTool nc(&grid);
  bool x_dim_exists, y_dim_exists, t_exists;
  ierr = nc.open_for_reading(filename); CHKERRQ(ierr);

  ierr = nc.find_dimension("x", NULL, x_dim_exists); CHKERRQ(ierr);
  ierr = nc.find_dimension("y", NULL, y_dim_exists); CHKERRQ(ierr);
  ierr = nc.find_variable("t", NULL, t_exists); CHKERRQ(ierr);
  ierr = nc.get_grid_info(gi);
  ierr = nc.close(); CHKERRQ(ierr);

  // if the horizontal dimensions are absent then we can not proceed
  if (!x_dim_exists) {
    ierr = PetscPrintf(grid.com,"bootstrapping file '%s' has no horizontal dimension 'x'\n",filename);
    CHKERRQ(ierr);
    PetscEnd();
  }
  if (!y_dim_exists) {
    ierr = PetscPrintf(grid.com,"bootstrapping file '%s' has no horizontal dimension 'y'\n",filename);
    CHKERRQ(ierr);
    PetscEnd();
  }

  // Set the grid center and horizontal extent:
  grid.x0 = gi.x0;
  grid.y0 = gi.y0;
  grid.Lx = gi.Lx;
  grid.Ly = gi.Ly;

  if (t_exists) {
    grid.year = gi.time / secpera; // set year from read-in time variable
    ierr = verbPrintf(2, grid.com, 
		      "  time t = %5.4f years found; setting current year\n",
		      grid.year); CHKERRQ(ierr);
  } else {
    grid.year = 0.0;
    ierr = verbPrintf(2, grid.com, 
		      "  time dimension was not found; setting current year to t = 0.0 years\n",
		      grid.year); CHKERRQ(ierr);
  }

  // Grid dimensions and its vertical extent should not be deduced from a
  // bootstrapping file, so we check if these options are set and stop if they
  // are not.
  // Note that here interpreting "-Mx 0" as "-Mx was not set" is OK.
  ierr = check_option("-Mx", Mx_set); CHKERRQ(ierr);
  ierr = check_option("-My", My_set); CHKERRQ(ierr);
  ierr = check_option("-Mz", Mz_set); CHKERRQ(ierr);
  ierr = check_option("-Lz", Lz_set); CHKERRQ(ierr);
  if ( !(Mx_set && My_set && Mz_set && Lz_set) ) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: All of -boot_from, -Mx, -My, -Mz, -Lz, are required for bootstrapping.\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

  return 0;
}

//! Initalizes the grid from options.
/*! 
  Reads all of -Mx, -My, -Mz, -Mbz, -Lx, -Ly, -Lz, -quadZ and -chebZ. Sets
  corresponding grid parameters.
 */
PetscErrorCode IceModel::set_grid_from_options() {
  PetscErrorCode ierr;
  PetscTruth Mx_set, My_set, Mz_set, Mbz_set, Lx_set, Ly_set, Lz_set,
    quadZ_set, chebZ_set;
  PetscScalar x_scale, y_scale, z_scale;
  int Mx, My, Mz, Mbz;

  // Process the options:
  ierr = PetscOptionsBegin(grid.com, PETSC_NULL,
			   "Options setting the computational grid extent and dimensions",
			   PETSC_NULL); CHKERRQ(ierr);

  // Read -Lx and -Ly. Note the transpose!
  ierr = PetscOptionsScalar("-Lx", "Half of the grid extent in the X direction, in km", "",
			    y_scale, &y_scale, &Ly_set); CHKERRQ(ierr);
  ierr = PetscOptionsScalar("-Ly", "Half of the grid extent in the Y direction, in km", "",
			    x_scale, &x_scale, &Lx_set); CHKERRQ(ierr);
  // Vertical extent (in the ice):
  ierr = PetscOptionsScalar("-Lz", "Grid extent in the Z (vertical) direction in the ice", "",
			    z_scale, &z_scale, &Lz_set); CHKERRQ(ierr);

  // Read -Mx, -My, -Mz and -Mbz. Note the transpose!
  ierr = PetscOptionsInt("-Mx", "Number of grid points in the X direction", "",
			 grid.My, &My, &My_set); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-My", "Number of grid points in the Y direction", "",
			 grid.Mx, &Mx, &Mx_set); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-Mz", "Number of grid points in the Z (vertical) direction in the ice", "",
			 grid.Mz, &Mz, &Mz_set); CHKERRQ(ierr);
  ierr = PetscOptionsInt("-Mbz", "Number of grid points in the Z (vertical) direction in the bedrock", PETSC_NULL,
			 grid.Mbz, &Mbz, &Mbz_set); CHKERRQ(ierr);

  // Determine the vertical grid spacing in the ice:
  ierr = PetscOptionsName("-quadZ", "Chooses the quadratic vertical grid spacing",
			  PETSC_NULL, &quadZ_set); CHKERRQ(ierr);
  ierr = PetscOptionsName("-chebZ", "Chooses the Chebyshev vertical grid spacing",
			  PETSC_NULL, &chebZ_set); CHKERRQ(ierr);

  // Only one of -quadZ and -chebZ is allowed.
  if (quadZ_set && chebZ_set) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: at most one of -quadZ and -chebZ is allowed.\n"); CHKERRQ(ierr);
    PetscEnd();
  }

  // Done with the options.
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  // Use the information obtained above:
  if (Lx_set)    grid.Lx = x_scale * 1000.0; // convert to meters
  if (Ly_set)    grid.Ly = y_scale * 1000.0; // convert to meters
  if (Lz_set)    grid.Lz = z_scale;	     // in meters already
  if (Mx_set)    grid.Mx = Mx;
  if (My_set)    grid.My = My;
  if (Mz_set)    grid.Mz = Mz;
  if (Mbz_set)   grid.Mbz = Mbz;
  if (quadZ_set) grid.vertical_spacing = QUADRATIC;
  if (chebZ_set) grid.vertical_spacing = CHEBYSHEV;

  ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
  ierr = grid.compute_vertical_levels();    CHKERRQ(ierr);

  // At this point all the fields except for da2, xs, xm, ys, ym should be
  // filled. We're ready to call grid.createDA().
  return 0;
}

//! Sets up the computational grid.
/*!
  There are two cases here:

  1) Initializing from a PISM ouput file, in which case all the options
  influencing the grid (currently: -Mx, -My, -Mz, -Mbz, -Lx, -Ly, -Lz, -quadZ,
  -chebZ) are ignored.

  2) Initializing using defaults, command-line options and (possibly) a
  bootstrapping file. Derived classes requiring special grid setup should
  reimplement IceGrid::set_grid_from_options().

  No memory allocation should happen here.
 */
PetscErrorCode IceModel::grid_setup() {
  PetscErrorCode ierr;
  PetscTruth i_set;
  char filename[PETSC_MAX_PATH_LEN];

  ierr = verbPrintf(3, grid.com,
		    "Setting up the computational grid...\n"); CHKERRQ(ierr);

  // Check if we are initializing from a PISM output file:
  ierr = check_old_option_and_stop(grid.com, "-if", "-i"); CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL, "-i",
			       filename, PETSC_MAX_PATH_LEN, &i_set); CHKERRQ(ierr);

  if (i_set) {
    NCTool nc(&grid);
    string source;

    // Get the 'source' global attribute to check if we are given a PISM output
    // file:
    ierr = nc.open_for_reading(filename); CHKERRQ(ierr);
    ierr = nc.get_att_text(NC_GLOBAL, "source", source);
    ierr = nc.close();

    // If it's missing, print a warning
    if (source.empty()) {
      ierr = verbPrintf(1, grid.com,
			"PISM WARNING: file '%s' does not have the 'source' global attribute.\n"
			"     If '%s' is a PISM output file, please run the following to get rid of this warning:\n"
			"     ncatted -a source,global,c,c,PISM %s\n",
			filename, filename, filename); CHKERRQ(ierr);
    } else if (source.find("PISM") == string::npos) {
      // If the 'source' attribute does not contain the string "PISM", then print
      // a message and stop:
      ierr = verbPrintf(1, grid.com,
			"PISM WARNING: '%s' does not seem to be a PISM output file.\n"
			"     If it is, please make sure that the 'source' global attribute contains the string \"PISM\".\n",
			filename); CHKERRQ(ierr);
    }

    ierr = nc.get_grid(filename);   CHKERRQ(ierr);

    // These options are ignored because we're getting *all* the grid
    // parameters from a file.
    ierr = ignore_option(grid.com, "-Mx");    CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-My");    CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-Mz");    CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-Mbz");   CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-Lx");    CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-Ly");    CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-Lz");    CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-chebZ"); CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-quadZ"); CHKERRQ(ierr);
  } else {
    ierr = set_grid_defaults(); CHKERRQ(ierr);
    ierr = set_grid_from_options(); CHKERRQ(ierr);
  }

  ierr = grid.createDA(); CHKERRQ(ierr);
  
  return 0;
}

//! Sets the starting values of model state variables.
/*!
  There are two cases:
  
  1) Initializing from a PISM output file.

  2) Setting the values using command-line options only (verification and
  simplified geometry runs, for example) or from a bootstrapping file, using
  heuristics to fill in missing and 3D fields.

  This function is called after all the memory allocation is done and all the
  physical parameters are set.
 */
PetscErrorCode IceModel::model_state_setup() {
  PetscErrorCode ierr;
  PetscTruth i_set;
  char filename[PETSC_MAX_PATH_LEN];
  
  // Check if we are initializing from a PISM output file:
  ierr = PetscOptionsGetString(PETSC_NULL, "-i",
			       filename, PETSC_MAX_PATH_LEN, &i_set); CHKERRQ(ierr);

  if (i_set) {
    ierr = initFromFile(filename); CHKERRQ(ierr);
  } else {
    ierr = set_vars_from_options(); CHKERRQ(ierr);
  }

  return 0;
}

//! Sets starting values of model state variables using command-line options.
/*!
  Sets starting values of model state variables using command-line options and
  (possibly) a bootstrapping file.

  In the base class there is only one case: bootstrapping.
 */
PetscErrorCode IceModel::set_vars_from_options() {
  PetscErrorCode ierr;
  PetscTruth boot_from_set;
  char filename[PETSC_MAX_PATH_LEN];

  ierr = verbPrintf(3, grid.com,
		    "Setting initial values of model state variables...\n"); CHKERRQ(ierr);

  ierr = PetscOptionsGetString(PETSC_NULL, "-boot_from",
			       filename, PETSC_MAX_PATH_LEN, &boot_from_set); CHKERRQ(ierr);
  
  if (boot_from_set) {
    ierr = bootstrapFromFile(filename); CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(grid.com, "PISM ERROR: No input file specified.\n"); CHKERRQ(ierr);
    PetscEnd();
  }
  
  return 0;
}

//! Initialize some physical parameters.
/*!
  This is the place for all non-trivial initialization of physical parameters
  (non-trivial meaning requiring more than just setting a value of a
  parameter).

  This method is called after memory allocation but before filling any of
  IceModelVecs.

  Rationale: all the physical parameters should be initialized before setting
  up the coupling or filling model-state variables.

  In the base class we just initialize the IceType and the shelf extension.

  Also, this is the good place for setting parameters that a user should not be
  able to override using a command-line option.
 */
PetscErrorCode IceModel::init_physics() {
  PetscErrorCode ierr;

  ierr = verbPrintf(3, grid.com,
		    "Initializing IceType and shelfExtension...\n"); CHKERRQ(ierr);

  ierr = iceFactory.setFromOptions(); CHKERRQ(ierr);
  // Initialize the IceType object:
  if (ice == PETSC_NULL) {
    ierr = iceFactory.create(&ice); CHKERRQ(ierr);
    ierr = ice->setFromOptions(); CHKERRQ(ierr); // Set options specific to this particular ice type
  }

  return 0;
}

//! Miscellaneous initialization tasks plus tasks that need the fields that can come from regridding.
PetscErrorCode IceModel::misc_setup() {
  PetscErrorCode ierr;

  ierr = verbPrintf(3, grid.com, "Finishing initialization...\n"); CHKERRQ(ierr);

  ierr = init_snapshots_from_options(); CHKERRQ(ierr);
  ierr = stampHistoryCommand(); CHKERRQ(ierr);
  ierr = createViewers(); CHKERRQ(ierr);

  // consistency of geometry after initialization;
  ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);

  // allocate and setup bed deformation model
  ierr = bedDefSetup(); CHKERRQ(ierr);

  // init basal till model, possibly inverting for phi, if desired;
  //   reads options "-topg_to_phi phi_min,phi_max,phi_ocean,topg_min,topg_max"
  //   or "-surf_vel_to_phi foo.nc";
  //   initializes PlasticBasalType* basal; sets fields vtauc, vtillphi
  ierr = initBasalTillModel(); CHKERRQ(ierr);
  
  return 0;
}

//! Initializes atmosphere and ocean couplers.
PetscErrorCode IceModel::init_couplers() {
  PetscErrorCode ierr;

  // so that we can let atmosPCC, oceanPCC know about these fields in IceModel state
  info_coupler.lat = &vLatitude;
  info_coupler.lon = &vLongitude;  
  info_coupler.mask = &vMask;
  info_coupler.thk = &vH;
  info_coupler.surfelev = &vh;
  info_coupler.topg = &vbed;

  ierr = verbPrintf(3, grid.com,
		    "Initializing atmosphere and ocean couplers...\n"); CHKERRQ(ierr);

  if (atmosPCC != PETSC_NULL) {
    ierr = atmosPCC->initFromOptions(&grid); CHKERRQ(ierr);
  } else {  SETERRQ(1,"PISM ERROR: atmosPCC == PETSC_NULL");  }
 
 if (oceanPCC != PETSC_NULL) {
    if (isDrySimulation == PETSC_TRUE) {  oceanPCC->reportInitializationToStdOut = false;  }
    ierr = oceanPCC->initFromOptions(&grid); CHKERRQ(ierr);
  } else {  SETERRQ(2,"PISM ERROR: oceanPCC == PETSC_NULL");  }


  return 0;
}


//! Allocates SSA tools and work vectors.
PetscErrorCode IceModel::allocate_internal_objects() {
  PetscErrorCode ierr;

  // a global Vec is needed for things like viewers and comm to proc zero
  ierr = DACreateGlobalVector(grid.da2, &g2); CHKERRQ(ierr);

  // setup (classical) SSA tools
  const PetscInt M = 2 * grid.Mx * grid.My;
  ierr = MatCreateMPIAIJ(grid.com, PETSC_DECIDE, PETSC_DECIDE, M, M,
                         13, PETSC_NULL, 13, PETSC_NULL,
                         &SSAStiffnessMatrix); CHKERRQ(ierr);
  ierr = VecCreateMPI(grid.com, PETSC_DECIDE, M, &SSAX); CHKERRQ(ierr);
  ierr = VecDuplicate(SSAX, &SSARHS); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_SELF, M, &SSAXLocal);
  ierr = VecScatterCreate(SSAX, PETSC_NULL, SSAXLocal, PETSC_NULL,
                          &SSAScatterGlobalToLocal); CHKERRQ(ierr);
  ierr = KSPCreate(grid.com, &SSAKSP); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(SSAKSP); CHKERRQ(ierr);

  // various internal quantities
  // 2d work vectors
  for (int j = 0; j < nWork2d; j++) {
    ierr = vWork2d[j].create(grid, "a_work_vector", true); CHKERRQ(ierr);
  }

  // 3d dedicated work vectors
  ierr = Tnew3.create(grid,"temp_new",false); CHKERRQ(ierr);
  ierr = Tnew3.set_attrs("internal", "ice temperature; temporary during update",
                         "K", ""); CHKERRQ(ierr);
  ierr = taunew3.create(grid,"age_new",false); CHKERRQ(ierr);
  ierr = taunew3.set_attrs("internal", "age of ice; temporary during update",
                           "s", ""); CHKERRQ(ierr);
  ierr = Sigmastag3[0].create(grid,"Sigma_stagx",true); CHKERRQ(ierr);
  ierr = Sigmastag3[0].set_attrs("internal",
             "rate of strain heating; on staggered grid offset in X direction",
	     "J s-1 m-3", ""); CHKERRQ(ierr);
  ierr = Sigmastag3[1].create(grid,"Sigma_stagy",true); CHKERRQ(ierr);
  ierr = Sigmastag3[1].set_attrs("internal",
             "rate of strain heating; on staggered grid offset in Y direction",
	     "J s-1 m-3", ""); CHKERRQ(ierr);
  ierr = Istag3[0].create(grid,"I_stagx",true); CHKERRQ(ierr);
  ierr = Istag3[0].set_attrs("internal","","",""); CHKERRQ(ierr);
  ierr = Istag3[1].create(grid,"I_stagy",true); CHKERRQ(ierr);
  ierr = Istag3[1].set_attrs("internal","","",""); CHKERRQ(ierr);

  return 0;
}

