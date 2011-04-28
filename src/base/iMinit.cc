// Copyright (C) 2009--2011 Ed Bueler and Constantine Khroulev
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
#include "PISMIO.hh"
#include "SIAFD.hh"
#include "SSAFD.hh"
#include "SSAFEM.hh"

//! Set default values of grid parameters.
/*!
  Derived classes (IceCompModel, for example) reimplement this to change the
  grid initialization when no -i option is set.
 */
PetscErrorCode IceModel::set_grid_defaults() {
  PetscErrorCode ierr;
  bool Mx_set, My_set, Mz_set, Lz_set, boot_file_set;
  string filename;
  grid_info gi;

  // Get the bootstrapping file name:
  
  ierr = PISMOptionsString("-boot_file", "Specifies the file to bootstrap from",
			   filename, boot_file_set); CHKERRQ(ierr);

  if (!boot_file_set) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: Please specify an input file using -i or -boot_file.\n");
    CHKERRQ(ierr);
    PISMEnd();
  }

  // Use a bootstrapping file to set some grid parameters (they can be
  // overridden later, in IceModel::set_grid_from_options()).

  // Determine the grid extent from a bootstrapping file:
  PISMIO nc(&grid);
  bool x_dim_exists, y_dim_exists, t_exists, time_exists;
  ierr = nc.open_for_reading(filename.c_str()); CHKERRQ(ierr);

  ierr = nc.find_dimension("x", NULL, x_dim_exists); CHKERRQ(ierr);
  ierr = nc.find_dimension("y", NULL, y_dim_exists); CHKERRQ(ierr);
  ierr = nc.find_variable("t", NULL, t_exists); CHKERRQ(ierr);
  ierr = nc.find_variable("time", NULL, time_exists); CHKERRQ(ierr);

  // Try to deduce grid information from present spatial fields. This is bad,
  // because theoretically these fields may use different grids. We need a
  // better way of specifying PISM's computational grid at bootstrapping.
  vector<string> names;
  names.push_back("land_ice_thickness");
  names.push_back("bedrock_altitude");
  names.push_back("thk");
  names.push_back("topg");
  for (unsigned int i = 0; i < names.size(); ++i) {
    ierr = nc.get_grid_info(names[i], gi);
    if (ierr == 0) break;
  }

  if (ierr != 0) {
    PetscPrintf(grid.com, "ERROR: no geometry information found in '%s'.\n",
                filename.c_str());
    PISMEnd();
  }

  bool mapping_exists;
  ierr = nc.find_variable("mapping", NULL, mapping_exists); CHKERRQ(ierr);
  if (mapping_exists) {
    ierr = mapping.read(filename.c_str()); CHKERRQ(ierr);
    ierr = mapping.print(); CHKERRQ(ierr);
  }

  ierr = nc.close(); CHKERRQ(ierr);

  // Set the grid center and horizontal extent:
  grid.x0 = gi.x0;
  grid.y0 = gi.y0;
  grid.Lx = gi.Lx;
  grid.Ly = gi.Ly;

  // read grid.year if no option overrides it (avoids unnecessary reporting)
  bool ys_set;
  ierr = PISMOptionsIsSet("-ys", ys_set); CHKERRQ(ierr);
  if (!ys_set) {
    if (t_exists || time_exists) {
      grid.year = gi.time / secpera; // set year from read-in time variable
      ierr = verbPrintf(2, grid.com, 
  		      "  time t = %5.4f years found; setting current year\n",
		      grid.year); CHKERRQ(ierr);
    } else {
      grid.year = 0.0;
      ierr = verbPrintf(2, grid.com, 
		      "  time dimension was not found; setting current year to 0.0 years\n",
		      grid.year); CHKERRQ(ierr);
    }
  }
  grid.start_year = grid.year;

  // Grid dimensions should not be deduced from a bootstrapping file, so we
  // check if these options are set and stop if they are not.
  ierr = PISMOptionsIsSet("-Mx", Mx_set); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-My", My_set); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-Mz", Mz_set); CHKERRQ(ierr);
  if ( !(Mx_set && My_set && Mz_set) ) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: All of -boot_file, -Mx, -My, -Mz are required for bootstrapping.\n");
    CHKERRQ(ierr);
    PISMEnd();
  }

  ierr = PISMOptionsIsSet("-Lz", Lz_set); CHKERRQ(ierr);
  if ( !Lz_set ) {
    ierr = verbPrintf(2, grid.com,
		      "PISM WARNING: -Lz is not set; trying to deduce it using the bootstrapping file...\n");
    CHKERRQ(ierr);
  }

  return 0;
}

//! Initalizes the grid from options.
/*! Reads all of -Mx, -My, -Mz, -Mbz, -Lx, -Ly, -Lz, -Lbz, -z_spacing and
    -zb_spacing. Sets corresponding grid parameters.
 */
PetscErrorCode IceModel::set_grid_from_options() {
  PetscErrorCode ierr;
  bool Mx_set, My_set, Mz_set, Lx_set, Ly_set, Lz_set,
    z_spacing_set, periodicity_set;
  PetscReal x_scale = grid.Lx / 1000.0,
    y_scale = grid.Ly / 1000.0,
    z_scale = grid.Lz;

  // Process the options:

  // Read -Lx and -Ly.
  ierr = PISMOptionsReal("-Ly", "Half of the grid extent in the X direction, in km",
			 y_scale,  Ly_set); CHKERRQ(ierr);
  ierr = PISMOptionsReal("-Lx", "Half of the grid extent in the Y direction, in km",
			 x_scale,  Lx_set); CHKERRQ(ierr);
  // Vertical extent (in the ice and bedrock, correspondingly):
  ierr = PISMOptionsReal("-Lz", "Grid extent in the Z (vertical) direction in the ice, in meters",
			 z_scale,  Lz_set); CHKERRQ(ierr);

  // Read -Mx, -My, -Mz and -Mbz.
  ierr = PISMOptionsInt("-My", "Number of grid points in the X direction",
			grid.My, My_set); CHKERRQ(ierr);
  ierr = PISMOptionsInt("-Mx", "Number of grid points in the Y direction",
			grid.Mx, Mx_set); CHKERRQ(ierr);
  ierr = PISMOptionsInt("-Mz", "Number of grid points in the Z (vertical) direction in the ice",
			grid.Mz, Mz_set); CHKERRQ(ierr);

  string keyword;
  set<string> z_spacing_choices;
  z_spacing_choices.insert("quadratic");
  z_spacing_choices.insert("equal");
  // Determine the vertical grid spacing in the ice:
  ierr = PISMOptionsList(grid.com, "-z_spacing", "Vertical spacing in the ice.",
			 z_spacing_choices, "quadratic", keyword, z_spacing_set); CHKERRQ(ierr);
			 
  if (keyword == "quadratic") {
    grid.ice_vertical_spacing = QUADRATIC;
  } else {
    grid.ice_vertical_spacing = EQUAL;
  }

  // Determine grid periodicity:
  set<string> periodicity_choices;
  periodicity_choices.insert("none");
  periodicity_choices.insert("x");
  periodicity_choices.insert("y");
  periodicity_choices.insert("xy");
  ierr = PISMOptionsList(grid.com, "-periodicity", "Horizontal grid periodicity.",
			 periodicity_choices, "none", keyword, periodicity_set); CHKERRQ(ierr);
  if (periodicity_set) {
    if (keyword == "none")
      grid.periodicity = NONE;
    else if (keyword == "x")
      grid.periodicity = X_PERIODIC;
    else if (keyword == "y")
      grid.periodicity = Y_PERIODIC;
    else if (keyword == "xy")
      grid.periodicity = XY_PERIODIC;
  }

  // Use the information obtained above:
  if (Lx_set)    grid.Lx  = x_scale * 1000.0; // convert to meters
  if (Ly_set)    grid.Ly  = y_scale * 1000.0; // convert to meters
  if (Lz_set)    grid.Lz  = z_scale;	      // in meters already

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
  influencing the grid (currently: -Mx, -My, -Mz, -Mbz, -Lx, -Ly, -Lz, -z_spacing,
  -zb_spacing) are ignored.

  2) Initializing using defaults, command-line options and (possibly) a
  bootstrapping file. Derived classes requiring special grid setup should
  reimplement IceGrid::set_grid_from_options().

  No memory allocation should happen here.
 */
PetscErrorCode IceModel::grid_setup() {
  PetscErrorCode ierr;
  bool i_set;
  string filename;

  ierr = PetscOptionsBegin(grid.com, "",
                           "Options controlling input and computational grid parameters",
                           ""); CHKERRQ(ierr);

  ierr = verbPrintf(3, grid.com,
		    "Setting up the computational grid...\n"); CHKERRQ(ierr);

  // Check if we are initializing from a PISM output file:
  ierr = PISMOptionsString("-i", "Specifies a PISM input file",
			   filename, i_set); CHKERRQ(ierr);

  if (i_set) {
    PISMIO nc(&grid);
    string source;

    // Get the 'source' global attribute to check if we are given a PISM output
    // file:
    ierr = nc.open_for_reading(filename.c_str()); CHKERRQ(ierr);
    ierr = nc.get_att_text(NC_GLOBAL, "source", source); CHKERRQ(ierr);

    bool mapping_exists;
    ierr = nc.find_variable("mapping", NULL, mapping_exists); CHKERRQ(ierr);
    if (mapping_exists) {
      ierr = mapping.read(filename.c_str()); CHKERRQ(ierr);
      ierr = mapping.print(); CHKERRQ(ierr);
    }

    ierr = nc.close(); CHKERRQ(ierr);

    // If it's missing, print a warning
    if (source.empty()) {
      ierr = verbPrintf(1, grid.com,
			"PISM WARNING: file '%s' does not have the 'source' global attribute.\n"
			"     If '%s' is a PISM output file, please run the following to get rid of this warning:\n"
			"     ncatted -a source,global,c,c,PISM %s\n",
			filename.c_str(), filename.c_str(), filename.c_str()); CHKERRQ(ierr);
    } else if (source.find("PISM") == string::npos) {
      // If the 'source' attribute does not contain the string "PISM", then print
      // a message and stop:
      ierr = verbPrintf(1, grid.com,
			"PISM WARNING: '%s' does not seem to be a PISM output file.\n"
			"     If it is, please make sure that the 'source' global attribute contains the string \"PISM\".\n",
			filename.c_str()); CHKERRQ(ierr);
    }

    vector<string> names;
    names.push_back("enthalpy");
    names.push_back("temp");
    for (unsigned int i = 0; i < names.size(); ++i) {
      ierr = nc.get_grid(filename, names[i]);
      if (ierr == 0) break;
    }

    if (ierr != 0) {
      PetscPrintf(grid.com, "PISM ERROR: file %s has neither enthalpy nor temperature in it!\n",
                  filename.c_str()); CHKERRQ(ierr);
      PISMEnd();
    }

    grid.start_year = grid.year; // can be overridden using the -ys option

    // These options are ignored because we're getting *all* the grid
    // parameters from a file.
    ierr = ignore_option(grid.com, "-Mx");    CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-My");    CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-Mz");    CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-Mbz");   CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-Lx");    CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-Ly");    CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-Lz");    CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-z_spacing"); CHKERRQ(ierr);
    ierr = ignore_option(grid.com, "-zb_spacing"); CHKERRQ(ierr);
  } else {
    ierr = set_grid_defaults(); CHKERRQ(ierr);
    ierr = set_grid_from_options(); CHKERRQ(ierr);
  }

  bool Nx_set, Ny_set;
  ierr = PISMOptionsInt("-Nx", "Number of processors in the x direction",
			grid.Nx, Nx_set); CHKERRQ(ierr);
  ierr = PISMOptionsInt("-Ny", "Number of processors in the y direction",
			grid.Ny, Ny_set); CHKERRQ(ierr);

  if (Nx_set ^ Ny_set) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: Please set both -Nx and -Ny.\n");
    CHKERRQ(ierr);
    PISMEnd();
  }

  if ((!Nx_set) && (!Ny_set)) {
    grid.compute_nprocs();
    grid.compute_ownership_ranges();
  } else {

    if ((grid.Mx / grid.Nx) < 2) {
      ierr = PetscPrintf(grid.com,
			 "PISM ERROR: Can't split %d grid points between %d processors.\n",
			 grid.Mx, grid.Nx);
      CHKERRQ(ierr);
      PISMEnd();
    }

    if ((grid.My / grid.Ny) < 2) {
      ierr = PetscPrintf(grid.com,
			 "PISM ERROR: Can't split %d grid points between %d processors.\n",
			 grid.My, grid.Ny);
      CHKERRQ(ierr);
      PISMEnd();
    }

    if (grid.Nx * grid.Ny != grid.size) {
      ierr = PetscPrintf(grid.com,
			 "PISM ERROR: Nx * Ny has to be equal to %d.\n",
			 grid.size);
      CHKERRQ(ierr);
      PISMEnd();
    }

    bool procs_x_set, procs_y_set;
    vector<PetscInt> tmp_x, tmp_y;
    ierr = PISMOptionsIntArray("-procs_x", "Processor ownership ranges (x direction)",
			       tmp_x, procs_x_set); CHKERRQ(ierr);
    ierr = PISMOptionsIntArray("-procs_y", "Processor ownership ranges (y direction)",
			       tmp_y, procs_y_set); CHKERRQ(ierr);

    if (procs_x_set ^ procs_y_set) {
      ierr = PetscPrintf(grid.com,
			 "PISM ERROR: Please set both -procs_x and -procs_y.\n");
      CHKERRQ(ierr);
      PISMEnd();
    }

    if (procs_x_set && procs_y_set) {
      if (tmp_x.size() != (unsigned int)grid.Nx) {
	ierr = PetscPrintf(grid.com,
			   "PISM ERROR: -Nx has to be equal to the -procs_x size.\n");
	CHKERRQ(ierr);
	PISMEnd();
      }

      if (tmp_y.size() != (unsigned int)grid.Ny) {
	ierr = PetscPrintf(grid.com,
			   "PISM ERROR: -Ny has to be equal to the -procs_y size.\n");
	CHKERRQ(ierr);
	PISMEnd();
      }

      grid.procs_x.resize(grid.Nx);
      grid.procs_y.resize(grid.Ny);

      for (PetscInt j=0; j < grid.Nx; j++)
	grid.procs_x[j] = tmp_x[j];

      for (PetscInt j=0; j < grid.Ny; j++)
	grid.procs_y[j] = tmp_y[j];
    } else {
      grid.compute_ownership_ranges();
    }
  } // -Nx and -Ny set

  // Process -y, -ys, -ye. We are reading these options here because couplers
  // might need to know what year it is.
  ierr = set_time_from_options(); CHKERRQ(ierr);

  grid.check_parameters();

  ierr = grid.createDA(); CHKERRQ(ierr);

  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  return 0;
}

//! Sets the starting values of model state variables.
/*!
  There are two cases:
  
  1) Initializing from a PISM output file.

  2) Setting the values using command-line options only (verification and
  simplified geometry runs, for example) or from a bootstrapping file, using
  heuristics to fill in missing and 3D fields.

  Calls IceModel::regrid().

  This function is called after all the memory allocation is done and all the
  physical parameters are set.

  Calling this method should be all one needs to set model state variables.
  Please avoid modifying them in other parts of the initialization sequence.

  Also, please avoid operations that would make it unsafe to call this more
  than once (memory allocation is one example).
 */
PetscErrorCode IceModel::model_state_setup() {
  PetscErrorCode ierr;
  bool i_set;
  string filename;
  
  // Check if we are initializing from a PISM output file:
  ierr = PISMOptionsString("-i", "Specifies a PISM input file",
			   filename, i_set); CHKERRQ(ierr);

  if (i_set) {
    ierr = initFromFile(filename.c_str()); CHKERRQ(ierr);

    ierr = regrid(0); CHKERRQ(ierr);
    // Check consistency of geometry after initialization:
    ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);
  } else {
    ierr = set_vars_from_options(); CHKERRQ(ierr);
  }

  // Initialize a bed deformation model (if needed); this should go after
  // the regrid() call.
  if (beddef) {
    ierr = beddef->init(variables); CHKERRQ(ierr);
    last_bed_def_update = grid.year;
  }

  if (btu) {
    ierr = get_bed_top_temp(bedtoptemp); CHKERRQ(ierr);
    ierr = btu->init(variables); CHKERRQ(ierr);
  }

  // init basal till model, possibly inverting for phi, if desired;
  //   reads options "-topg_to_phi phi_min,phi_max,phi_ocean,topg_min,topg_max"
  //   or "-surf_vel_to_phi foo.nc";
  //   initializes IceBasalResistancePlasticLaw* basal; sets fields vtauc, vtillphi
  ierr = initBasalTillModel(); CHKERRQ(ierr);

  ierr = stampHistoryCommand(); CHKERRQ(ierr);

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
  bool boot_file_set;
  string filename;

  ierr = verbPrintf(3, grid.com,
		    "Setting initial values of model state variables...\n"); CHKERRQ(ierr);

  ierr = PISMOptionsString("-boot_file", "Specifies the file to bootstrap from",
			   filename, boot_file_set); CHKERRQ(ierr);
  
  if (boot_file_set) {
    ierr = bootstrapFromFile(filename.c_str()); CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(grid.com, "PISM ERROR: No input file specified.\n"); CHKERRQ(ierr);
    PISMEnd();
  }
  
  return 0;
}

//! Initialize some physical parameters.
/*!
  This is the place for all non-trivial initialization of physical parameters.
  ("Non-trivial" means that the initialization requires more than just setting
  a value of a parameter.  Such trivial changes can go here too, or earlier.)
  Also, this is the good place to set those parameters that a user should not be
  able to override using a command-line option.

  This method is called after memory allocation but before filling any of
  IceModelVecs because all the physical parameters should be initialized before
  setting up the coupling or filling model-state variables.
 */
PetscErrorCode IceModel::init_physics() {
  PetscErrorCode ierr;

  if (ice == NULL) {
    // Initialize the IceFlowLaw object:
    if (!config.get_flag("do_cold_ice_methods")) {
      ierr = verbPrintf(2, grid.com,
                        "  setting flow law to polythermal type ...\n"); CHKERRQ(ierr);
      ierr = verbPrintf(3, grid.com,
                        "      (= Glen-Paterson-Budd-Lliboutry-Duval type)\n"); CHKERRQ(ierr);

      // new flowlaw which has dependence on enthalpy, not temperature
      ice = new GPBLDIce(grid.com, "", config);
    } else {
      ierr = verbPrintf(2, grid.com,
                        "  doing cold ice methods ...\n"); CHKERRQ(ierr);

      ierr = iceFactory.setFromOptions(); CHKERRQ(ierr);

      // FIXME:  the semantics of IceFlowLaw should be cleared up; lots of PISM
      //   (e.g. verification and EISMINT II and EISMINT-Greenland) are cold,
      //   but the really important cases (e.g. SeaRISE-Greenland) are polythermal
      // in cold case we may have various IceFlowLaw s, e.g. set by derived classes
      if (ice == PETSC_NULL) {
        ierr = iceFactory.create(&ice); CHKERRQ(ierr);
      }
    }

    // set options specific to this particular ice type:
    ierr = ice->setFromOptions(); CHKERRQ(ierr);
  }

  // Create the stress balance object:
  PetscScalar pseudo_plastic_q = config.get("pseudo_plastic_q"),
    pseudo_plastic_uthreshold = config.get("pseudo_plastic_uthreshold") / secpera,
    plastic_regularization = config.get("plastic_regularization") / secpera;

  bool do_pseudo_plastic_till = config.get_flag("do_pseudo_plastic_till"),
    use_ssa_velocity = config.get_flag("use_ssa_velocity"),
    do_sia = config.get_flag("do_sia");
  
  if (basal == NULL)
    basal = new IceBasalResistancePlasticLaw(plastic_regularization, do_pseudo_plastic_till, 
                                             pseudo_plastic_q, pseudo_plastic_uthreshold);

  if (EC == NULL) {
    EC = new EnthalpyConverter(config);
    if (getVerbosityLevel() > 3) {
      PetscViewer viewer;
      ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer); CHKERRQ(ierr);
      ierr = EC->viewConstants(viewer); CHKERRQ(ierr);
    }
  }

  if (btu == NULL) {
    btu = new PISMBedThermalUnit(grid, config);
  }

  // If both SIA and SSA are "on", the SIA and SSA velocities are always added
  // up (there is no switch saying "do the hybrid").
  if (stress_balance == NULL) {
    ShallowStressBalance *my_stress_balance;

    SSB_Modifier *modifier;
    if (do_sia) {
      modifier = new SIAFD(grid, *ice, *EC, config);
    } else {
      modifier = new SSBM_Trivial(grid, *ice, *EC, config);
    }

    if (use_ssa_velocity) {
      string ssa_method = config.get_string("ssa_method");
      if( ssa_method == "fd" ) {
        my_stress_balance = new SSAFD(grid, *basal, *ice, *EC, config);
      } else if(ssa_method == "fem") {
        my_stress_balance = new SSAFEM(grid, *basal, *ice, *EC, config);
      } else {
        SETERRQ(1,"SSA algorithm flag should be one of \"fd\" or \"fem\"");
      }
    } else {
      my_stress_balance = new SSB_Trivial(grid, *basal, *ice, *EC, config);
    }
  
    // ~PISMStressBalance() will de-allocate my_stress_balance and modifier.
    stress_balance = new PISMStressBalance(grid, my_stress_balance,
                                           modifier, ocean, config);

    // Note that in PISM stress balance computations are diagnostic, i.e. do not
    // have a state that changes in time. This means that this call can be here
    // and not in model_state_setup() and we don't need to re-initialize after
    // the "diagnostic time step".
    ierr = stress_balance->init(variables); CHKERRQ(ierr);

    if (config.get_flag("include_bmr_in_continuity")) {
      ierr = stress_balance->set_basal_melt_rate(&vbmr); CHKERRQ(ierr);
    }
  }

  ierr = bed_def_setup(); CHKERRQ(ierr);

  return 0;
}


//! Initializes atmosphere and ocean couplers.
PetscErrorCode IceModel::init_couplers() {
  PetscErrorCode ierr;

  ierr = verbPrintf(3, grid.com,
		    "Initializing boundary models...\n"); CHKERRQ(ierr);

  if (surface != PETSC_NULL) {
    ierr = surface->init(variables); CHKERRQ(ierr);
  } else {  SETERRQ(2,"PISM ERROR: surface == PETSC_NULL");  }

  if (ocean != PETSC_NULL) {
    ierr = ocean->init(variables); CHKERRQ(ierr);
  } else {  SETERRQ(2,"PISM ERROR: ocean == PETSC_NULL");  }

  return 0;
}


//! Allocates work vectors.
PetscErrorCode IceModel::allocate_internal_objects() {
  PetscErrorCode ierr;
  PetscInt WIDE_STENCIL = grid.max_stencil_width;

  // various internal quantities
  // 2d work vectors
  for (int j = 0; j < nWork2d; j++) {
    char namestr[30];
    snprintf(namestr, sizeof(namestr), "work_vector_%d", j);
    ierr = vWork2d[j].create(grid, namestr, true, WIDE_STENCIL); CHKERRQ(ierr);
  }

  ierr = vWork2dV.create(grid, "vWork2dV", true); CHKERRQ(ierr);
  ierr = vWork2dV.set_attrs("internal", "velocity work vector", "", ""); CHKERRQ(ierr);

  // 3d work vectors
  ierr = vWork3d.create(grid,"work_vector_3d",false); CHKERRQ(ierr);
  ierr = vWork3d.set_attrs(
           "internal", 
           "e.g. new values of temperature or age or enthalpy during time step",
           "", ""); CHKERRQ(ierr);

  return 0;
}


//! Miscellaneous initialization tasks plus tasks that need the fields that can come from regridding.
PetscErrorCode IceModel::misc_setup() {
  PetscErrorCode ierr;

  ierr = verbPrintf(3, grid.com, "Finishing initialization...\n"); CHKERRQ(ierr);

  ierr = set_output_size("-o_size", "Sets the 'size' of an output file.",
			 "medium", output_vars); CHKERRQ(ierr);

  ierr = init_diagnostics(); CHKERRQ(ierr); 
  ierr = init_snapshots(); CHKERRQ(ierr);
  ierr = init_backups(); CHKERRQ(ierr);
  ierr = init_timeseries(); CHKERRQ(ierr);
  ierr = init_extras(); CHKERRQ(ierr);
  ierr = init_viewers(); CHKERRQ(ierr);

  // compute (possibly corrected) cell areas:
  ierr = compute_cell_areas(); CHKERRQ(ierr);

  event_step      = grid.profiler->create("step",     "time spent doing time-stepping");
  event_velocity  = grid.profiler->create("velocity", "time spent updating ice velocity");

  event_energy  = grid.profiler->create("energy",   "time spent inside energy time-stepping");
  event_age     = grid.profiler->create("age",      "time spent inside age time-stepping");
  event_mass    = grid.profiler->create("masscont", "time spent inside mass continuity time-stepping");

  event_beddef  = grid.profiler->create("bed_def",  "time spent updating the bed deformation model");

  event_output    = grid.profiler->create("output", "time spent writing an output file");
  event_snapshots = grid.profiler->create("snapshots", "time spent writing snapshots");
  event_backups   = grid.profiler->create("backups", "time spent writing backups");

  return 0;
}

//! Determine the run length, starting and ending years using command-line options.
PetscErrorCode  IceModel::set_time_from_options() {
  PetscErrorCode ierr;

  // read options about year of start, year of end, number of run years;
  // note grid.year has already been set from input file or defaults
  PetscReal usrStartYear = grid.start_year,
    usrEndYear = grid.end_year,
    usrRunYears = grid.end_year - grid.start_year;
  bool ysSet = false, yeSet = false, ySet = false;

  ierr = PetscOptionsBegin(grid.com, "", "Time: start year, end year and run length", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsReal("-ys", "Start year",        usrStartYear, ysSet); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-ye", "End year",          usrEndYear,   yeSet); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-y",  "years; Run length", usrRunYears,  ySet);  CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);


  if (ysSet && yeSet && ySet) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: all of -y, -ys, -ye are set. Exiting...\n");
    CHKERRQ(ierr);
    PISMEnd();
  }
  if (ySet && yeSet) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: using -y and -ye together is not allowed. Exiting...\n"); CHKERRQ(ierr);
    PISMEnd();
  }

  // Set the start year if -ys is set, use the default (stored in
  // grid.start_year) otherwise.
  if (ysSet == PETSC_TRUE) {
    grid.start_year = usrStartYear;
    grid.year = usrStartYear;
  } else {
    grid.year = grid.start_year;
  }

  if (yeSet == PETSC_TRUE) {
    if (usrEndYear < grid.start_year) {
      ierr = PetscPrintf(grid.com,
			"PISM ERROR: -ye (%3.3f) is less than -ys (%3.3f) (or input file year or default).\n"
			"PISM cannot run backward in time.\n",
			 usrEndYear, grid.start_year); CHKERRQ(ierr);
      PISMEnd();
    }
    grid.end_year = usrEndYear;
  } else if (ySet == PETSC_TRUE) {
    grid.end_year = grid.start_year + usrRunYears;
  } else {
    grid.end_year = grid.start_year + config.get("run_length_years");
  }

  t_years_TempAge = grid.year;

  return 0;
}
