// Copyright (C) 2009, 2010 Ed Bueler and Constantine Khroulev
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

//! Set default values of grid parameters.
/*!
  Derived classes (IceCompModel, for example) reimplement this to change the
  grid initialization when no -i option is set.
 */
PetscErrorCode IceModel::set_grid_defaults() {
  PetscErrorCode ierr;
  bool Mx_set, My_set, Mz_set, Lz_set, boot_from_set;
  string filename;
  grid_info gi;

  // Get the bootstrapping file name:
  
  ierr = PISMOptionsString("-boot_from", "Specifies the file to bootstrap from",
			   filename, boot_from_set); CHKERRQ(ierr);

  if (!boot_from_set) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: Please specify an input file using -i or -boot_from.\n");
    CHKERRQ(ierr);
    PetscEnd();
  }

  // Use a bootstrapping file to set some grid parameters (they can be
  // overridden later, in IceModel::set_grid_from_options()).

  // Determine the grid extent from a bootstrapping file:
  PISMIO nc(&grid);
  bool x_dim_exists, y_dim_exists, t_exists;
  ierr = nc.open_for_reading(filename.c_str()); CHKERRQ(ierr);

  ierr = nc.find_dimension("x", NULL, x_dim_exists); CHKERRQ(ierr);
  ierr = nc.find_dimension("y", NULL, y_dim_exists); CHKERRQ(ierr);
  ierr = nc.find_variable("t", NULL, t_exists); CHKERRQ(ierr);
  ierr = nc.get_grid_info(gi);

  bool mapping_exists;
  ierr = nc.find_variable("mapping", NULL, mapping_exists); CHKERRQ(ierr);
  if (mapping_exists) {
    ierr = mapping.read(filename.c_str()); CHKERRQ(ierr);
    ierr = mapping.print(); CHKERRQ(ierr);
  }

  ierr = nc.close(); CHKERRQ(ierr);

  // if the horizontal dimensions are absent then we can not proceed
  if (!x_dim_exists) {
    ierr = PetscPrintf(grid.com,"bootstrapping file '%s' has no horizontal dimension 'x'\n",
		       filename.c_str());
    CHKERRQ(ierr);
    PetscEnd();
  }
  if (!y_dim_exists) {
    ierr = PetscPrintf(grid.com,"bootstrapping file '%s' has no horizontal dimension 'y'\n",
		       filename.c_str());
    CHKERRQ(ierr);
    PetscEnd();
  }

  // Set the grid center and horizontal extent:
  grid.x0 = gi.x0;
  grid.y0 = gi.y0;
  grid.Lx = gi.Lx;
  grid.Ly = gi.Ly;

  // read grid.year if no option overrides it (avoids unnecessary reporting)
  bool ys_set;
  ierr = PISMOptionsIsSet("-ys", ys_set); CHKERRQ(ierr);
  if (!ys_set) {
    if (t_exists) {
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

  // Grid dimensions should not be deduced from a bootstrapping file, so we
  // check if these options are set and stop if they are not.
  ierr = PISMOptionsIsSet("-Mx", Mx_set); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-My", My_set); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-Mz", Mz_set); CHKERRQ(ierr);
  if ( !(Mx_set && My_set && Mz_set) ) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: All of -boot_from, -Mx, -My, -Mz are required for bootstrapping.\n");
    CHKERRQ(ierr);
    PetscEnd();
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
  bool Mx_set, My_set, Mz_set, Mbz_set, Lx_set, Ly_set, Lz_set, Lbz_set,
    z_spacing_set, zb_spacing_set;
  PetscReal x_scale = grid.Lx / 1000.0,
    y_scale = grid.Ly / 1000.0,
    z_scale = grid.Lz,
    zb_scale = grid.Lbz;

  // Process the options:

  // Read -Lx and -Ly.
  ierr = PISMOptionsReal("-Ly", "Half of the grid extent in the X direction, in km",
			 y_scale,  Ly_set); CHKERRQ(ierr);
  ierr = PISMOptionsReal("-Lx", "Half of the grid extent in the Y direction, in km",
			 x_scale,  Lx_set); CHKERRQ(ierr);
  // Vertical extent (in the ice and bedrock, correspondingly):
  ierr = PISMOptionsReal("-Lz", "Grid extent in the Z (vertical) direction in the ice, in meters",
			 z_scale,  Lz_set); CHKERRQ(ierr);
  ierr = PISMOptionsReal("-Lbz", "Grid extent in the Z (vertical) direction in the bedrock, in meters",
			 zb_scale, Lbz_set); CHKERRQ(ierr);

  // Read -Mx, -My, -Mz and -Mbz.
  ierr = PISMOptionsInt("-My", "Number of grid points in the X direction",
			grid.My, My_set); CHKERRQ(ierr);
  ierr = PISMOptionsInt("-Mx", "Number of grid points in the Y direction",
			grid.Mx, Mx_set); CHKERRQ(ierr);
  ierr = PISMOptionsInt("-Mz", "Number of grid points in the Z (vertical) direction in the ice",
			grid.Mz, Mz_set); CHKERRQ(ierr);
  ierr = PISMOptionsInt("-Mbz", "Number of grid points in the Z (vertical) direction in the bedrock",
			grid.Mbz, Mbz_set); CHKERRQ(ierr);

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

  // Determine the vertical grid spacing in the bedrock:
  ierr = PISMOptionsList(grid.com, "-zb_spacing", "Vertical spacing in the bedrock thermal layer.",
			 z_spacing_choices, "quadratic", keyword, zb_spacing_set); CHKERRQ(ierr);
  if (keyword == "quadratic") {
    grid.bed_vertical_spacing = QUADRATIC;
  } else {
    grid.bed_vertical_spacing = EQUAL;
  }

  if (Mbz_set) {
    if ((grid.Mbz > 1) && !Lbz_set) {
      ierr = PetscPrintf(grid.com,
			 "PISM ERROR: Please specify bedrock layer thickness using -Lbz.\n"); CHKERRQ(ierr);
      PetscEnd();
    }
  }
  
  // Use the information obtained above:
  if (Lx_set)    grid.Lx  = x_scale * 1000.0; // convert to meters
  if (Ly_set)    grid.Ly  = y_scale * 1000.0; // convert to meters
  if (Lz_set)    grid.Lz  = z_scale;	      // in meters already
  if (Lbz_set)   grid.Lbz = zb_scale;	      // in meters already

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

  ierr = PetscOptionsBegin(grid.com, "", "Options controlling input and computational grid parameters", ""); CHKERRQ(ierr);

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

    ierr = nc.get_grid(filename.c_str());   CHKERRQ(ierr);
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
    PetscEnd();
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
      PetscEnd();
    }

    if ((grid.My / grid.Ny) < 2) {
      ierr = PetscPrintf(grid.com,
			 "PISM ERROR: Can't split %d grid points between %d processors.\n",
			 grid.My, grid.Ny);
      CHKERRQ(ierr);
      PetscEnd();
    }

    if (grid.Nx * grid.Ny != grid.size) {
      ierr = PetscPrintf(grid.com,
			 "PISM ERROR: Nx * Ny has to be equal to %d.\n",
			 grid.size);
      CHKERRQ(ierr);
      PetscEnd();
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
      PetscEnd();
    }

    if (procs_x_set && procs_y_set) {
      if (tmp_x.size() != (unsigned int)grid.Nx) {
	ierr = PetscPrintf(grid.com,
			   "PISM ERROR: -Nx has to be equal to the -procs_x size.\n");
	CHKERRQ(ierr);
	PetscEnd();
      }

      if (tmp_y.size() != (unsigned int)grid.Ny) {
	ierr = PetscPrintf(grid.com,
			   "PISM ERROR: -Ny has to be equal to the -procs_y size.\n");
	CHKERRQ(ierr);
	PetscEnd();
      }

      delete [] grid.procs_x;
      delete [] grid.procs_y;

      grid.procs_x = new PetscInt[grid.Nx];
      grid.procs_y = new PetscInt[grid.Ny];

      for (PetscInt j=0; j < grid.Nx; j++)
	grid.procs_x[j] = tmp_x[j];

      for (PetscInt j=0; j < grid.Ny; j++)
	grid.procs_y[j] = tmp_y[j];
    } else {
      grid.compute_ownership_ranges();
    }
  }

  // Process -y, -ys, -ye. We are reading these options here because couplers
  // might need to know what year it is.
  ierr = set_time_from_options(); CHKERRQ(ierr);

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
  } else {
    ierr = set_vars_from_options(); CHKERRQ(ierr);
  }

  ierr = regrid(); CHKERRQ(ierr);

  // Check consistency of geometry after initialization:
  ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);

  // Initialize a bed deformation model (if needed); this should go after
  // the regrid() call.
  if (beddef) {
    ierr = beddef->init(variables); CHKERRQ(ierr);
    last_bed_def_update = grid.year;
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
  bool boot_from_set;
  string filename;

  ierr = verbPrintf(3, grid.com,
		    "Setting initial values of model state variables...\n"); CHKERRQ(ierr);

  ierr = PISMOptionsString("-boot_from", "Specifies the file to bootstrap from",
			   filename, boot_from_set); CHKERRQ(ierr);
  
  if (boot_from_set) {
    ierr = bootstrapFromFile(filename.c_str()); CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(grid.com, "PISM ERROR: No input file specified.\n"); CHKERRQ(ierr);
    PetscEnd();
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

  In the base class IceModel we just initialize the IceFlowLaw and the
  EnthalpyConverter.
 */
PetscErrorCode IceModel::init_physics() {
  PetscErrorCode ierr;

  ierr = verbPrintf(3, grid.com,
		    "initializing IceFlowLaw and EnthalpyConverter ...\n"); CHKERRQ(ierr);

  if (EC == NULL)
    EC = new EnthalpyConverter(config);

  ierr = iceFactory.setFromOptions(); CHKERRQ(ierr);

  // Initialize the IceFlowLaw object:
  if (!config.get_flag("do_cold_ice_methods")) {
    ierr = verbPrintf(2, grid.com,
      "  setting flow law to polythermal type ...\n"); CHKERRQ(ierr);
    ierr = verbPrintf(3, grid.com,
      "      (= Glen-Paterson-Budd-Lliboutry-Duval type)\n"); CHKERRQ(ierr);
    if (ice != NULL)  delete ice;  // kill choice already made
    iceFactory.setType(ICE_GPBLD); // new flowlaw which has dependence on enthalpy
                                   //   not temperature
    iceFactory.create(&ice);
    PolyThermalGPBLDIce *gpbldi = dynamic_cast<PolyThermalGPBLDIce*>(ice);
    if (gpbldi == NULL) {
      ThermoGlenIce *tgi = dynamic_cast<ThermoGlenIce*>(ice);
      if (tgi) {
        ierr = verbPrintf(2, grid.com,
          "    [flow law was actually set to ThermoGlenIce]\n");
          CHKERRQ(ierr);
      } else {
        ierr = verbPrintf(1, grid.com,
          "PISM WARNING: flow law unclear ...\n"); CHKERRQ(ierr);
      }
    }
  } else {
    ierr = verbPrintf(2, grid.com,
      "  doing cold ice methods ...\n"); CHKERRQ(ierr);

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
  ierr = ice->printInfo(4);CHKERRQ(ierr);

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


//! Allocates work vectors (and calls more).
PetscErrorCode IceModel::allocate_internal_objects() {
  PetscErrorCode ierr;
  PetscInt WIDE_STENCIL = 2;

  // since SSA tools are part of IceModel, allocate them here
  ierr = allocateSSAobjects(); CHKERRQ(ierr);

  // SIA needs a Schoof (2003)-type smoother, allocate it here
  sia_bed_smoother = new PISMBedSmoother(grid, config, WIDE_STENCIL);

  // various internal quantities
  // 2d work vectors
  for (int j = 0; j < nWork2d; j++) {
    char namestr[30];
    snprintf(namestr, sizeof(namestr), "work_vector_%d", j);
    ierr = vWork2d[j].create(grid, namestr, true, WIDE_STENCIL); CHKERRQ(ierr);
  }

  ierr = vel_ssa_old.create(grid, "bar_ssa_old", true, WIDE_STENCIL); CHKERRQ(ierr);
  // components are ubar_ssa_old and vbar_ssa_old

  // 3d work vectors
  ierr = vWork3d.create(grid,"work_vector_3d",false); CHKERRQ(ierr);

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

  ierr = hardav.create(grid, "averaged_hardness", false); CHKERRQ(ierr);
  ierr = hardav.set_attrs("internal", "vertically-averaged ice hardness",
                          "", ""); CHKERRQ(ierr);

  return 0;
}


//! Miscellaneous initialization tasks plus tasks that need the fields that can come from regridding.
PetscErrorCode IceModel::misc_setup() {
  PetscErrorCode ierr;

  ierr = verbPrintf(3, grid.com, "Finishing initialization...\n"); CHKERRQ(ierr);

  ierr = set_output_size("-o_size", "Sets the 'size' of an output file.",
			 "medium", output_vars); CHKERRQ(ierr);

  ierr = init_snapshots(); CHKERRQ(ierr);
  ierr = init_timeseries(); CHKERRQ(ierr);
  ierr = init_extras(); CHKERRQ(ierr);
  ierr = init_viewers(); CHKERRQ(ierr);

  // by now we already know if SSA velocities in the output will be valid:
  global_attributes.set_flag("pism_ssa_velocities_are_valid",
			     config.get_flag("use_ssa_velocity"));

  // set info in bed smoother based on initial bed
  ierr = vbed.beginGhostComm(); CHKERRQ(ierr);
  ierr = vbed.endGhostComm(); CHKERRQ(ierr);
  ierr = sia_bed_smoother->preprocess_bed(vbed,
               config.get("Glen_exponent"), config.get("bed_smoother_range") );
               CHKERRQ(ierr);

  // compute (possibly corrected) cell areas:
  ierr = compute_cell_areas(); CHKERRQ(ierr);

  prof = new PISMProf(&grid);

  event_step     = prof->create("step",     "time stepping (total)");
  event_velocity = prof->create("velocity", "velocity computation");
  event_vel_inc  = prof->create("vel_inc",  "vert. velocity using incompressibility");
  event_sia      = prof->create("vel_sia",  "SIA velocity computation");
  event_ssa      = prof->create("vel_ssa",  "SSA velocity computation");
  event_energy   = prof->create("energy",   "energy balance computation");
  event_vel_com  = prof->create("vel_com",  "velocity ghost points communication");
  event_thk_com  = prof->create("thk_com",  "ice thickness ghost points communication");
  event_mass     = prof->create("mass",     "mass conservation computation");
  event_age      = prof->create("age",      "age computation");
  event_beddef   = prof->create("beddef",   "bed deformation computation");
  event_output   = prof->create("output",   "writing the output file");

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
    PetscEnd();
  }
  if (ySet && yeSet) {
    ierr = PetscPrintf(grid.com, "PISM ERROR: using -y and -ye together is not allowed. Exiting...\n"); CHKERRQ(ierr);
    PetscEnd();
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
      PetscEnd();
    }
    grid.end_year = usrEndYear;
  } else if (ySet == PETSC_TRUE) {
    grid.end_year = grid.start_year + usrRunYears;
  } else {
    grid.end_year = grid.start_year + config.get("run_length_years");
  }
  return 0;
}
