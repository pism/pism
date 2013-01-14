// Copyright (C) 2009--2013 Ed Bueler and Constantine Khroulev
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

#include <petscdmda.h>

#include "iceModel.hh"
#include "PIO.hh"
#include "SIAFD.hh"
#include "SSAFD.hh"
#include "SSAFEM.hh"
#include "PISMStressBalance.hh"
#include "Mask.hh"
#include "enthalpyConverter.hh"
#include "varcEnthalpyConverter.hh"
#include "PISMBedDef.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "PISMHydrology.hh"
#include "PISMMohrCoulombYieldStress.hh"
#include "PISMConstantYieldStress.hh"
#include "bedrockThermalUnit.hh"
#include "flowlaw_factory.hh"
#include "basal_resistance.hh"
#include "PISMProf.hh"
#include "LocalInterpCtx.hh"
#include "pism_options.hh"

//! Set default values of grid parameters.
/*!
  Derived classes (IceCompModel, for example) reimplement this to change the
  grid initialization when no -i option is set.
 */
PetscErrorCode IceModel::set_grid_defaults() {
  PetscErrorCode ierr;
  bool Mx_set, My_set, Mz_set, Lz_set, boot_file_set;
  string filename;
  grid_info input;

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
  PIO nc(grid, "netcdf3"); // OK to use netcdf3, we read very little data here.
  bool x_dim_exists, y_dim_exists, t_exists;
  ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);

  ierr = nc.inq_dim("x", x_dim_exists); CHKERRQ(ierr);
  ierr = nc.inq_dim("y", y_dim_exists); CHKERRQ(ierr);
  ierr = nc.inq_var(config.get_string("time_dimension_name"), t_exists); CHKERRQ(ierr);

  // Try to deduce grid information from present spatial fields. This is bad,
  // because theoretically these fields may use different grids. We need a
  // better way of specifying PISM's computational grid at bootstrapping.
  vector<string> names;
  names.push_back("land_ice_thickness");
  names.push_back("bedrock_altitude");
  names.push_back("thk");
  names.push_back("topg");
  bool grid_info_found = false;
  for (unsigned int i = 0; i < names.size(); ++i) {

    ierr = nc.inq_var(names[i], grid_info_found); CHKERRQ(ierr);
    if (grid_info_found == false) {
      string dummy1;
      bool dummy2;
      ierr = nc.inq_var("dummy", names[i], grid_info_found, dummy1, dummy2); CHKERRQ(ierr);
    }

    if (grid_info_found) {
      ierr = nc.inq_grid_info(names[i], input); CHKERRQ(ierr);
      break;
    }
  }

  if (grid_info_found == false) {
    PetscPrintf(grid.com, "ERROR: no geometry information found in '%s'.\n",
                filename.c_str());
    PISMEnd();
  }

  bool mapping_exists;
  ierr = nc.inq_var("mapping", mapping_exists); CHKERRQ(ierr);
  if (mapping_exists) {
    ierr = mapping.read(filename); CHKERRQ(ierr);
    ierr = mapping.print(); CHKERRQ(ierr);
  }

  ierr = nc.close(); CHKERRQ(ierr);

  // Set the grid center and horizontal extent:
  grid.x0 = (input.x_max + input.x_min) / 2.0;
  grid.y0 = (input.y_max + input.y_min) / 2.0;
  grid.Lx = (input.x_max - input.x_min) / 2.0;
  grid.Ly = (input.y_max - input.y_min) / 2.0;

  // read current time if no option overrides it (avoids unnecessary reporting)
  bool ys_set;
  ierr = PISMOptionsIsSet("-ys", ys_set); CHKERRQ(ierr);
  if (!ys_set) {
    if (t_exists) {
      grid.time->set_start(input.time);
      ierr = verbPrintf(2, grid.com,
  		      "  time t = %s found; setting current time\n",
                        grid.time->date().c_str()); CHKERRQ(ierr);
    }
  }

  ierr =  grid.time->init(); CHKERRQ(ierr);

  // Grid dimensions should not be deduced from a bootstrapping file, so we
  // check if these options are set and stop if they are not.
  ierr = PISMOptionsIsSet("-Mx", Mx_set); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-My", My_set); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-Mz", Mz_set); CHKERRQ(ierr);
  ierr = PISMOptionsIsSet("-Lz", Lz_set); CHKERRQ(ierr);
  if ( !(Mx_set && My_set && Mz_set && Lz_set) ) {
    ierr = PetscPrintf(grid.com,
		       "PISM ERROR: All of -boot_file, -Mx, -My, -Mz, -Lz are required for bootstrapping.\n");
    CHKERRQ(ierr);
    PISMEnd();
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
  // Vertical extent (in the ice):
  ierr = PISMOptionsReal("-Lz", "Grid extent in the Z (vertical) direction in the ice, in meters",
			 z_scale,  Lz_set); CHKERRQ(ierr);

  // Read -Mx, -My, -Mz and -Mbz.
  ierr = PISMOptionsInt("-My", "Number of grid points in the X direction",
			grid.My, My_set); CHKERRQ(ierr);
  ierr = PISMOptionsInt("-Mx", "Number of grid points in the Y direction",
			grid.Mx, Mx_set); CHKERRQ(ierr);
  ierr = PISMOptionsInt("-Mz", "Number of grid points in the Z (vertical) direction in the ice",
			grid.Mz, Mz_set); CHKERRQ(ierr);

  vector<double> x_range, y_range;
  bool x_range_set, y_range_set;
  ierr = PISMOptionsRealArray("-x_range", "min,max x coordinate values",
                              x_range, x_range_set); CHKERRQ(ierr);
  ierr = PISMOptionsRealArray("-y_range", "min,max y coordinate values",
                              y_range, y_range_set); CHKERRQ(ierr);

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

  if (x_range_set && y_range_set) {
    if (x_range.size() != 2 || y_range.size() != 2) {
      PetscPrintf(grid.com, "PISM ERROR: -x_range and/or -y_range argument is invalid.\n");
      PISMEnd();
    }

    grid.x0 = (x_range[0] + x_range[1]) / 2.0;
    grid.y0 = (y_range[0] + y_range[1]) / 2.0;
    grid.Lx = (x_range[1] - x_range[0]) / 2.0;
    grid.Ly = (y_range[1] - y_range[0]) / 2.0;
  }

  grid.check_parameters();
  ierr = grid.compute_horizontal_spacing(); CHKERRQ(ierr);
  ierr = grid.compute_vertical_levels();    CHKERRQ(ierr);

  // At this point all the fields except for da2, xs, xm, ys, ym should be
  // filled. We're ready to call grid.allocate().
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
    PIO nc(grid, "guess_mode");
    string source;

    // Get the 'source' global attribute to check if we are given a PISM output
    // file:
    ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
    ierr = nc.get_att_text("PISM_GLOBAL", "source", source); CHKERRQ(ierr);

    bool mapping_exists;
    ierr = nc.inq_var("mapping", mapping_exists); CHKERRQ(ierr);
    if (mapping_exists) {
      ierr = mapping.read(filename); CHKERRQ(ierr);
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

    ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);

    bool var_exists = false;
    for (unsigned int i = 0; i < names.size(); ++i) {
      ierr = nc.inq_var(names[i], var_exists); CHKERRQ(ierr);

      if (var_exists == true) {
        ierr = nc.inq_grid(names[i], &grid, NOT_PERIODIC); CHKERRQ(ierr);
        break;
      }
    }

    if (var_exists == false) {
      PetscPrintf(grid.com, "PISM ERROR: file %s has neither enthalpy nor temperature in it!\n",
                  filename.c_str());

      ierr = nc.close(); CHKERRQ(ierr);

      PISMEnd();
    }

    ierr = nc.close(); CHKERRQ(ierr);

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

  grid.check_parameters();

  ierr = grid.allocate(); CHKERRQ(ierr);

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
  ierr = PISMOptionsString("-i", "Specifies the PISM input file",
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
  }

  if (btu) {
    PetscReal max_dt = 0;
    bool restrict = false;
    // FIXME: this will break if a surface or an ocean model requires
    // contiguous update intervals
    ierr = surface->max_timestep(grid.time->start(), max_dt, restrict); CHKERRQ(ierr);

    if (restrict == false)
      max_dt = convert(1, "year", "seconds");

    ierr = surface->update(grid.time->start(), max_dt); CHKERRQ(ierr);

    ierr = ocean->max_timestep(grid.time->start(), max_dt, restrict); CHKERRQ(ierr);

    if (restrict == false)
      max_dt = convert(1, "year", "seconds");

    ierr = ocean->update(grid.time->start(), max_dt); CHKERRQ(ierr);
    ierr = get_bed_top_temp(bedtoptemp); CHKERRQ(ierr);
    ierr = btu->init(variables); CHKERRQ(ierr);
  }

  if (subglacial_hydrology) {
    ierr = subglacial_hydrology->init(variables); CHKERRQ(ierr);
  }

  // basal_yield_stress->init() needs bwat so this must happen after subglacial_hydrology->init()
  if (basal_yield_stress) {
    ierr = basal_yield_stress->init(variables); CHKERRQ(ierr);
  }

  if (climatic_mass_balance_cumulative.was_created()) {
    if (i_set) {
      ierr = verbPrintf(2, grid.com,
                        "* Trying to read cumulative climatic mass balance from '%s'...\n",
                        filename.c_str()); CHKERRQ(ierr);
      ierr = climatic_mass_balance_cumulative.regrid(filename, 0.0); CHKERRQ(ierr);
    } else {
      ierr = climatic_mass_balance_cumulative.set(0.0); CHKERRQ(ierr);
    }
  }

  if (ocean_kill_flux_2D_cumulative.was_created()) {
    if (i_set) {
      ierr = verbPrintf(2, grid.com,
                        "* Trying to read cumulative ocean kill flux from '%s'...\n",
                        filename.c_str()); CHKERRQ(ierr);
      ierr = ocean_kill_flux_2D_cumulative.regrid(filename, 0.0); CHKERRQ(ierr);
    } else {
      ierr = ocean_kill_flux_2D_cumulative.set(0.0); CHKERRQ(ierr);
    }
  }

  ierr = compute_cell_areas(); CHKERRQ(ierr);

  // a report on whether PISM-PIK modifications of IceModel are in use
  const bool pg   = config.get_flag("part_grid"),
    pr   = config.get_flag("part_redist"),
    ki   = config.get_flag("kill_icebergs");
  if (pg || pr || ki) {
    ierr = verbPrintf(2, grid.com,
                      "* PISM-PIK mass/geometry methods are in use:  "); CHKERRQ(ierr);

    if (pg)   { ierr = verbPrintf(2, grid.com, "part_grid,"); CHKERRQ(ierr); }
    if (pr)   { ierr = verbPrintf(2, grid.com, "part_redist,"); CHKERRQ(ierr); }
    if (ki)   { ierr = verbPrintf(2, grid.com, "kill_icebergs"); CHKERRQ(ierr); }

    ierr = verbPrintf(2, grid.com, "\n"); CHKERRQ(ierr);
  }

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

//! \brief Decide which ice flow law to use.
PetscErrorCode IceModel::allocate_flowlaw() {
  PetscErrorCode ierr;
  string sia_flow_law = config.get_string("sia_flow_law"),
    ssa_flow_law = config.get_string("ssa_flow_law");

  if (config.get_flag("do_cold_ice_methods") == false) {
    ierr = verbPrintf(2, grid.com,
                      "  setting flow law to polythermal type ...\n"); CHKERRQ(ierr);
    ierr = verbPrintf(3, grid.com,
                      "      (= Glen-Paterson-Budd-Lliboutry-Duval type)\n"); CHKERRQ(ierr);

    // new flowlaw which has dependence on enthalpy, not temperature
    sia_flow_law = "gpbld";
    ssa_flow_law = "gpbld";
  } else {
    ierr = verbPrintf(2, grid.com,
                      "  doing cold ice methods ...\n"); CHKERRQ(ierr);

    sia_flow_law = "pb";
    ssa_flow_law = "pb";
  }

  config.set_string("sia_flow_law", sia_flow_law);
  config.set_string("ssa_flow_law", ssa_flow_law);

  return 0;
}

//! \brief Decide which enthalpy converter to use.
PetscErrorCode IceModel::allocate_enthalpy_converter() {
  PetscErrorCode ierr;

  if (EC != NULL)
    return 0;

  if (config.get_flag("use_linear_in_temperature_heat_capacity"))
    EC = new varcEnthalpyConverter(config);
  else
    EC = new EnthalpyConverter(config);

  if (getVerbosityLevel() > 3) {
    PetscViewer viewer;
    ierr = PetscViewerASCIIGetStdout(PETSC_COMM_WORLD,&viewer); CHKERRQ(ierr);
    ierr = EC->viewConstants(viewer); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Decide which stress balance model to use.
PetscErrorCode IceModel::allocate_stressbalance() {
  PetscErrorCode ierr;

  if (stress_balance != NULL)
    return 0;

  bool
    use_ssa_velocity = config.get_flag("use_ssa_velocity"),
    do_blatter = config.get_flag("do_blatter"),
    do_sia = config.get_flag("do_sia");

  // If both SIA and SSA are "on", the SIA and SSA velocities are always added
  // up (there is no switch saying "do the hybrid").
  if (stress_balance == NULL) {
    if (do_blatter) {
      PetscPrintf(grid.com, "Blatter solver is disabled for now.\n");
      PISMEnd();
      // stress_balance = new BlatterStressBalance(grid, ocean, config);
    } else {
      ShallowStressBalance *my_stress_balance;
      if (use_ssa_velocity) {
        string ssa_method = config.get_string("ssa_method");
        if( ssa_method == "fd" ) {
          my_stress_balance = new SSAFD(grid, *basal, *EC, config);
        } else if(ssa_method == "fem") {
          my_stress_balance = new SSAFEM(grid, *basal, *EC, config);
        } else {
          SETERRQ(grid.com, 1,"SSA algorithm flag should be one of \"fd\" or \"fem\"");
        }
      } else {
        my_stress_balance = new SSB_Trivial(grid, *basal, *EC, config);
      }
      SSB_Modifier *my_modifier;
      if (do_sia) {
        my_modifier = new SIAFD(grid, *EC, config);
      } else {
        my_modifier = new SSBM_Trivial(grid, *EC, config);
      }
      // ~PISMStressBalance() will de-allocate my_stress_balance and modifier.
      stress_balance = new PISMStressBalance(grid, my_stress_balance,
                                             my_modifier, ocean, config);
    }

    // PISM stress balance computations are diagnostic, i.e. do not
    // have a state that changes in time.  Therefore this call can be here
    // and not in model_state_setup().  We don't need to re-initialize after
    // the "diagnostic time step".
    ierr = stress_balance->init(variables); CHKERRQ(ierr);

    if (config.get_flag("include_bmr_in_continuity")) {
      ierr = stress_balance->set_basal_melt_rate(&vbmr); CHKERRQ(ierr);
    }
  }

  return 0;
}

//! \brief Decide which bedrock thermal unit to use.
PetscErrorCode IceModel::allocate_bedrock_thermal_unit() {

  if (btu != NULL)
    return 0;

  btu = new PISMBedThermalUnit(grid, config);

  return 0;
}

//! \brief Decide which subglacial hydrology model to use.
PetscErrorCode IceModel::allocate_subglacial_hydrology() {
  string hydrology_model = config.get_string("hydrology_model");

  if (subglacial_hydrology != NULL) // indicates it has already been allocated
    return 0;
  if      (hydrology_model == "tillcan")
    subglacial_hydrology = new PISMTillCanHydrology(grid, config, false);
  else if (hydrology_model == "diffuseonly")
    subglacial_hydrology = new PISMDiffuseOnlyHydrology(grid, config);
  else if (hydrology_model == "lakes")
    subglacial_hydrology = new PISMLakesHydrology(grid, config);
  else if (hydrology_model == "distributed")
    subglacial_hydrology = new PISMDistributedHydrology(grid, config, stress_balance);
  else { SETERRQ(grid.com,1,"unknown value for 'hydrology_model'"); }

  return 0;
}

//! \brief Decide which basal yield stress model to use.
PetscErrorCode IceModel::allocate_basal_yield_stress() {
  PetscErrorCode ierr;

  if (basal_yield_stress != NULL)
    return 0;

  bool use_ssa_velocity = config.get_flag("use_ssa_velocity"),
    do_blatter = config.get_flag("do_blatter");

  if (use_ssa_velocity || do_blatter) {
    bool hold_tauc;
    ierr = PISMOptionsIsSet("-hold_tauc", hold_tauc); CHKERRQ(ierr);

    if (hold_tauc) {
      basal_yield_stress = new PISMConstantYieldStress(grid, config);
    } else {
      basal_yield_stress = new PISMMohrCoulombYieldStress(grid, config, subglacial_hydrology);
    }
  }

  return 0;
}

//! \brief Decide which basal resistance law to use.
PetscErrorCode IceModel::allocate_basal_resistance_law() {

  if (basal != NULL)
    return 0;

  if (config.get_flag("do_pseudo_plastic_till") == true)
    basal = new IceBasalResistancePseudoPlasticLaw(config);
  else
    basal = new IceBasalResistancePlasticLaw(config);

  return 0;
}

//! Allocate PISM's sub-models implementing some physical processes.
/*!
  This method is called after memory allocation but before filling any of
  IceModelVecs because all the physical parameters should be initialized before
  setting up the coupling or filling model-state variables.
 */
PetscErrorCode IceModel::allocate_submodels() {
  PetscErrorCode ierr;

  // this has to go first:
  ierr = allocate_enthalpy_converter(); CHKERRQ(ierr);
  // then this:
  ierr = allocate_flowlaw(); CHKERRQ(ierr);

  // this has to happen before allocate_stressbalance() is called
  ierr = allocate_basal_resistance_law(); CHKERRQ(ierr);

  ierr = allocate_stressbalance(); CHKERRQ(ierr);

  // this has to happen *after* allocate_stressbalance()
  ierr = allocate_subglacial_hydrology(); CHKERRQ(ierr);

  // this has to happen *after* allocate_subglacial_hydrology()
  ierr = allocate_basal_yield_stress(); CHKERRQ(ierr);

  ierr = allocate_bedrock_thermal_unit(); CHKERRQ(ierr);

  ierr = allocate_bed_deformation(); CHKERRQ(ierr);

  return 0;
}


//! Initializes atmosphere and ocean couplers.
PetscErrorCode IceModel::init_couplers() {
  PetscErrorCode ierr;

  ierr = verbPrintf(3, grid.com,
		    "Initializing boundary models...\n"); CHKERRQ(ierr);

  if (surface != PETSC_NULL) {
    ierr = surface->init(variables); CHKERRQ(ierr);
  } else {  SETERRQ(grid.com, 2,"PISM ERROR: surface == PETSC_NULL");  }

  if (ocean != PETSC_NULL) {
    ierr = ocean->init(variables); CHKERRQ(ierr);
  } else {  SETERRQ(grid.com, 2,"PISM ERROR: ocean == PETSC_NULL");  }

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

  ierr = init_ocean_kill(); CHKERRQ(ierr);
  ierr = init_diagnostics(); CHKERRQ(ierr);
  ierr = init_snapshots(); CHKERRQ(ierr);
  ierr = init_backups(); CHKERRQ(ierr);
  ierr = init_timeseries(); CHKERRQ(ierr);
  ierr = init_extras(); CHKERRQ(ierr);
  ierr = init_viewers(); CHKERRQ(ierr);

  // Make sure that we use the output_variable_order that works with NetCDF-4,
  // "quilt", and HDF5 parallel I/O. (For different reasons, but mainly because
  // it is faster.)
  string o_format = config.get_string("output_format");
  if ((o_format == "netcdf4_parallel" || o_format == "quilt" || o_format == "hdf5") &&
      config.get_string("output_variable_order") != "xyz") {
    PetscPrintf(grid.com,
                "PISM ERROR: output formats netcdf4_parallel, quilt, and hdf5 require -o_order xyz.\n");
    PISMEnd();
  }

  event_step      = grid.profiler->create("step",     "time spent doing time-stepping");
  event_velocity  = grid.profiler->create("velocity", "time spent updating ice velocity");

  event_energy  = grid.profiler->create("energy",   "time spent inside energy time-stepping");
  event_hydrology = grid.profiler->create("hydrology",   "time spent inside hydrology time-stepping");
  event_age     = grid.profiler->create("age",      "time spent inside age time-stepping");
  event_mass    = grid.profiler->create("masscont", "time spent inside mass continuity time-stepping");

  event_beddef  = grid.profiler->create("bed_def",  "time spent updating the bed deformation model");

  event_output    = grid.profiler->create("output", "time spent writing output files");
  event_output_define = grid.profiler->create("output_define", "time spent defining variables");
  event_snapshots = grid.profiler->create("snapshots", "time spent writing snapshots");
  event_backups   = grid.profiler->create("backups", "time spent writing backups");

  return 0;
}

//! \brief Initialize the mask used by the -ocean_kill code.
PetscErrorCode IceModel::init_ocean_kill() {
  PetscErrorCode ierr;
  string filename;
  bool flag;

  if (!config.get_flag("ocean_kill"))
    return 0;

  ierr = PetscOptionsBegin(grid.com, "", "Fixed calving front options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-ocean_kill", "Specifies a file to get -ocean_kill thickness from",
                             filename, flag, true); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  MaskQuery m(vMask);

  IceModelVec2S thickness, *tmp;

  if (filename.empty()) {
    ierr = verbPrintf(2, grid.com,
       "* Option -ocean_kill seen: using ice thickness at the beginning of the run\n"
       "  to set the fixed calving front location.\n"); CHKERRQ(ierr);
    tmp = &vH;
  } else {
    ierr = verbPrintf(2, grid.com,
       "* Option -ocean_kill seen: setting fixed calving front location using\n"
       "  ice thickness from '%s'.\n",filename.c_str()); CHKERRQ(ierr);

    ierr = thickness.create(grid, "thk", false); CHKERRQ(ierr);
    ierr = thickness.set_attrs("temporary", "land ice thickness",
                               "m", "land_ice_thickness"); CHKERRQ(ierr);
    ierr = thickness.set_attr("valid_min", 0.0); CHKERRQ(ierr);

    ierr = thickness.regrid(filename, true); CHKERRQ(ierr);

    tmp = &thickness;
  }

  ierr = ocean_kill_mask.begin_access(); CHKERRQ(ierr);
  ierr = tmp->begin_access(); CHKERRQ(ierr);
  ierr = vMask.begin_access(); CHKERRQ(ierr);

  for (PetscInt   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (PetscInt j = grid.ys; j < grid.ys+grid.ym; ++j) {
      if ((*tmp)(i, j) > 0 || m.grounded(i, j) ) // FIXME: use GeometryCalculator
        ocean_kill_mask(i, j) = 0;
      else
        ocean_kill_mask(i, j) = 1;
    }
  }

  ierr = vMask.end_access(); CHKERRQ(ierr);
  ierr = tmp->end_access(); CHKERRQ(ierr);
  ierr = ocean_kill_mask.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceModel::allocate_bed_deformation() {
  PetscErrorCode ierr;
  string model = config.get_string("bed_deformation_model");
  set<string> choices;

  ierr = check_old_option_and_stop(grid.com, "-bed_def_iso", "-bed_def"); CHKERRQ(ierr);
  ierr = check_old_option_and_stop(grid.com, "-bed_def_lc",  "-bed_def"); CHKERRQ(ierr);

  choices.insert("none");
  choices.insert("iso");
#if (PISM_USE_FFTW==1)
  choices.insert("lc");
#endif

  ierr = PetscOptionsBegin(grid.com, "", "Bed deformation model", ""); CHKERRQ(ierr);
  {
    bool dummy;
    ierr = PISMOptionsList(grid.com, "-bed_def", "Specifies a bed deformation model.",
			 choices, model, model, dummy); CHKERRQ(ierr);

  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (model == "none")
    return 0;

  if ((model == "iso") && (beddef == NULL)) {
    beddef = new PBPointwiseIsostasy(grid, config);
    return 0;
  }

#if (PISM_USE_FFTW==1)
  if ((model == "lc") && (beddef == NULL)) {
    beddef = new PBLingleClark(grid, config);
    return 0;
  }
#endif

  return 0;
}
