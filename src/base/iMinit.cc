// Copyright (C) 2009--2014 Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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
#include <assert.h>
#include <algorithm>

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
#include "PBLingleClark.hh"
#include "PAFactory.hh"
#include "POFactory.hh"
#include "PSFactory.hh"
#include "PISMSurface.hh"
#include "PISMOcean.hh"
#include "PISMHydrology.hh"
#include "PISMMohrCoulombYieldStress.hh"
#include "PISMConstantYieldStress.hh"
#include "bedrockThermalUnit.hh"
#include "flowlaw_factory.hh"
#include "basal_resistance.hh"
#include "pism_options.hh"
#include "PISMIcebergRemover.hh"
#include "PISMOceanKill.hh"
#include "PISMCalvingAtThickness.hh"
#include "PISMEigenCalving.hh"
#include "PISMFloatKill.hh"

//! Set default values of grid parameters.
/*!
  Derived classes (IceCompModel, for example) reimplement this to change the
  grid initialization when no -i option is set.
 */
PetscErrorCode IceModel::set_grid_defaults() {
  PetscErrorCode ierr;
  bool Mx_set, My_set, Mz_set, Lz_set, boot_file_set;
  std::string filename;
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
  std::vector<std::string> names;
  names.push_back("land_ice_thickness");
  names.push_back("bedrock_altitude");
  names.push_back("thk");
  names.push_back("topg");
  bool grid_info_found = false;
  for (unsigned int i = 0; i < names.size(); ++i) {

    ierr = nc.inq_var(names[i], grid_info_found); CHKERRQ(ierr);
    if (grid_info_found == false) {
      std::string dummy1;
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
    ierr = nc.read_attributes(mapping.get_name(), mapping); CHKERRQ(ierr);
    ierr = mapping.report_to_stdout(grid.com, 4); CHKERRQ(ierr);
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
  double x_scale = grid.Lx / 1000.0,
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
  int tmp_Mx = grid.Mx, tmp_My = grid.My, tmp_Mz = grid.Mz;
  ierr = PISMOptionsInt("-My", "Number of grid points in the X direction",
                        tmp_My, My_set); CHKERRQ(ierr);
  ierr = PISMOptionsInt("-Mx", "Number of grid points in the Y direction",
                        tmp_Mx, Mx_set); CHKERRQ(ierr);
  ierr = PISMOptionsInt("-Mz", "Number of grid points in the Z (vertical) direction in the ice",
                        tmp_Mz, Mz_set); CHKERRQ(ierr);


  if (tmp_Mx > 0 && tmp_My > 0 && tmp_Mz > 0) {
    grid.Mx = tmp_Mx;
    grid.My = tmp_My;
    grid.Mz = tmp_Mz;
  } else {
    PetscPrintf(grid.com, "PISM ERROR: -Mx %d -My %d -Mz %d is invalid"
                " (have to have a positive number of grid points).\n",
                tmp_Mx, tmp_My, tmp_Mz);
    PISMEnd();
  }

  std::vector<double> x_range, y_range;
  bool x_range_set, y_range_set;
  ierr = PISMOptionsRealArray("-x_range", "min,max x coordinate values",
                              x_range, x_range_set); CHKERRQ(ierr);
  ierr = PISMOptionsRealArray("-y_range", "min,max y coordinate values",
                              y_range, y_range_set); CHKERRQ(ierr);

  std::string keyword;
  std::set<std::string> z_spacing_choices;
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
  std::set<std::string> periodicity_choices;
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
  if (Lz_set)    grid.Lz  = z_scale;          // in meters already

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
  std::string filename;

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
    std::string source;

    // Get the 'source' global attribute to check if we are given a PISM output
    // file:
    ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
    ierr = nc.get_att_text("PISM_GLOBAL", "source", source); CHKERRQ(ierr);

    bool mapping_exists;
    ierr = nc.inq_var("mapping", mapping_exists); CHKERRQ(ierr);
    if (mapping_exists) {
      ierr = nc.read_attributes(mapping.get_name(), mapping); CHKERRQ(ierr);
      ierr = mapping.report_to_stdout(grid.com, 4); CHKERRQ(ierr);
    }

    ierr = nc.close(); CHKERRQ(ierr);

    // If it's missing, print a warning
    if (source.empty()) {
      ierr = verbPrintf(1, grid.com,
                        "PISM WARNING: file '%s' does not have the 'source' global attribute.\n"
                        "     If '%s' is a PISM output file, please run the following to get rid of this warning:\n"
                        "     ncatted -a source,global,c,c,PISM %s\n",
                        filename.c_str(), filename.c_str(), filename.c_str()); CHKERRQ(ierr);
    } else if (source.find("PISM") == std::string::npos) {
      // If the 'source' attribute does not contain the string "PISM", then print
      // a message and stop:
      ierr = verbPrintf(1, grid.com,
                        "PISM WARNING: '%s' does not seem to be a PISM output file.\n"
                        "     If it is, please make sure that the 'source' global attribute contains the string \"PISM\".\n",
                        filename.c_str()); CHKERRQ(ierr);
    }

    std::vector<std::string> names;
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
    std::vector<int> tmp_x, tmp_y;
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

      for (int j=0; j < grid.Nx; j++)
        grid.procs_x[j] = tmp_x[j];

      for (int j=0; j < grid.Ny; j++)
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
  std::string filename;

  reset_counters();

  // Initialize (or re-initialize) boundary models.
  ierr = init_couplers(); CHKERRQ(ierr);

  // Check if we are initializing from a PISM output file:
  ierr = PISMOptionsString("-i", "Specifies the PISM input file",
                           filename, i_set); CHKERRQ(ierr);

  if (i_set) {
    ierr = initFromFile(filename); CHKERRQ(ierr);

    ierr = regrid(0); CHKERRQ(ierr);
    // Check consistency of geometry after initialization:
    ierr = updateSurfaceElevationAndMask(); CHKERRQ(ierr);
  } else {
    ierr = set_vars_from_options(); CHKERRQ(ierr);
  }

  // Initialize a bed deformation model (if needed); this should go after
  // the regrid(0) call.
  if (beddef) {
    ierr = beddef->init(variables); CHKERRQ(ierr);
  }

  if (btu) {
    // update surface and ocean models so that we can get the
    // temperature at the top of the bedrock
    ierr = init_step_couplers(); CHKERRQ(ierr);

    ierr = get_bed_top_temp(bedtoptemp); CHKERRQ(ierr);
    ierr = btu->init(variables); CHKERRQ(ierr);
  }

  if (subglacial_hydrology) {
    ierr = subglacial_hydrology->init(variables); CHKERRQ(ierr);
  }

  // basal_yield_stress_model->init() needs bwat so this must happen after subglacial_hydrology->init()
  if (basal_yield_stress_model) {
    ierr = basal_yield_stress_model->init(variables); CHKERRQ(ierr);
  }

  if (climatic_mass_balance_cumulative.was_created()) {
    if (i_set) {
      ierr = verbPrintf(2, grid.com,
                        "* Trying to read cumulative climatic mass balance from '%s'...\n",
                        filename.c_str()); CHKERRQ(ierr);
      ierr = climatic_mass_balance_cumulative.regrid(filename, OPTIONAL, 0.0); CHKERRQ(ierr);
    } else {
      ierr = climatic_mass_balance_cumulative.set(0.0); CHKERRQ(ierr);
    }
  }

  if (grounded_basal_flux_2D_cumulative.was_created()) {
    if (i_set) {
      ierr = verbPrintf(2, grid.com,
                        "* Trying to read cumulative grounded basal flux from '%s'...\n",
                        filename.c_str()); CHKERRQ(ierr);
      ierr = grounded_basal_flux_2D_cumulative.regrid(filename, OPTIONAL, 0.0); CHKERRQ(ierr);
    } else {
      ierr = grounded_basal_flux_2D_cumulative.set(0.0); CHKERRQ(ierr);
    }
  }

  if (floating_basal_flux_2D_cumulative.was_created()) {
    if (i_set) {
      ierr = verbPrintf(2, grid.com,
                        "* Trying to read cumulative floating basal flux from '%s'...\n",
                        filename.c_str()); CHKERRQ(ierr);
      ierr = floating_basal_flux_2D_cumulative.regrid(filename, OPTIONAL, 0.0); CHKERRQ(ierr);
    } else {
      ierr = floating_basal_flux_2D_cumulative.set(0.0); CHKERRQ(ierr);
    }
  }

  if (nonneg_flux_2D_cumulative.was_created()) {
    if (i_set) {
      ierr = verbPrintf(2, grid.com,
                        "* Trying to read cumulative nonneg flux from '%s'...\n",
                        filename.c_str()); CHKERRQ(ierr);
      ierr = nonneg_flux_2D_cumulative.regrid(filename, OPTIONAL, 0.0); CHKERRQ(ierr);
    } else {
      ierr = nonneg_flux_2D_cumulative.set(0.0); CHKERRQ(ierr);
    }
  }

  if (i_set) {
    PIO nc(grid.com, "netcdf3", grid.get_unit_system());
    bool run_stats_exists;

    ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
    ierr = nc.inq_var("run_stats", run_stats_exists); CHKERRQ(ierr);
    if (run_stats_exists) {
      ierr = nc.read_attributes(run_stats.get_name(), run_stats); CHKERRQ(ierr);
    }
    ierr = nc.close(); CHKERRQ(ierr);

    if (run_stats.has_attribute("grounded_basal_ice_flux_cumulative"))
      grounded_basal_ice_flux_cumulative = run_stats.get_double("grounded_basal_ice_flux_cumulative");

    if (run_stats.has_attribute("nonneg_rule_flux_cumulative"))
      nonneg_rule_flux_cumulative = run_stats.get_double("nonneg_rule_flux_cumulative");

    if (run_stats.has_attribute("sub_shelf_ice_flux_cumulative"))
      sub_shelf_ice_flux_cumulative = run_stats.get_double("sub_shelf_ice_flux_cumulative");

    if (run_stats.has_attribute("surface_ice_flux_cumulative"))
      surface_ice_flux_cumulative = run_stats.get_double("surface_ice_flux_cumulative");

    if (run_stats.has_attribute("sum_divQ_SIA_cumulative"))
      sum_divQ_SIA_cumulative = run_stats.get_double("sum_divQ_SIA_cumulative");

    if (run_stats.has_attribute("sum_divQ_SSA_cumulative"))
      sum_divQ_SSA_cumulative = run_stats.get_double("sum_divQ_SSA_cumulative");

    if (run_stats.has_attribute("Href_to_H_flux_cumulative"))
      Href_to_H_flux_cumulative = run_stats.get_double("Href_to_H_flux_cumulative");

    if (run_stats.has_attribute("H_to_Href_flux_cumulative"))
      H_to_Href_flux_cumulative = run_stats.get_double("H_to_Href_flux_cumulative");

    if (run_stats.has_attribute("discharge_flux_cumulative"))
      discharge_flux_cumulative = run_stats.get_double("discharge_flux_cumulative");
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
  std::string filename;

  ierr = verbPrintf(3, grid.com,
                    "Setting initial values of model state variables...\n"); CHKERRQ(ierr);

  ierr = PISMOptionsString("-boot_file", "Specifies the file to bootstrap from",
                           filename, boot_file_set); CHKERRQ(ierr);

  if (boot_file_set) {
    ierr = bootstrapFromFile(filename); CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(grid.com, "PISM ERROR: No input file specified.\n"); CHKERRQ(ierr);
    PISMEnd();
  }

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

  std::string model = config.get_string("stress_balance_model");

  ShallowStressBalance *sliding = NULL;
  if (model == "none" || model == "sia") {
    sliding = new ZeroSliding(grid, *EC, config);
  } else if (model == "prescribed_sliding" || model == "prescribed_sliding+sia") {
    sliding = new PrescribedSliding(grid, *EC, config);
  } else if (model == "ssa" || model == "ssa+sia") {
    std::string method = config.get_string("ssa_method");

    if (method == "fem") {
      sliding = new SSAFEM(grid, *EC, config);
    } else if (method == "fd") {
      sliding = new SSAFD(grid, *EC, config);
    } else {
      SETERRQ(grid.com, 1, "invalid ssa method");
    }

  } else {
    SETERRQ(grid.com, 1, "invalid stress balance model");
  }

  SSB_Modifier *modifier = NULL;
  if (model == "none" || model == "ssa" || model == "prescribed_sliding") {
    modifier = new ConstantInColumn(grid, *EC, config);
  } else if (model == "prescribed_sliding+sia" || "ssa+sia") {
    modifier = new SIAFD(grid, *EC, config);
  } else {
    SETERRQ(grid.com, 1, "invalid stress balance model");
  }

  // ~PISMStressBalance() will de-allocate sliding and modifier.
  stress_balance = new PISMStressBalance(grid, sliding, modifier, config);

  // PISM stress balance computations are diagnostic, i.e. do not
  // have a state that changes in time.  Therefore this call can be here
  // and not in model_state_setup().  We don't need to re-initialize after
  // the "diagnostic time step".
  ierr = stress_balance->init(variables); CHKERRQ(ierr);

  if (config.get_flag("include_bmr_in_continuity")) {
    ierr = stress_balance->set_basal_melt_rate(&basal_melt_rate); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceModel::allocate_iceberg_remover() {
  PetscErrorCode ierr;

  if (iceberg_remover != NULL)
    return 0;

  if (config.get_flag("kill_icebergs")) {
    iceberg_remover = new PISMIcebergRemover(grid, config);

    if (iceberg_remover == NULL) {
      PetscPrintf(grid.com, "PISM ERROR: failed to allocate the 'iceberg remover' object.\n");
      PISMEnd();
    }

    // Iceberg Remover does not have a state, so it is OK to
    // initialize here.
    ierr = iceberg_remover->init(variables); CHKERRQ(ierr);
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
  std::string hydrology_model = config.get_string("hydrology_model");

  if (subglacial_hydrology != NULL) // indicates it has already been allocated
    return 0;
  if      (hydrology_model == "null")
    subglacial_hydrology = new PISMNullTransportHydrology(grid, config);
  else if (hydrology_model == "routing")
    subglacial_hydrology = new PISMRoutingHydrology(grid, config);
  else if (hydrology_model == "distributed")
    subglacial_hydrology = new PISMDistributedHydrology(grid, config, stress_balance);
  else {
    SETERRQ1(grid.com,1,"unknown value for configuration string 'hydrology_model':\n"
             "  has value '%s'\n", hydrology_model.c_str());
  }

  return 0;
}

//! \brief Decide which basal yield stress model to use.
PetscErrorCode IceModel::allocate_basal_yield_stress() {

  if (basal_yield_stress_model != NULL)
    return 0;

  std::string model = config.get_string("stress_balance_model");

  // only these two use the yield stress (so far):
  if (model == "ssa" || model == "ssa+sia") {
    std::string yield_stress_model = config.get_string("yield_stress_model");

    if (yield_stress_model == "constant") {
      basal_yield_stress_model = new PISMConstantYieldStress(grid, config);
    } else if (yield_stress_model == "mohr_coulomb") {
      basal_yield_stress_model = new PISMMohrCoulombYieldStress(grid, config, subglacial_hydrology);
    } else {
      PetscPrintf(grid.com, "PISM ERROR: yield stress model \"%s\" is not supported.\n",
                  yield_stress_model.c_str());
      PISMEnd();
    }
  }

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

  // FIXME: someday we will have an "energy balance" sub-model...
  if (config.get_flag("do_energy") == true) {
    if (config.get_flag("do_cold_ice_methods") == false) {
      ierr = verbPrintf(2, grid.com,
                        "* Using the enthalpy-based energy balance model...\n"); CHKERRQ(ierr);
    } else {
      ierr = verbPrintf(2, grid.com,
                        "* Using the temperature-based energy balance model...\n"); CHKERRQ(ierr);
    }
  }

  // this has to go first:
  ierr = allocate_enthalpy_converter(); CHKERRQ(ierr);

  ierr = allocate_iceberg_remover(); CHKERRQ(ierr);

  ierr = allocate_stressbalance(); CHKERRQ(ierr);

  // this has to happen *after* allocate_stressbalance()
  ierr = allocate_subglacial_hydrology(); CHKERRQ(ierr);

  // this has to happen *after* allocate_subglacial_hydrology()
  ierr = allocate_basal_yield_stress(); CHKERRQ(ierr);

  ierr = allocate_bedrock_thermal_unit(); CHKERRQ(ierr);

  ierr = allocate_bed_deformation(); CHKERRQ(ierr);

  ierr = allocate_couplers(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IceModel::allocate_couplers() {
  PetscErrorCode ierr;
  // Initialize boundary models:
  PAFactory pa(grid, config);
  PSFactory ps(grid, config);
  POFactory po(grid, config);
  PISMAtmosphereModel *atmosphere;

  ierr = PetscOptionsBegin(grid.com, "", "Options choosing PISM boundary models", ""); CHKERRQ(ierr);

  if (surface == NULL) {
    ps.create(surface);
    external_surface_model = false;

    pa.create(atmosphere);
    surface->attach_atmosphere_model(atmosphere);
  }

  if (ocean == NULL) {
    po.create(ocean);
    external_ocean_model = false;
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  return 0;
}

//! Initializes atmosphere and ocean couplers.
PetscErrorCode IceModel::init_couplers() {
  PetscErrorCode ierr;

  ierr = verbPrintf(3, grid.com,
                    "Initializing boundary models...\n"); CHKERRQ(ierr);

  assert(surface != PETSC_NULL);
  ierr = surface->init(variables); CHKERRQ(ierr);

  assert(ocean != PETSC_NULL);
  ierr = ocean->init(variables); CHKERRQ(ierr);

  return 0;
}


//! Some sub-models need fields provided by surface and ocean models
//! for initialization, so here we call update() to make sure that
//! surface and ocean models report a decent state
PetscErrorCode IceModel::init_step_couplers() {
  PetscErrorCode ierr;

  assert(surface != NULL);
  assert(ocean != NULL);

  double max_dt = 0.0;
  bool restrict_dt = false;
  const double current_time = grid.time->current();
  std::vector<double> dt_restrictions;

  // Take a one year long step if we can:
  double one_year_from_now = grid.time->increment_date(current_time, 1.0);
  dt_restrictions.push_back(one_year_from_now - current_time);

  double apcc_dt = 0.0;
  ierr = surface->max_timestep(current_time, apcc_dt, restrict_dt); CHKERRQ(ierr);
  if (restrict_dt)
    dt_restrictions.push_back(apcc_dt);

  double opcc_dt = 0.0;
  ierr = ocean->max_timestep(current_time, opcc_dt, restrict_dt); CHKERRQ(ierr);
  if (restrict_dt)
    dt_restrictions.push_back(opcc_dt);

  // find the smallest of the max. time-steps reported by boundary models:
  if (dt_restrictions.empty() == false)
    max_dt = *std::min_element(dt_restrictions.begin(), dt_restrictions.end());

  // Do not take time-steps shorter than 1 second
  if (max_dt < 1.0)
    max_dt = 1.0;

  ierr = surface->update(current_time, max_dt); CHKERRQ(ierr);
  ierr = ocean->update(current_time, max_dt); CHKERRQ(ierr);

  return 0;
}


//! Allocates work vectors.
PetscErrorCode IceModel::allocate_internal_objects() {
  PetscErrorCode ierr;
  int WIDE_STENCIL = grid.max_stencil_width;

  // various internal quantities
  // 2d work vectors
  for (int j = 0; j < nWork2d; j++) {
    char namestr[30];
    snprintf(namestr, sizeof(namestr), "work_vector_%d", j);
    ierr = vWork2d[j].create(grid, namestr, WITH_GHOSTS, WIDE_STENCIL); CHKERRQ(ierr);
  }

  // 3d work vectors
  ierr = vWork3d.create(grid,"work_vector_3d",WITHOUT_GHOSTS); CHKERRQ(ierr);
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

  // Quietly re-initialize couplers (they might have done one
  // time-step during initialization)
  {
    int user_verbosity = getVerbosityLevel();
    ierr = setVerbosityLevel(1); CHKERRQ(ierr);
    ierr = init_couplers(); CHKERRQ(ierr);
    ierr = setVerbosityLevel(user_verbosity); CHKERRQ(ierr);
  }

  ierr = init_calving(); CHKERRQ(ierr);
  ierr = init_diagnostics(); CHKERRQ(ierr);
  ierr = init_snapshots(); CHKERRQ(ierr);
  ierr = init_backups(); CHKERRQ(ierr);
  ierr = init_timeseries(); CHKERRQ(ierr);
  ierr = init_extras(); CHKERRQ(ierr);
  ierr = init_viewers(); CHKERRQ(ierr);

  // Make sure that we use the output_variable_order that works with NetCDF-4,
  // "quilt", and HDF5 parallel I/O. (For different reasons, but mainly because
  // it is faster.)
  std::string o_format = config.get_string("output_format");
  if ((o_format == "netcdf4_parallel" || o_format == "quilt" || o_format == "hdf5") &&
      config.get_string("output_variable_order") != "xyz") {
    PetscPrintf(grid.com,
                "PISM ERROR: output formats netcdf4_parallel, quilt, and hdf5 require -o_order xyz.\n");
    PISMEnd();
  }

  return 0;
}

//! \brief Initialize calving mechanisms.
PetscErrorCode IceModel::init_calving() {
  PetscErrorCode ierr;

  std::istringstream arg(config.get_string("calving_methods"));
  std::string method_name;
  std::set<std::string> methods;

    while (getline(arg, method_name, ','))
      methods.insert(method_name);

  if (methods.find("ocean_kill") != methods.end()) {

    if (ocean_kill_calving == NULL) {
      ocean_kill_calving = new PISMOceanKill(grid, config);
    }

    ierr = ocean_kill_calving->init(variables); CHKERRQ(ierr);
    methods.erase("ocean_kill");
  }

  if (methods.find("thickness_calving") != methods.end()) {

    if (thickness_threshold_calving == NULL) {
      thickness_threshold_calving = new PISMCalvingAtThickness(grid, config);
    }

    ierr = thickness_threshold_calving->init(variables); CHKERRQ(ierr);
    methods.erase("thickness_calving");
  }


  if (methods.find("eigen_calving") != methods.end()) {

    if (eigen_calving == NULL) {
      eigen_calving = new PISMEigenCalving(grid, config,
                                           stress_balance);
    }

    ierr = eigen_calving->init(variables); CHKERRQ(ierr);
    methods.erase("eigen_calving");
  }

  if (methods.find("float_kill") != methods.end()) {
    if (float_kill_calving == NULL) {
      float_kill_calving = new PISMFloatKill(grid, config);
    }

    ierr = float_kill_calving->init(variables); CHKERRQ(ierr);
    methods.erase("float_kill");
  }

  std::set<std::string>::iterator j = methods.begin();
  std::string unused;
  while (j != methods.end()) {
    unused += (*j + ",");
    ++j;
  }

  if (unused.empty() == false) {
    ierr = verbPrintf(2, grid.com,
                      "PISM ERROR: calving method(s) [%s] are unknown and are ignored.\n",
                      unused.c_str()); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode IceModel::allocate_bed_deformation() {
  PetscErrorCode ierr;
  std::string model = config.get_string("bed_deformation_model");
  std::set<std::string> choices;

  choices.insert("none");
  choices.insert("iso");
  choices.insert("lc");

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

  if ((model == "lc") && (beddef == NULL)) {
    beddef = new PBLingleClark(grid, config);
    return 0;
  }

  return 0;
}
