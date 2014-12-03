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
#include "error_handling.hh"

namespace pism {

//! Set default values of grid parameters.
/*!
  Derived classes (IceCompModel, for example) reimplement this to change the
  grid initialization when no -i option is set.
 */
void IceModel::set_grid_defaults() {
  bool Mx_set, My_set, Mz_set, Lz_set, boot_file_set;
  std::string filename;
  grid_info input;

  // Get the bootstrapping file name:

  OptionsString("-boot_file", "Specifies the file to bootstrap from",
                filename, boot_file_set);

  if (!boot_file_set) {
    throw RuntimeError("Please specify an input file using -i or -boot_file.");
  }

  // Use a bootstrapping file to set some grid parameters (they can be
  // overridden later, in IceModel::set_grid_from_options()).

  // Determine the grid extent from a bootstrapping file:
  PIO nc(grid, "netcdf3"); // OK to use netcdf3, we read very little data here.
  bool x_dim_exists, y_dim_exists, t_exists;
  nc.open(filename, PISM_READONLY);

  x_dim_exists = nc.inq_dim("x");
  y_dim_exists = nc.inq_dim("y");
  t_exists = nc.inq_var(config.get_string("time_dimension_name"));

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

    grid_info_found = nc.inq_var(names[i]);
    if (grid_info_found == false) {
      std::string dummy1;
      bool dummy2;
      nc.inq_var("dummy", names[i], grid_info_found, dummy1, dummy2);
    }

    if (grid_info_found) {
      input = nc.inq_grid_info(names[i], grid.periodicity);
      break;
    }
  }

  if (grid_info_found == false) {
    throw RuntimeError::formatted("no geometry information found in '%s'",
                                  filename.c_str());
  }

  std::string proj4_string = nc.get_att_text("PISM_GLOBAL", "proj4");
  if (proj4_string.empty() == false) {
    global_attributes.set_string("proj4", proj4_string);
  }

  bool mapping_exists = nc.inq_var("mapping");
  if (mapping_exists) {
    nc.read_attributes(mapping.get_name(), mapping);
    mapping.report_to_stdout(grid.com, 4);
  }

  nc.close();

  // Set the grid center and horizontal extent:
  grid.x0 = input.x0;
  grid.y0 = input.y0;
  grid.Lx = input.Lx;
  grid.Ly = input.Ly;

  // read current time if no option overrides it (avoids unnecessary reporting)
  bool ys_set;
  OptionsIsSet("-ys", ys_set);
  if (!ys_set) {
    if (t_exists) {
      grid.time->set_start(input.time);
      verbPrintf(2, grid.com,
                 "  time t = %s found; setting current time\n",
                 grid.time->date().c_str());
    }
  }

  grid.time->init();

  // Grid dimensions should not be deduced from a bootstrapping file, so we
  // check if these options are set and stop if they are not.
  OptionsIsSet("-Mx", Mx_set);
  OptionsIsSet("-My", My_set);
  OptionsIsSet("-Mz", Mz_set);
  OptionsIsSet("-Lz", Lz_set);
  if (not (Mx_set && My_set && Mz_set && Lz_set)) {
    throw RuntimeError("All of -boot_file, -Mx, -My, -Mz, -Lz are required for bootstrapping.");
  }
}

//! Initalizes the grid from options.
/*! Reads all of -Mx, -My, -Mz, -Mbz, -Lx, -Ly, -Lz, -Lbz, -z_spacing and
    -zb_spacing. Sets corresponding grid parameters.
 */
void IceModel::set_grid_from_options() {
  bool Mx_set, My_set, Mz_set, Lx_set, Ly_set, Lz_set,
    z_spacing_set;
  double x_scale = grid.Lx / 1000.0,
    y_scale = grid.Ly / 1000.0,
    z_scale = grid.Lz;

  // Process the options:

  // Read -Lx and -Ly.
  OptionsReal("-Ly", "Half of the grid extent in the X direction, in km",
              y_scale,  Ly_set);
  OptionsReal("-Lx", "Half of the grid extent in the Y direction, in km",
              x_scale,  Lx_set);
  // Vertical extent (in the ice):
  OptionsReal("-Lz", "Grid extent in the Z (vertical) direction in the ice, in meters",
              z_scale,  Lz_set);

  // Read -Mx, -My, -Mz and -Mbz.
  int tmp_Mx = grid.Mx(), tmp_My = grid.My(), tmp_Mz = grid.Mz;
  OptionsInt("-My", "Number of grid points in the X direction",
             tmp_My, My_set);
  OptionsInt("-Mx", "Number of grid points in the Y direction",
             tmp_Mx, Mx_set);
  OptionsInt("-Mz", "Number of grid points in the Z (vertical) direction in the ice",
             tmp_Mz, Mz_set);


  if (tmp_Mx > 0 && tmp_My > 0 && tmp_Mz > 0) {
    grid.set_Mx(tmp_Mx);
    grid.set_My(tmp_My);
    grid.Mz = tmp_Mz;
  } else {
    throw RuntimeError::formatted("-Mx %d -My %d -Mz %d is invalid\n"
                                  "(have to have a positive number of grid points).",
                                  tmp_Mx, tmp_My, tmp_Mz);
  }

  std::vector<double> x_range, y_range;
  bool x_range_set, y_range_set;
  OptionsRealArray("-x_range", "min,max x coordinate values",
                   x_range, x_range_set);
  OptionsRealArray("-y_range", "min,max y coordinate values",
                   y_range, y_range_set);

  std::string keyword;
  std::set<std::string> z_spacing_choices;
  z_spacing_choices.insert("quadratic");
  z_spacing_choices.insert("equal");
  // Determine the vertical grid spacing in the ice:
  OptionsList("-z_spacing", "Vertical spacing in the ice.",
              z_spacing_choices, "quadratic", keyword, z_spacing_set);

  if (keyword == "quadratic") {
    grid.ice_vertical_spacing = QUADRATIC;
  } else {
    grid.ice_vertical_spacing = EQUAL;
  }

  // Use the information obtained above:
  if (Lx_set) {
    grid.Lx  = x_scale * 1000.0; // convert to meters
  }
  if (Ly_set) {
    grid.Ly  = y_scale * 1000.0; // convert to meters
  }
  if (Lz_set) {
    grid.Lz  = z_scale;          // in meters already
  }

  if (x_range_set && y_range_set) {
    if (x_range.size() != 2 || y_range.size() != 2) {
      throw RuntimeError("-x_range and/or -y_range argument is invalid.");
    }

    grid.x0 = (x_range[0] + x_range[1]) / 2.0;
    grid.y0 = (y_range[0] + y_range[1]) / 2.0;
    grid.Lx = (x_range[1] - x_range[0]) / 2.0;
    grid.Ly = (y_range[1] - y_range[0]) / 2.0;
  }

  grid.check_parameters();
  grid.compute_horizontal_spacing();
  grid.compute_vertical_levels();

  // At this point all the fields except for da2, xs, xm, ys, ym should be
  // filled. We're ready to call grid.allocate().
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
void IceModel::grid_setup() {
  bool i_set;
  std::string filename;

  verbPrintf(3, grid.com,
             "Setting up the computational grid...\n");

  // Check if we are initializing from a PISM output file:
  OptionsString("-i", "Specifies a PISM input file",
                filename, i_set);

  if (i_set) {
    PIO nc(grid, "guess_mode");

    // Get the 'source' global attribute to check if we are given a PISM output
    // file:
    nc.open(filename, PISM_READONLY);
    std::string source = nc.get_att_text("PISM_GLOBAL", "source");

    std::string proj4_string = nc.get_att_text("PISM_GLOBAL", "proj4");
    if (proj4_string.empty() == false) {
      global_attributes.set_string("proj4", proj4_string);
    }

    bool mapping_exists = nc.inq_var("mapping");
    if (mapping_exists) {
      nc.read_attributes(mapping.get_name(), mapping);
      mapping.report_to_stdout(grid.com, 4);
    }

    nc.close();

    // If it's missing, print a warning
    if (source.empty()) {
      verbPrintf(1, grid.com,
                 "PISM WARNING: file '%s' does not have the 'source' global attribute.\n"
                 "     If '%s' is a PISM output file, please run the following to get rid of this warning:\n"
                 "     ncatted -a source,global,c,c,PISM %s\n",
                 filename.c_str(), filename.c_str(), filename.c_str());
    } else if (source.find("PISM") == std::string::npos) {
      // If the 'source' attribute does not contain the string "PISM", then print
      // a message and stop:
      verbPrintf(1, grid.com,
                 "PISM WARNING: '%s' does not seem to be a PISM output file.\n"
                 "     If it is, please make sure that the 'source' global attribute contains the string \"PISM\".\n",
                 filename.c_str());
    }

    std::vector<std::string> names;
    names.push_back("enthalpy");
    names.push_back("temp");

    nc.open(filename, PISM_READONLY);

    bool var_exists = false;
    for (unsigned int i = 0; i < names.size(); ++i) {
      var_exists = nc.inq_var(names[i]);

      if (var_exists == true) {
        nc.inq_grid(names[i], &grid, grid.periodicity);
        break;
      }
    }

    if (var_exists == false) {
      nc.close();
      throw RuntimeError::formatted("file %s has neither enthalpy nor temperature in it",
                                    filename.c_str());
    }

    nc.close();

    // These options are ignored because we're getting *all* the grid
    // parameters from a file.
    ignore_option(grid.com, "-Mx");
    ignore_option(grid.com, "-My");
    ignore_option(grid.com, "-Mz");
    ignore_option(grid.com, "-Mbz");
    ignore_option(grid.com, "-Lx");
    ignore_option(grid.com, "-Ly");
    ignore_option(grid.com, "-Lz");
    ignore_option(grid.com, "-z_spacing");
  } else {
    set_grid_defaults();
    set_grid_from_options();
  }

  grid.allocate();
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
void IceModel::model_state_setup() {
  bool i_set;
  std::string filename;

  reset_counters();

  // Initialize (or re-initialize) boundary models.
  init_couplers();

  // Check if we are initializing from a PISM output file:
  OptionsString("-i", "Specifies the PISM input file",
                filename, i_set);

  if (i_set) {
    initFromFile(filename);

    regrid(0);
    // Check consistency of geometry after initialization:
    updateSurfaceElevationAndMask();
  } else {
    set_vars_from_options();
  }

  // Initialize a bed deformation model (if needed); this should go after
  // the regrid(0) call.
  if (beddef) {
    beddef->init(variables);
  }

  if (btu) {
    bool bootstrapping_needed = false;
    btu->init(variables, bootstrapping_needed);

    if (bootstrapping_needed == true) {
      // update surface and ocean models so that we can get the
      // temperature at the top of the bedrock
      verbPrintf(2, grid.com,
                 "getting surface B.C. from couplers...\n");
      init_step_couplers();

      get_bed_top_temp(bedtoptemp);

      btu->bootstrap();
    }
  }

  if (subglacial_hydrology) {
    subglacial_hydrology->init(variables);
  }

  // basal_yield_stress_model->init() needs bwat so this must happen after subglacial_hydrology->init()
  if (basal_yield_stress_model) {
    basal_yield_stress_model->init(variables);
  }

  if (climatic_mass_balance_cumulative.was_created()) {
    if (i_set) {
      verbPrintf(2, grid.com,
                 "* Trying to read cumulative climatic mass balance from '%s'...\n",
                 filename.c_str());
      climatic_mass_balance_cumulative.regrid(filename, OPTIONAL, 0.0);
    } else {
      climatic_mass_balance_cumulative.set(0.0);
    }
  }

  if (grounded_basal_flux_2D_cumulative.was_created()) {
    if (i_set) {
      verbPrintf(2, grid.com,
                 "* Trying to read cumulative grounded basal flux from '%s'...\n",
                 filename.c_str());
      grounded_basal_flux_2D_cumulative.regrid(filename, OPTIONAL, 0.0);
    } else {
      grounded_basal_flux_2D_cumulative.set(0.0);
    }
  }

  if (floating_basal_flux_2D_cumulative.was_created()) {
    if (i_set) {
      verbPrintf(2, grid.com,
                 "* Trying to read cumulative floating basal flux from '%s'...\n",
                 filename.c_str());
      floating_basal_flux_2D_cumulative.regrid(filename, OPTIONAL, 0.0);
    } else {
      floating_basal_flux_2D_cumulative.set(0.0);
    }
  }

  if (nonneg_flux_2D_cumulative.was_created()) {
    if (i_set) {
      verbPrintf(2, grid.com,
                 "* Trying to read cumulative nonneg flux from '%s'...\n",
                 filename.c_str());
      nonneg_flux_2D_cumulative.regrid(filename, OPTIONAL, 0.0);
    } else {
      nonneg_flux_2D_cumulative.set(0.0);
    }
  }

  if (i_set) {
    PIO nc(grid.com, "netcdf3", grid.get_unit_system());

    nc.open(filename, PISM_READONLY);
    bool run_stats_exists = nc.inq_var("run_stats");
    if (run_stats_exists) {
      nc.read_attributes(run_stats.get_name(), run_stats);
    }
    nc.close();

    if (run_stats.has_attribute("grounded_basal_ice_flux_cumulative")) {
      grounded_basal_ice_flux_cumulative = run_stats.get_double("grounded_basal_ice_flux_cumulative");
    }

    if (run_stats.has_attribute("nonneg_rule_flux_cumulative")) {
      nonneg_rule_flux_cumulative = run_stats.get_double("nonneg_rule_flux_cumulative");
    }

    if (run_stats.has_attribute("sub_shelf_ice_flux_cumulative")) {
      sub_shelf_ice_flux_cumulative = run_stats.get_double("sub_shelf_ice_flux_cumulative");
    }

    if (run_stats.has_attribute("surface_ice_flux_cumulative")) {
      surface_ice_flux_cumulative = run_stats.get_double("surface_ice_flux_cumulative");
    }

    if (run_stats.has_attribute("sum_divQ_SIA_cumulative")) {
      sum_divQ_SIA_cumulative = run_stats.get_double("sum_divQ_SIA_cumulative");
    }

    if (run_stats.has_attribute("sum_divQ_SSA_cumulative")) {
      sum_divQ_SSA_cumulative = run_stats.get_double("sum_divQ_SSA_cumulative");
    }

    if (run_stats.has_attribute("Href_to_H_flux_cumulative")) {
      Href_to_H_flux_cumulative = run_stats.get_double("Href_to_H_flux_cumulative");
    }

    if (run_stats.has_attribute("H_to_Href_flux_cumulative")) {
      H_to_Href_flux_cumulative = run_stats.get_double("H_to_Href_flux_cumulative");
    }

    if (run_stats.has_attribute("discharge_flux_cumulative")) {
      discharge_flux_cumulative = run_stats.get_double("discharge_flux_cumulative");
    }
  }

  compute_cell_areas();

  // a report on whether PISM-PIK modifications of IceModel are in use
  const bool pg   = config.get_flag("part_grid"),
    pr   = config.get_flag("part_redist"),
    ki   = config.get_flag("kill_icebergs");
  if (pg || pr || ki) {
    verbPrintf(2, grid.com,
               "* PISM-PIK mass/geometry methods are in use:  ");

    if (pg)   {
      verbPrintf(2, grid.com, "part_grid,");
    }
    if (pr)   {
      verbPrintf(2, grid.com, "part_redist,");
    }
    if (ki)   {
      verbPrintf(2, grid.com, "kill_icebergs");
    }

    verbPrintf(2, grid.com, "\n");
  }

  stampHistoryCommand();
}

//! Sets starting values of model state variables using command-line options.
/*!
  Sets starting values of model state variables using command-line options and
  (possibly) a bootstrapping file.

  In the base class there is only one case: bootstrapping.
 */
void IceModel::set_vars_from_options() {
  bool boot_file_set;
  std::string filename;

  verbPrintf(3, grid.com,
             "Setting initial values of model state variables...\n");

  OptionsString("-boot_file", "Specifies the file to bootstrap from",
                filename, boot_file_set);

  if (boot_file_set) {
    bootstrapFromFile(filename);
  } else {
    throw RuntimeError("No input file specified.");
  }
}

//! \brief Decide which enthalpy converter to use.
void IceModel::allocate_enthalpy_converter() {

  if (EC != NULL) {
    return;
  }

  if (config.get_flag("use_linear_in_temperature_heat_capacity")) {
    EC = new varcEnthalpyConverter(config);
  } else {
    EC = new EnthalpyConverter(config);
  }
}

//! \brief Decide which stress balance model to use.
void IceModel::allocate_stressbalance() {

  if (stress_balance != NULL) {
    return;
  }

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
      throw RuntimeError::formatted("invalid ssa method: %s", method.c_str());
    }

  } else {
    throw RuntimeError::formatted("invalid stress balance model: %s", model.c_str());
  }

  SSB_Modifier *modifier = NULL;
  if (model == "none" || model == "ssa" || model == "prescribed_sliding") {
    modifier = new ConstantInColumn(grid, *EC, config);
  } else if (model == "prescribed_sliding+sia" || "ssa+sia") {
    modifier = new SIAFD(grid, *EC, config);
  } else {
    throw RuntimeError::formatted("invalid stress balance model: %s", model.c_str());
  }

  // ~StressBalance() will de-allocate sliding and modifier.
  stress_balance = new StressBalance(grid, sliding, modifier, config);

  // PISM stress balance computations are diagnostic, i.e. do not
  // have a state that changes in time.  Therefore this call can be here
  // and not in model_state_setup().  We don't need to re-initialize after
  // the "diagnostic time step".
  stress_balance->init(variables);

  if (config.get_flag("include_bmr_in_continuity")) {
    stress_balance->set_basal_melt_rate(&basal_melt_rate);
  }
}

void IceModel::allocate_iceberg_remover() {

  if (iceberg_remover != NULL) {
    return;
  }

  if (config.get_flag("kill_icebergs")) {

    // this will throw an exception on failure
    iceberg_remover = new IcebergRemover(grid, config);

    // Iceberg Remover does not have a state, so it is OK to
    // initialize here.
    iceberg_remover->init(variables);
  }
}

//! \brief Decide which bedrock thermal unit to use.
void IceModel::allocate_bedrock_thermal_unit() {

  if (btu != NULL) {
    return;
  }

  btu = new BedThermalUnit(grid, config);
}

//! \brief Decide which subglacial hydrology model to use.
void IceModel::allocate_subglacial_hydrology() {
  std::string hydrology_model = config.get_string("hydrology_model");

  if (subglacial_hydrology != NULL) { // indicates it has already been allocated
    return;
  }

  if (hydrology_model == "null") {
    subglacial_hydrology = new NullTransportHydrology(grid, config);
  } else if (hydrology_model == "routing") {
    subglacial_hydrology = new RoutingHydrology(grid, config);
  } else if (hydrology_model == "distributed") {
    subglacial_hydrology = new DistributedHydrology(grid, config, stress_balance);
  } else {
    throw RuntimeError::formatted("unknown value for configuration string 'hydrology_model':\n"
                                  "has value '%s'", hydrology_model.c_str());
  }
}

//! \brief Decide which basal yield stress model to use.
void IceModel::allocate_basal_yield_stress() {

  if (basal_yield_stress_model != NULL) {
    return;
  }

  std::string model = config.get_string("stress_balance_model");

  // only these two use the yield stress (so far):
  if (model == "ssa" || model == "ssa+sia") {
    std::string yield_stress_model = config.get_string("yield_stress_model");

    if (yield_stress_model == "constant") {
      basal_yield_stress_model = new ConstantYieldStress(grid, config);
    } else if (yield_stress_model == "mohr_coulomb") {
      basal_yield_stress_model = new MohrCoulombYieldStress(grid, config, subglacial_hydrology);
    } else {
      throw RuntimeError::formatted("yield stress model '%s' is not supported.",
                                    yield_stress_model.c_str());
    }
  }
}

//! Allocate PISM's sub-models implementing some physical processes.
/*!
  This method is called after memory allocation but before filling any of
  IceModelVecs because all the physical parameters should be initialized before
  setting up the coupling or filling model-state variables.
 */
void IceModel::allocate_submodels() {

  // FIXME: someday we will have an "energy balance" sub-model...
  if (config.get_flag("do_energy") == true) {
    if (config.get_flag("do_cold_ice_methods") == false) {
      verbPrintf(2, grid.com,
                 "* Using the enthalpy-based energy balance model...\n");
    } else {
      verbPrintf(2, grid.com,
                 "* Using the temperature-based energy balance model...\n");
    }
  }

  // this has to go first:
  allocate_enthalpy_converter();

  allocate_iceberg_remover();

  allocate_stressbalance();

  // this has to happen *after* allocate_stressbalance()
  allocate_subglacial_hydrology();

  // this has to happen *after* allocate_subglacial_hydrology()
  allocate_basal_yield_stress();

  allocate_bedrock_thermal_unit();

  allocate_bed_deformation();

  allocate_couplers();
}


void IceModel::allocate_couplers() {
  // Initialize boundary models:
  PAFactory pa(grid, config);
  PSFactory ps(grid, config);
  POFactory po(grid, config);
  AtmosphereModel *atmosphere;

  if (surface == NULL) {
    surface = ps.create();
    external_surface_model = false;

    atmosphere = pa.create();
    surface->attach_atmosphere_model(atmosphere);
  }

  if (ocean == NULL) {
    ocean = po.create();
    external_ocean_model = false;
  }
}

//! Initializes atmosphere and ocean couplers.
void IceModel::init_couplers() {

  verbPrintf(3, grid.com,
             "Initializing boundary models...\n");

  assert(surface != NULL);
  surface->init(variables);

  assert(ocean != NULL);
  ocean->init(variables);
}


//! Some sub-models need fields provided by surface and ocean models
//! for initialization, so here we call update() to make sure that
//! surface and ocean models report a decent state
void IceModel::init_step_couplers() {

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
  surface->max_timestep(current_time, apcc_dt, restrict_dt);
  if (restrict_dt) {
    dt_restrictions.push_back(apcc_dt);
  }

  double opcc_dt = 0.0;
  ocean->max_timestep(current_time, opcc_dt, restrict_dt);
  if (restrict_dt) {
    dt_restrictions.push_back(opcc_dt);
  }

  // find the smallest of the max. time-steps reported by boundary models:
  if (dt_restrictions.empty() == false) {
    max_dt = *std::min_element(dt_restrictions.begin(), dt_restrictions.end());
  }

  // Do not take time-steps shorter than 1 second
  if (max_dt < 1.0) {
    max_dt = 1.0;
  }

  surface->update(current_time, max_dt);
  ocean->update(current_time, max_dt);
}


//! Allocates work vectors.
void IceModel::allocate_internal_objects() {
  const unsigned int WIDE_STENCIL = config.get("grid_max_stencil_width");

  // various internal quantities
  // 2d work vectors
  for (int j = 0; j < nWork2d; j++) {
    char namestr[30];
    snprintf(namestr, sizeof(namestr), "work_vector_%d", j);
    vWork2d[j].create(grid, namestr, WITH_GHOSTS, WIDE_STENCIL);
  }

  // 3d work vectors
  vWork3d.create(grid,"work_vector_3d",WITHOUT_GHOSTS);
  vWork3d.set_attrs("internal",
                    "e.g. new values of temperature or age or enthalpy during time step",
                    "", "");
}


//! Miscellaneous initialization tasks plus tasks that need the fields that can come from regridding.
void IceModel::misc_setup() {

  verbPrintf(3, grid.com, "Finishing initialization...\n");

  output_size_from_option("-o_size", "Sets the 'size' of an output file.",
                          "medium", output_vars);

  // Quietly re-initialize couplers (they might have done one
  // time-step during initialization)
  {
    int user_verbosity = getVerbosityLevel();
    setVerbosityLevel(1);
    init_couplers();
    setVerbosityLevel(user_verbosity);
  }

  init_calving();
  init_diagnostics();
  init_snapshots();
  init_backups();
  init_timeseries();
  init_extras();
  init_viewers();

  // Make sure that we use the output_variable_order that works with NetCDF-4,
  // "quilt", and HDF5 parallel I/O. (For different reasons, but mainly because
  // it is faster.)
  std::string o_format = config.get_string("output_format");
  if ((o_format == "netcdf4_parallel" || o_format == "quilt" || o_format == "hdf5") &&
      config.get_string("output_variable_order") != "xyz") {
    throw RuntimeError("output formats netcdf4_parallel, quilt, and hdf5 require -o_order xyz.");
  }
}

//! \brief Initialize calving mechanisms.
void IceModel::init_calving() {

  std::istringstream arg(config.get_string("calving_methods"));
  std::string method_name;
  std::set<std::string> methods;

    while (getline(arg, method_name, ',')) {
      methods.insert(method_name);
    }

  if (methods.find("ocean_kill") != methods.end()) {

    if (ocean_kill_calving == NULL) {
      ocean_kill_calving = new OceanKill(grid, config);
    }

    ocean_kill_calving->init(variables);
    methods.erase("ocean_kill");
  }

  if (methods.find("thickness_calving") != methods.end()) {

    if (thickness_threshold_calving == NULL) {
      thickness_threshold_calving = new CalvingAtThickness(grid, config);
    }

    thickness_threshold_calving->init(variables);
    methods.erase("thickness_calving");
  }


  if (methods.find("eigen_calving") != methods.end()) {

    if (eigen_calving == NULL) {
      eigen_calving = new EigenCalving(grid, config,
                                           stress_balance);
    }

    eigen_calving->init(variables);
    methods.erase("eigen_calving");
  }

  if (methods.find("float_kill") != methods.end()) {
    if (float_kill_calving == NULL) {
      float_kill_calving = new FloatKill(grid, config);
    }

    float_kill_calving->init(variables);
    methods.erase("float_kill");
  }

  std::set<std::string>::iterator j = methods.begin();
  std::string unused;
  while (j != methods.end()) {
    unused += (*j + ",");
    ++j;
  }

  if (unused.empty() == false) {
    verbPrintf(2, grid.com,
               "PISM ERROR: calving method(s) [%s] are unknown and are ignored.\n",
               unused.c_str());
  }
}

void IceModel::allocate_bed_deformation() {
  std::string model = config.get_string("bed_deformation_model");

  if (beddef == NULL) {
    if (model == "none") {
      beddef = NULL;
      return;
    }

    if (model == "iso") {
      beddef = new PBPointwiseIsostasy(grid, config);
      return;
    }

    if (model == "lc") {
      beddef = new PBLingleClark(grid, config);
      return;
    }
  }
}

} // end of namespace pism
