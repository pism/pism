// Copyright (C) 2009--2018 Ed Bueler and Constantine Khroulev
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

#include "IceModel.hh"
#include "pism/basalstrength/ConstantYieldStress.hh"
#include "pism/basalstrength/MohrCoulombYieldStress.hh"
#include "pism/basalstrength/basal_resistance.hh"
#include "pism/calving/CalvingAtThickness.hh"
#include "pism/calving/EigenCalving.hh"
#include "pism/calving/vonMisesCalving.hh"
#include "pism/calving/FloatKill.hh"
#include "pism/calving/IcebergRemover.hh"
#include "pism/calving/OceanKill.hh"
#include "pism/calving/FrontalMelt.hh"
#include "pism/energy/BedThermalUnit.hh"
#include "pism/hydrology/NullTransport.hh"
#include "pism/hydrology/Routing.hh"
#include "pism/hydrology/Distributed.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/stressbalance/sia/SIAFD.hh"
#include "pism/stressbalance/ssa/SSAFD.hh"
#include "pism/stressbalance/ssa/SSAFEM.hh"
#include "pism/util/Mask.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/pism_options.hh"
#include "pism/coupler/OceanModel.hh"
#include "pism/coupler/SurfaceModel.hh"
#include "pism/coupler/atmosphere/Factory.hh"
#include "pism/coupler/ocean/Factory.hh"
#include "pism/coupler/ocean/Initialization.hh"
#include "pism/coupler/ocean/sea_level/Factory.hh"
#include "pism/coupler/ocean/sea_level/Initialization.hh"
#include "pism/coupler/surface/Factory.hh"
#include "pism/coupler/surface/Initialization.hh"
#include "pism/earth/LingleClark.hh"
#include "pism/earth/BedDef.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/util/Vars.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/projection.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/age/AgeModel.hh"
#include "pism/energy/EnthalpyModel.hh"
#include "pism/energy/TemperatureModel.hh"

namespace pism {

//! Initialize time from an input file or command-line options.
void IceModel::time_setup() {
  initialize_time(m_grid->com,
                  m_config->get_string("time.dimension_name"),
                  *m_log, *m_time);

  bool use_calendar = m_config->get_boolean("output.runtime.time_use_calendar");

  if (use_calendar) {
    m_log->message(2,
                   "* Run time: [%s, %s]  (%s years, using the '%s' calendar)\n",
                   m_time->start_date().c_str(),
                   m_time->end_date().c_str(),
                   m_time->run_length().c_str(),
                   m_time->calendar().c_str());
  } else {
    std::string time_units = m_config->get_string("output.runtime.time_unit_name");

    double
      start  = m_time->convert_time_interval(m_time->start(), time_units),
      end    = m_time->convert_time_interval(m_time->end(), time_units),
      length = end - start;

    m_log->message(2,
                   "* Run time: [%f %s, %f %s]  (%f %s)\n",
                   start, time_units.c_str(),
                   end, time_units.c_str(),
                   length, time_units.c_str());
  }
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

  // Check if we are initializing from a PISM output file:
  InputOptions input = process_input_options(m_ctx->com(), m_config);

  const bool use_input_file = input.type == INIT_BOOTSTRAP or input.type == INIT_RESTART;

  std::unique_ptr<PIO> input_file;

  if (use_input_file) {
    input_file.reset(new PIO(m_grid->com, "guess_mode", input.filename, PISM_READONLY));
  }

  // Get projection information and compute latitudes and longitudes *before* a component
  // decides to use them...
  {
    if (use_input_file) {
      std::string mapping_name = m_grid->get_mapping_info().mapping.get_name();
      MappingInfo info = get_projection_info(*input_file, mapping_name,
                                             m_ctx->unit_system());

      if (not info.proj4.empty()) {
        m_log->message(2, "* Got projection parameters \"%s\" from \"%s\".\n",
                       info.proj4.c_str(), input.filename.c_str());
      }

      m_output_global_attributes.set_string("proj4", info.proj4);
      m_grid->set_mapping_info(info);

      std::string history = input_file->get_att_text("PISM_GLOBAL", "history");
      m_output_global_attributes.set_string("history",
                                            history + m_output_global_attributes.get_string("history"));

    }

    compute_lat_lon();
  }

  // Initialize 2D fields owned by IceModel (ice geometry, etc)
  {
    switch (input.type) {
    case INIT_RESTART:
      restart_2d(*input_file, input.record);
      break;
    case INIT_BOOTSTRAP:
      bootstrap_2d(*input_file);
      break;
    case INIT_OTHER:
    default:
      initialize_2d();
    }

    regrid();
  }

  m_sea_level->init(m_geometry);

  // Initialize a bed deformation model.
  if (m_beddef) {
    m_beddef->init(input, m_geometry.ice_thickness, m_sea_level->elevation());
    m_grid->variables().add(m_beddef->bed_elevation());
    m_grid->variables().add(m_beddef->uplift());
  }

  // Now ice thickness, bed elevation, and sea level are available, so we can compute the ice
  // surface elevation and the cell type mask. This also ensures consistency of ice geometry.
  //
  // The stress balance code is executed near the beginning of a step and ice geometry is
  // updated near the end, so during the second time step the stress balance code is
  // guaranteed not to see "icebergs". Here we make sure that the first time step is OK
  // too.
  enforce_consistency_of_geometry(REMOVE_ICEBERGS);

  // By now ice geometry is set (including regridding) and so we can initialize the ocean model,
  // which may need ice thickness, bed topography, and the cell type mask.
  {
    m_ocean->init(m_geometry);
  }

  // Now surface elevation is initialized, so we can initialize surface models (some use
  // elevation-based parameterizations of surface temperature and/or mass balance).
  m_surface->init(m_geometry);

  if (m_subglacial_hydrology) {

    switch (input.type) {
    case INIT_RESTART:
      m_subglacial_hydrology->restart(*input_file, input.record);
      break;
    case INIT_BOOTSTRAP:
      m_subglacial_hydrology->bootstrap(*input_file,
                                        m_geometry.ice_thickness);
      break;
    case INIT_OTHER:
      {
        IceModelVec2S
          &W_till = m_work2d[0],
          &W      = m_work2d[1],
          &P      = m_work2d[2];

        W_till.set(m_config->get_double("bootstrapping.defaults.tillwat"));
        W.set(m_config->get_double("bootstrapping.defaults.bwat"));
        P.set(m_config->get_double("bootstrapping.defaults.bwp"));

        m_subglacial_hydrology->initialize(W_till, W, P);
        break;
      }
    }
  }

  // basal_yield_stress_model->init() needs bwat so this must happen
  // after subglacial_hydrology->init()
  if (m_basal_yield_stress_model) {
    assert((bool)m_subglacial_hydrology);
    m_basal_yield_stress_model->init(m_geometry,
                                     m_subglacial_hydrology->till_water_thickness(),
                                     m_subglacial_hydrology->overburden_pressure());
  }

  // Initialize the bedrock thermal layer model.
  //
  // If
  // - PISM is bootstrapping and
  // - we are using a non-zero-thickness thermal layer
  //
  // initialization of m_btu requires the temperature at the top of the bedrock. This is a problem
  // because get_bed_top_temp() uses the enthalpy field that is not initialized until later and
  // bootstrapping enthalpy uses the flux through the bottom surface of the ice (top surface of the
  // bedrock) provided by m_btu.
  //
  // We get out of this by using the fact that the full state of m_btu is not needed and
  // bootstrapping of the temperature field can be delayed.
  //
  // Note that to bootstrap m_btu we use the steady state solution of the heat equation in columns
  // of the bedrock (a straight line at each column), so the flux through the top surface of the
  // bedrock after bootstrapping is the same as the time-independent geothermal flux applied at the
  // BOTTOM surface of the bedrock layer.
  //
  // The code then delays bootstrapping of the thickness field until the first time step.
  if (m_btu) {
    m_btu->init(input);
  }

  if (m_age_model != NULL) {
    m_age_model->init(input);
    m_grid->variables().add(m_age_model->age());
  }

  // Initialize the energy balance sub-model.
  {
    switch (input.type) {
    case INIT_RESTART:
      {
        m_energy_model->restart(*input_file, input.record);
        break;
      }
    case INIT_BOOTSTRAP:
      {

        m_energy_model->bootstrap(*input_file,
                                  m_geometry.ice_thickness,
                                  m_surface->temperature(),
                                  m_surface->mass_flux(),
                                  m_btu->flux_through_top_surface());
        break;
      }
    case INIT_OTHER:
    default:
      {
        m_basal_melt_rate.set(m_config->get_double("bootstrapping.defaults.bmelt"));

        m_energy_model->initialize(m_basal_melt_rate,
                                   m_geometry.ice_thickness,
                                   m_surface->temperature(),
                                   m_surface->mass_flux(),
                                   m_btu->flux_through_top_surface());

      }
    }
    m_grid->variables().add(m_energy_model->enthalpy());
  }

  // this has to go after we add enthalpy to m_grid->variables()
  if (m_stress_balance) {
    m_stress_balance->init();
  }

  // miscellaneous steps
  {
    reset_counters();

    auto startstr = pism::printf("PISM (%s) started on %d procs.",
                                 PISM_Revision, (int)m_grid->size());
    prepend_history(startstr + args_string());
  }
}

//! Initialize 2D model state fields managed by IceModel from a file (for re-starting).
/*!
 * This method should eventually go away as IceModel turns into a "coupler" and all physical
 * processes are handled by sub-models.
 */
void IceModel::restart_2d(const PIO &input_file, unsigned int last_record) {
  std::string filename = input_file.inq_filename();

  m_log->message(2, "initializing 2D fields from NetCDF file '%s'...\n", filename.c_str());

  for (auto variable : m_model_state) {
    variable->read(input_file, last_record);
  }
}

void IceModel::bootstrap_2d(const PIO &input_file) {

  m_log->message(2, "bootstrapping from file '%s'...\n", input_file.inq_filename().c_str());

  std::string usurf_name;
  bool usurf_found = false, mask_found = false, usurf_found_by_std_name = false;
  input_file.inq_var("usurf", "surface_altitude",
                     usurf_found, usurf_name, usurf_found_by_std_name);
  mask_found = input_file.inq_var("mask");

  // now work through all the 2d variables, regridding if present and otherwise
  // setting to default values appropriately

  if (mask_found) {
    m_log->message(2, "  WARNING: 'mask' found; IGNORING IT!\n");
  }

  if (usurf_found) {
    m_log->message(2, "  WARNING: surface elevation 'usurf' found; IGNORING IT!\n");
  }

  m_log->message(2, "  reading 2D model state variables by regridding ...\n");

  // longitude
  {
    m_geometry.longitude.regrid(input_file, OPTIONAL);

    std::string lon_name;
    bool lon_found = false, lon_found_by_std_name = false;
    input_file.inq_var("lon", "longitude", lon_found, lon_name, lon_found_by_std_name);

    if (not lon_found) {
      m_geometry.longitude.metadata().set_string("missing_at_bootstrap", "true");
    }
  }

  // latitude
  {
    m_geometry.latitude.regrid(input_file, OPTIONAL);

    std::string lat_name;
    bool lat_found = false, lat_found_by_std_name = false;
    input_file.inq_var("lat", "latitude",  lat_found, lat_name, lat_found_by_std_name);

    if (not lat_found) {
      m_geometry.latitude.metadata().set_string("missing_at_bootstrap", "true");
    }
  }

  m_geometry.ice_thickness.regrid(input_file, OPTIONAL,
                                  m_config->get_double("bootstrapping.defaults.ice_thickness"));
  // check the range of the ice thickness
  {
    Range thk_range = m_geometry.ice_thickness.range();

    if (thk_range.max >= m_grid->Lz() + 1e-6) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Maximum ice thickness (%f meters)\n"
                                    "exceeds the height of the computational domain (%f meters).",
                                    thk_range.max, m_grid->Lz());
    }
  }

  if (m_config->get_boolean("geometry.part_grid.enabled")) {
    // Read the Href field from an input file. This field is
    // grid-dependent, so interpolating it from one grid to a
    // different one does not make sense in general.
    // (IceModel::Href_cleanup() will take care of the side effects of
    // such interpolation, though.)
    //
    // On the other hand, we need to read it in to be able to re-start
    // from a PISM output file using the -bootstrap option.
    m_geometry.ice_area_specific_volume.regrid(input_file, OPTIONAL, 0.0);
  }

  if (m_config->get_boolean("stress_balance.ssa.dirichlet_bc")) {
    // Do not use Dirichlet B.C. anywhere if bc_mask is not present.
    m_ssa_dirichlet_bc_mask.regrid(input_file, OPTIONAL, 0.0);
    // In the absence of u_ssa_bc and v_ssa_bc in the file the only B.C. that make sense are the
    // zero Dirichlet B.C.
    m_ssa_dirichlet_bc_values.regrid(input_file, OPTIONAL,  0.0);
  } else {
    m_ssa_dirichlet_bc_mask.set(0.0);
    m_ssa_dirichlet_bc_values.set(0.0);
  }

  // check if Lz is valid
  Range thk_range = m_geometry.ice_thickness.range();

  if (thk_range.max > m_grid->Lz()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Max. ice thickness (%3.3f m)\n"
                                  "exceeds the height of the computational domain (%3.3f m).",
                                  thk_range.max, m_grid->Lz());
  }
}

void IceModel::initialize_2d() {
  // This method should NOT have the "noreturn" attribute. (This attribute does not mix with virtual
  // methods).
  throw RuntimeError(PISM_ERROR_LOCATION, "cannot initialize IceModel without an input file");
}

//! Manage regridding based on user options.
void IceModel::regrid() {

  auto filename    = m_config->get_string("input.regrid.file");
  auto regrid_vars = set_split(m_config->get_string("input.regrid.vars"), ',');

  // Return if no regridding is requested:
  if (filename.empty()) {
     return;
  }

  m_log->message(2, "regridding from file %s ...\n", filename.c_str());

  {
    PIO regrid_file(m_grid->com, "guess_mode", filename, PISM_READONLY);
    for (auto v : m_model_state) {
      if (regrid_vars.find(v->get_name()) != regrid_vars.end()) {
        v->regrid(regrid_file, CRITICAL);
      }
    }

    // Check the range of the ice thickness.
    {
      double
        max_thickness = m_geometry.ice_thickness.range().max,
        Lz            = m_grid->Lz();

      if (max_thickness >= Lz + 1e-6) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION, "Maximum ice thickness (%f meters)\n"
                                      "exceeds the height of the computational domain (%f meters).",
                                      max_thickness, Lz);
      }
    }
  }
}

//! \brief Decide which stress balance model to use.
void IceModel::allocate_stressbalance() {

  if (m_stress_balance) {
    return;
  }

  m_log->message(2, "# Allocating a stress balance model...\n");

  // false means "not regional"
  m_stress_balance = stressbalance::create(m_config->get_string("stress_balance.model"),
                                           m_grid, false);

  m_submodels["stress balance"] = m_stress_balance.get();
}

void IceModel::allocate_geometry_evolution() {
  if (m_geometry_evolution) {
    return;
  }

  m_log->message(2, "# Allocating the geometry evolution model...\n");

  m_geometry_evolution.reset(new GeometryEvolution(m_grid));

  m_submodels["geometry_evolution"] = m_geometry_evolution.get();
}


void IceModel::allocate_iceberg_remover() {

  if (m_iceberg_remover != NULL) {
    return;
  }

  m_log->message(2,
             "# Allocating an iceberg remover (part of a calving model)...\n");

  if (m_config->get_boolean("geometry.remove_icebergs")) {

    // this will throw an exception on failure
    m_iceberg_remover = new calving::IcebergRemover(m_grid);

    // Iceberg Remover does not have a state, so it is OK to
    // initialize here.
    m_iceberg_remover->init();

    m_submodels["iceberg remover"] = m_iceberg_remover;
  }
}

void IceModel::allocate_age_model() {

  if (m_age_model != NULL) {
    return;
  }

  if (m_config->get_boolean("age.enabled")) {
    m_log->message(2, "# Allocating an ice age model...\n");

    if (m_stress_balance == NULL) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "Cannot allocate an age model: m_stress_balance == NULL.");
    }

    m_age_model = new AgeModel(m_grid, m_stress_balance.get());
    m_submodels["age model"] = m_age_model;
  }
}

void IceModel::allocate_energy_model() {

  if (m_energy_model != NULL) {
    return;
  }

  m_log->message(2, "# Allocating an energy balance model...\n");

  if (m_config->get_boolean("energy.enabled")) {
    if (m_config->get_boolean("energy.temperature_based")) {
      m_energy_model = new energy::TemperatureModel(m_grid, m_stress_balance.get());
    } else {
      m_energy_model = new energy::EnthalpyModel(m_grid, m_stress_balance.get());
    }
  } else {
    m_energy_model = new energy::DummyEnergyModel(m_grid, m_stress_balance.get());
  }

  m_submodels["energy balance model"] = m_energy_model;
}


//! \brief Decide which bedrock thermal unit to use.
void IceModel::allocate_bedrock_thermal_unit() {

  if (m_btu != NULL) {
    return;
  }

  m_log->message(2, "# Allocating a bedrock thermal layer model...\n");

  m_btu = energy::BedThermalUnit::FromOptions(m_grid, m_ctx);

  m_submodels["bedrock thermal model"] = m_btu;
}

//! \brief Decide which subglacial hydrology model to use.
void IceModel::allocate_subglacial_hydrology() {

  using namespace pism::hydrology;

  std::string hydrology_model = m_config->get_string("hydrology.model");

  if (m_subglacial_hydrology) { // indicates it has already been allocated
    return;
  }

  m_log->message(2, "# Allocating a subglacial hydrology model...\n");

  if (hydrology_model == "null") {
    m_subglacial_hydrology.reset(new NullTransport(m_grid));
  } else if (hydrology_model == "routing") {
    m_subglacial_hydrology.reset(new Routing(m_grid));
  } else if (hydrology_model == "distributed") {
    m_subglacial_hydrology.reset(new Distributed(m_grid));
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "unknown 'hydrology.model': %s", hydrology_model.c_str());
  }

  m_submodels["subglacial hydrology"] = m_subglacial_hydrology.get();
}

//! \brief Decide which basal yield stress model to use.
void IceModel::allocate_basal_yield_stress() {

  if (m_basal_yield_stress_model) {
    return;
  }

  m_log->message(2,
             "# Allocating a basal yield stress model...\n");

  std::string model = m_config->get_string("stress_balance.model");

  // only these two use the yield stress (so far):
  if (model == "ssa" || model == "ssa+sia") {
    std::string yield_stress_model = m_config->get_string("basal_yield_stress.model");

    if (yield_stress_model == "constant") {
      m_basal_yield_stress_model.reset(new ConstantYieldStress(m_grid));
    } else if (yield_stress_model == "mohr_coulomb") {
      m_basal_yield_stress_model.reset(new MohrCoulombYieldStress(m_grid));
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "yield stress model '%s' is not supported.",
                                    yield_stress_model.c_str());
    }

    m_submodels["basal yield stress"] = m_basal_yield_stress_model.get();
  }
}

//! Allocate PISM's sub-models implementing some physical processes.
/*!
  This method is called after memory allocation but before filling any of
  IceModelVecs because all the physical parameters should be initialized before
  setting up the coupling or filling model-state variables.
 */
void IceModel::allocate_submodels() {

  allocate_geometry_evolution();

  allocate_iceberg_remover();

  allocate_stressbalance();

  // this has to happen *after* allocate_stressbalance()
  {
    allocate_age_model();
    allocate_energy_model();
    allocate_subglacial_hydrology();
  }

  // this has to happen *after* allocate_subglacial_hydrology()
  allocate_basal_yield_stress();

  allocate_bedrock_thermal_unit();

  allocate_bed_deformation();

  allocate_couplers();
}


void IceModel::allocate_couplers() {
  // Initialize boundary models:

  if (not m_surface) {

    m_log->message(2, "# Allocating a surface process model or coupler...\n");

    surface::Factory ps(m_grid, atmosphere::Factory(m_grid).create());

    m_surface.reset(new surface::InitializationHelper(m_grid, ps.create()));

    m_submodels["surface process model"] = m_surface.get();
  }

  if (not m_sea_level) {
    m_log->message(2, "# Allocating sea level forcing...\n");

    using namespace ocean::sea_level;

    m_sea_level.reset(new InitializationHelper(m_grid, Factory(m_grid).create()));

    m_submodels["sea level forcing"] = m_sea_level.get();
  }

  if (not m_ocean) {
    m_log->message(2, "# Allocating an ocean model or coupler...\n");

    using namespace ocean;

    m_ocean.reset(new InitializationHelper(m_grid, Factory(m_grid).create()));

    m_submodels["ocean model"] = m_ocean.get();
  }
}

//! Miscellaneous initialization tasks plus tasks that need the fields that can come from regridding.
void IceModel::misc_setup() {

  m_log->message(3, "Finishing initialization...\n");
  InputOptions opts = process_input_options(m_ctx->com(), m_config);

  if (not (opts.type == INIT_OTHER)) {
    // initializing from a file
    PIO nc(m_grid->com, "guess_mode", opts.filename, PISM_READONLY);

    std::string source = nc.get_att_text("PISM_GLOBAL", "source");

    if (opts.type == INIT_RESTART) {
      // If it's missing, print a warning
      if (source.empty()) {
        m_log->message(1,
                       "PISM WARNING: file '%s' does not have the 'source' global attribute.\n"
                       "     If '%s' is a PISM output file, please run the following to get rid of this warning:\n"
                       "     ncatted -a source,global,c,c,PISM %s\n",
                       opts.filename.c_str(), opts.filename.c_str(), opts.filename.c_str());
      } else if (source.find("PISM") == std::string::npos) {
        // If the 'source' attribute does not contain the string "PISM", then print
        // a message and stop:
        m_log->message(1,
                       "PISM WARNING: '%s' does not seem to be a PISM output file.\n"
                       "     If it is, please make sure that the 'source' global attribute contains the string \"PISM\".\n",
                       opts.filename.c_str());
      }
    }
  }

  {
    // A single record of a time-dependent variable cannot exceed 2^32-4
    // bytes in size. See the NetCDF User's Guide
    // <http://www.unidata.ucar.edu/software/netcdf/docs/netcdf.html#g_t64-bit-Offset-Limitations>.
    // Here we use "long int" to avoid integer overflow.
    const long int two_to_thirty_two = 4294967296L;
    const long int
      Mx = m_grid->Mx(),
      My = m_grid->My(),
      Mz = m_grid->Mz();
    std::string output_format = m_config->get_string("output.format");
    if (Mx * My * Mz * sizeof(double) > two_to_thirty_two - 4 and
        (output_format == "netcdf3" or output_format == "pnetcdf")) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "The computational grid is too big to fit in a NetCDF-3 file.\n"
                                    "Each 3D variable requires %lu Mb.\n"
                                    "Please use re-build PISM with parallel NetCDF-4\n"
                                    "and use '-o_format netcdf4_parallel' to proceed.",
                                    Mx * My * Mz * sizeof(double) / (1024 * 1024));
    }
  }

  m_output_vars = output_variables(m_config->get_string("output.size"));

#if (PISM_USE_PROJ4==1)
  {
    std::string proj_string = m_grid->get_mapping_info().proj4;
    if (not proj_string.empty()) {
      m_output_vars.insert("lon_bnds");
      m_output_vars.insert("lat_bnds");
      m_geometry.latitude.metadata().set_string("bounds", "lat_bnds");
      m_geometry.longitude.metadata().set_string("bounds", "lon_bnds");
    }
  }
#endif

  init_calving();
  init_diagnostics();
  init_snapshots();
  init_backups();
  init_timeseries();
  init_extras();

  // Make sure that we use the output.variable_order that works with
  // NetCDF-4 I/O.
  std::string o_format = m_config->get_string("output.format");
  if (o_format == "netcdf4_parallel" and m_config->get_string("output.variable_order") != "yxz") {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "output format netcdf4_parallel requires -o_order yxz.");
  }

  // a report on whether PISM-PIK modifications of IceModel are in use
  {
    std::vector<std::string> pik_methods;
    if (m_config->get_boolean("geometry.part_grid.enabled")) {
      pik_methods.push_back("part_grid");
    }
    if (m_config->get_boolean("geometry.remove_icebergs")) {
      pik_methods.push_back("kill_icebergs");
    }

    if (not pik_methods.empty()) {
      m_log->message(2,
                     "* PISM-PIK mass/geometry methods are in use: %s\n",
                     join(pik_methods, ", ").c_str());
    }
  }

  // initialize diagnostics
  {
    // reset: this gives diagnostics a chance to capture the current state of the model at the
    // beginning of the run
    for (auto d : m_diagnostics) {
      d.second->reset();
    }

    // read in the state (accumulators) if we are re-starting a run
    if (opts.type == INIT_RESTART) {
      PIO file(m_grid->com, "guess_mode", opts.filename, PISM_READONLY);
      for (auto d : m_diagnostics) {
        d.second->init(file, opts.record);
      }
    }
  }

  if (m_surface_input_for_hydrology) {
    m_surface_input_for_hydrology->init(m_config->get_string("hydrology.surface_input_file"),
                                        0, 0);
  }
}

//! \brief Initialize calving mechanisms.
void IceModel::init_calving() {

  std::set<std::string> methods = set_split(m_config->get_string("calving.methods"), ',');

  if (methods.find("ocean_kill") != methods.end()) {

    if (not m_ocean_kill_calving) {
      m_ocean_kill_calving = new calving::OceanKill(m_grid);
    }

    m_ocean_kill_calving->init();
    methods.erase("ocean_kill");

    m_submodels["ocean kill calving"] = m_ocean_kill_calving;
  }

  if (methods.find("thickness_calving") != methods.end()) {

    if (not m_thickness_threshold_calving) {
      m_thickness_threshold_calving = new calving::CalvingAtThickness(m_grid);
    }

    m_thickness_threshold_calving->init();
    methods.erase("thickness_calving");

    m_submodels["thickness threshold calving"] = m_thickness_threshold_calving;
  }


  if (methods.find("eigen_calving") != methods.end()) {

    if (not m_eigen_calving) {
      m_eigen_calving = new calving::EigenCalving(m_grid);
    }

    m_eigen_calving->init();
    methods.erase("eigen_calving");

    m_submodels["eigen calving"] = m_eigen_calving;
  }

  if (methods.find("vonmises_calving") != methods.end()) {

    if (not m_vonmises_calving) {
      m_vonmises_calving = new calving::vonMisesCalving(m_grid,
                                                        m_stress_balance->shallow()->flow_law());
    }

    m_vonmises_calving->init();
    methods.erase("vonmises_calving");

    m_submodels["von Mises calving"] = m_vonmises_calving;
  }

  if (methods.find("frontal_melt") != methods.end()) {

    if (not m_frontal_melt) {
      m_frontal_melt = new FrontalMelt(m_grid);
    }

    m_frontal_melt->init();
    methods.erase("frontal_melt");

    m_submodels["frontal melt"] = m_frontal_melt;
  }

  if (methods.find("float_kill") != methods.end()) {
    if (not m_float_kill_calving) {
      m_float_kill_calving = new calving::FloatKill(m_grid);
    }

    m_float_kill_calving->init();
    methods.erase("float_kill");

    m_submodels["float kill calving"] = m_float_kill_calving;
  }

  if (not methods.empty()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "PISM ERROR: calving method(s) [%s] are not supported.\n",
                                  set_join(methods, ",").c_str());
  }
}

void IceModel::allocate_bed_deformation() {
  std::string model = m_config->get_string("bed_deformation.model");

  if (m_beddef != NULL) {
    return;
  }

  m_log->message(2,
                 "# Allocating a bed deformation model...\n");

  if (model == "none") {
    m_beddef = new bed::Null(m_grid);
  }
  else if (model == "iso") {
    m_beddef = new bed::PointwiseIsostasy(m_grid);
  }
  else if (model == "lc") {
    m_beddef = new bed::LingleClark(m_grid);
  }

  m_submodels["bed deformation"] = m_beddef;
}

//! Read some runtime (command line) options and alter the
//! corresponding parameters or flags as appropriate.
void IceModel::process_options() {

  m_log->message(3,
             "Processing physics-related command-line options...\n");

  set_config_from_options(*m_config);

  // Set global attributes using the config database:
  m_output_global_attributes.set_string("title", m_config->get_string("run_info.title"));
  m_output_global_attributes.set_string("institution", m_config->get_string("run_info.institution"));
  m_output_global_attributes.set_string("command", args_string());

  // warn about some option combinations

  if (m_config->get_double("time_stepping.maximum_time_step") <= 0) {
    throw RuntimeError(PISM_ERROR_LOCATION, "time_stepping.maximum_time_step has to be greater than 0.");
  }

  if (not m_config->get_boolean("geometry.update.enabled") &&
      m_config->get_boolean("time_stepping.skip.enabled")) {
    m_log->message(2,
               "PISM WARNING: Both -skip and -no_mass are set.\n"
               "              -skip only makes sense in runs updating ice geometry.\n");
  }

  if (m_config->get_string("calving.methods").find("thickness_calving") != std::string::npos &&
      not m_config->get_boolean("geometry.part_grid.enabled")) {
    m_log->message(2,
               "PISM WARNING: Calving at certain terminal ice thickness (-calving thickness_calving)\n"
               "              without application of partially filled grid cell scheme (-part_grid)\n"
               "              may lead to (incorrect) non-moving ice shelf front.\n");
  }
}

//! Assembles a list of diagnostics corresponding to an output file size.
std::set<std::string> IceModel::output_variables(const std::string &keyword) {

  std::string variables;

  if (keyword == "none" or
      keyword == "small") {
    variables = "";
  } else if (keyword == "medium") {
    variables = m_config->get_string("output.sizes.medium");
  } else if (keyword == "big_2d") {
    variables = (m_config->get_string("output.sizes.medium") + "," +
                 m_config->get_string("output.sizes.big_2d"));
  } else if (keyword == "big") {
    variables = (m_config->get_string("output.sizes.medium") + "," +
                 m_config->get_string("output.sizes.big_2d") + "," +
                 m_config->get_string("output.sizes.big"));
  }

  return set_split(variables, ',');
}

void IceModel::compute_lat_lon() {

  std::string projection = m_grid->get_mapping_info().proj4;

  if (m_config->get_boolean("grid.recompute_longitude_and_latitude") and
      not projection.empty()) {
    m_log->message(2,
                   "* Computing longitude and latitude using projection parameters...\n");

    compute_longitude(projection, m_geometry.longitude);
    m_geometry.longitude.metadata().set_string("missing_at_bootstrap", "");
    compute_latitude(projection, m_geometry.latitude);
    m_geometry.latitude.metadata().set_string("missing_at_bootstrap", "");
  }
}

} // end of namespace pism
