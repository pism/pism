// Copyright (C) 2009--2016 Ed Bueler and Constantine Khroulev
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

#include <algorithm>

#include "iceModel.hh"
#include "base/basalstrength/PISMConstantYieldStress.hh"
#include "base/basalstrength/PISMMohrCoulombYieldStress.hh"
#include "base/basalstrength/basal_resistance.hh"
#include "base/calving/CalvingAtThickness.hh"
#include "base/calving/EigenCalving.hh"
#include "base/calving/vonMisesCalving.hh"
#include "base/calving/FloatKill.hh"
#include "base/calving/IcebergRemover.hh"
#include "base/calving/OceanKill.hh"
#include "base/calving/FrontalMelt.hh"
#include "base/energy/BedThermalUnit.hh"
#include "base/hydrology/PISMHydrology.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/stressbalance/sia/SIAFD.hh"
#include "base/stressbalance/ssa/SSAFD.hh"
#include "base/stressbalance/ssa/SSAFEM.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMTime.hh"
#include "base/util/error_handling.hh"
#include "base/util/io/PIO.hh"
#include "base/util/pism_options.hh"
#include "coupler/PISMOcean.hh"
#include "coupler/PISMSurface.hh"
#include "coupler/atmosphere/PAFactory.hh"
#include "coupler/ocean/POFactory.hh"
#include "coupler/ocean/POInitialization.hh"
#include "coupler/surface/PSFactory.hh"
#include "coupler/surface/PSInitialization.hh"
#include "earth/PBLingleClark.hh"
#include "earth/PISMBedDef.hh"
#include "enthalpyConverter.hh"
#include "base/util/PISMVars.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/projection.hh"
#include "base/util/pism_utilities.hh"
#include "base/age/AgeModel.hh"
#include "base/energy/EnthalpyModel.hh"
#include "base/energy/TemperatureModel.hh"

namespace pism {

//! Initialize time from an input file or command-line options.
void IceModel::time_setup() {
  initialize_time(m_grid->com,
                  m_config->get_string("time.dimension_name"),
                  *m_log, *m_time);

  m_log->message(2,
             "* Run time: [%s, %s]  (%s years, using the '%s' calendar)\n",
             m_time->start_date().c_str(),
             m_time->end_date().c_str(),
             m_time->run_length().c_str(),
             m_time->calendar().c_str());
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
  InputOptions input = process_input_options(m_ctx->com());

  const bool use_input_file = input.type == INIT_BOOTSTRAP or input.type == INIT_RESTART;

  PISM_SHARED_PTR(PIO) input_file;

  if (use_input_file) {
    input_file.reset(new PIO(m_grid->com, "guess_mode", input.filename, PISM_READONLY));
  }

  // Get projection information and compute cell areas and lat/lon *before* a component decides to
  // use latitude or longitude...
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

      initialize_cumulative_fluxes(*input_file);
    } else {
      reset_cumulative_fluxes();
    }

    compute_cell_areas();
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

    regrid(2);
  }

  // By now ice geometry is set (including regridding) and so we can initialize the ocean model,
  // which may need ice thickness to bootstrap.
  {
    m_log->message(2, "* Initializing the ocean model...\n");
    m_ocean->init();
  }

  // Initialize a bed deformation model. This may use ice thickness initialized above.
  if (m_beddef) {
    m_beddef->init();
    m_grid->variables().add(m_beddef->bed_elevation());
    m_grid->variables().add(m_beddef->uplift());
  }

  // Now ice thickness, bed elevation, and sea level are available, so we can compute the ice
  // surface elevation and the cell type mask. This also ensures consistency of ice geometry.
  enforce_consistency_of_geometry();

  // Now surface elevation is initialized, so we can initialize surface models (some use
  // elevation-based parameterizations of surface temperature and/or mass balance).
  m_surface->init();

  if (m_subglacial_hydrology) {
    m_subglacial_hydrology->init();
  }

  // basal_yield_stress_model->init() needs bwat so this must happen
  // after subglacial_hydrology->init()
  if (m_basal_yield_stress_model) {
    m_basal_yield_stress_model->init();
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
    IceModelVec2S &ice_surface_temperature = m_work2d[0];
    IceModelVec2S &climatic_mass_balance   = m_work2d[1];

    switch (input.type) {
    case INIT_RESTART:
      {
        m_energy_model->restart(*input_file, input.record);
        break;
      }
    case INIT_BOOTSTRAP:
      {
        m_surface->ice_surface_temperature(ice_surface_temperature);
        m_surface->ice_surface_mass_flux(climatic_mass_balance);
        m_energy_model->bootstrap(*input_file,
                                  m_ice_thickness,
                                  ice_surface_temperature,
                                  climatic_mass_balance,
                                  m_btu->flux_through_top_surface());
        break;
      }
    case INIT_OTHER:
    default:
      {
        m_basal_melt_rate.set(m_config->get_double("bootstrapping.defaults.bmelt"));
        m_surface->ice_surface_temperature(ice_surface_temperature);
        m_surface->ice_surface_mass_flux(climatic_mass_balance);
        m_energy_model->initialize(m_basal_melt_rate,
                                   m_ice_thickness,
                                   ice_surface_temperature,
                                   climatic_mass_balance,
                                   m_btu->flux_through_top_surface());

      }
    }
    m_grid->variables().add(m_energy_model->enthalpy());
  }

  // this has to go after we add enthalpy to m_grid->variables()
  if (m_stress_balance) {
    m_stress_balance->init();

    if (m_config->get_boolean("geometry.update.use_basal_melt_rate")) {
      m_stress_balance->set_basal_melt_rate(m_basal_melt_rate);
    }
  }

  // miscellaneous steps
  {
    reset_counters();
    stampHistoryCommand();
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

  // Read the model state, mapping and climate_steady variables:
  std::set<std::string> vars = m_grid->variables().keys();

  std::set<std::string>::iterator i;
  for (i = vars.begin(); i != vars.end(); ++i) {
    // FIXME: remove const_cast. This is bad.
    IceModelVec *var = const_cast<IceModelVec*>(m_grid->variables().get(*i));
    SpatialVariableMetadata &m = var->metadata();

    std::string
      intent     = m.get_string("pism_intent"),
      short_name = m.get_string("short_name");

    if (intent == "model_state" ||
        intent == "mapping"     ||
        intent == "climate_steady") {

      // skip "enthalpy" and "Href" for now: we'll take care
      // of them a little later
      if (short_name == "enthalpy" ||
          short_name == "Href") {
        continue;
      }

      var->read(input_file, last_record);
    }
  }

  // check if the input file has Href; set to 0 if it is not present
  if (m_config->get_boolean("geometry.part_grid.enabled")) {

    if (input_file.inq_var("Href")) {
      m_Href.read(input_file, last_record);
    } else {
      m_log->message(2,
                     "PISM WARNING: Href for PISM-PIK -part_grid not found in '%s'."
                     " Setting it to zero...\n",
                     filename.c_str());
      m_Href.set(0.0);
    }
  }
}

void IceModel::initialize_2d() {
  throw RuntimeError(PISM_ERROR_LOCATION, "cannot initialize IceModel without an input file");
}

void IceModel::reset_cumulative_fluxes() {
  // 2D
  m_cumulative_flux_fields.reset();
  // scalar
  m_cumulative_fluxes = FluxCounters();
}


void IceModel::initialize_cumulative_fluxes(const PIO &input_file) {
  // 2D
  m_cumulative_flux_fields.regrid(input_file);

  // scalar, stored in run_stats
  if (input_file.inq_var("run_stats")) {
    io::read_attributes(input_file, m_run_stats.get_name(), m_run_stats);

    try {
      m_cumulative_fluxes.H_to_Href      = m_run_stats.get_double("H_to_Href_flux_cumulative");
      m_cumulative_fluxes.Href_to_H      = m_run_stats.get_double("Href_to_H_flux_cumulative");
      m_cumulative_fluxes.discharge      = m_run_stats.get_double("discharge_flux_cumulative");
      m_cumulative_fluxes.grounded_basal = m_run_stats.get_double("grounded_basal_ice_flux_cumulative");
      m_cumulative_fluxes.nonneg_rule    = m_run_stats.get_double("nonneg_rule_flux_cumulative");
      m_cumulative_fluxes.sub_shelf      = m_run_stats.get_double("sub_shelf_ice_flux_cumulative");
      m_cumulative_fluxes.sum_divQ_SIA   = m_run_stats.get_double("sum_divQ_SIA_cumulative");
      m_cumulative_fluxes.sum_divQ_SSA   = m_run_stats.get_double("sum_divQ_SSA_cumulative");
      m_cumulative_fluxes.surface        = m_run_stats.get_double("surface_ice_flux_cumulative");
    }
    catch (RuntimeError &e) {
      e.add_context("initializing cumulative flux counters from '%s'",
                    input_file.inq_filename().c_str());
      throw;
    }
  }
}


//! \brief Decide which stress balance model to use.
void IceModel::allocate_stressbalance() {

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  using namespace pism::stressbalance;

  if (m_stress_balance != NULL) {
    return;
  }

  m_log->message(2,
             "# Allocating a stress balance model...\n");

  std::string model = m_config->get_string("stress_balance.model");

  ShallowStressBalance *sliding = NULL;
  if (model == "none" || model == "sia") {
    sliding = new ZeroSliding(m_grid);
  } else if (model == "prescribed_sliding" || model == "prescribed_sliding+sia") {
    sliding = new PrescribedSliding(m_grid);
  } else if (model == "ssa" || model == "ssa+sia") {
    std::string method = m_config->get_string("stress_balance.ssa.method");

    if (method == "fem") {
      sliding = new SSAFEM(m_grid);
    } else if (method == "fd") {
      sliding = new SSAFD(m_grid);
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid ssa method: %s", method.c_str());
    }

  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid stress balance model: %s", model.c_str());
  }

  SSB_Modifier *modifier = NULL;
  if (model == "none" || model == "ssa" || model == "prescribed_sliding") {
    modifier = new ConstantInColumn(m_grid);
  } else if (model == "prescribed_sliding+sia" || model == "ssa+sia" || model == "sia") {
    modifier = new SIAFD(m_grid);
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid stress balance model: %s", model.c_str());
  }

  // ~StressBalance() will de-allocate sliding and modifier.
  m_stress_balance = new StressBalance(m_grid, sliding, modifier);

  m_submodels["stress balance"] = m_stress_balance;
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

    m_age_model = new AgeModel(m_grid, m_stress_balance);
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
      m_energy_model = new energy::TemperatureModel(m_grid, m_stress_balance);
    } else {
      m_energy_model = new energy::EnthalpyModel(m_grid, m_stress_balance);
    }
  } else {
    m_energy_model = new energy::DummyEnergyModel(m_grid, m_stress_balance);
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

  if (m_subglacial_hydrology != NULL) { // indicates it has already been allocated
    return;
  }

  m_log->message(2,
                 "# Allocating a subglacial hydrology model...\n");

  if (hydrology_model == "null") {
    m_subglacial_hydrology = new NullTransport(m_grid);
  } else if (hydrology_model == "routing") {
    m_subglacial_hydrology = new Routing(m_grid);
  } else if (hydrology_model == "distributed") {
    m_subglacial_hydrology = new Distributed(m_grid, m_stress_balance);
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "unknown value for configuration string 'hydrology.model':\n"
                                  "has value '%s'", hydrology_model.c_str());
  }

  m_submodels["subglacial hydrology"] = m_subglacial_hydrology;
}

//! \brief Decide which basal yield stress model to use.
void IceModel::allocate_basal_yield_stress() {

  if (m_basal_yield_stress_model != NULL) {
    return;
  }

  m_log->message(2,
             "# Allocating a basal yield stress model...\n");

  std::string model = m_config->get_string("stress_balance.model");

  // only these two use the yield stress (so far):
  if (model == "ssa" || model == "ssa+sia") {
    std::string yield_stress_model = m_config->get_string("basal_yield_stress.model");

    if (yield_stress_model == "constant") {
      m_basal_yield_stress_model = new ConstantYieldStress(m_grid);
    } else if (yield_stress_model == "mohr_coulomb") {
      m_basal_yield_stress_model = new MohrCoulombYieldStress(m_grid, m_subglacial_hydrology);
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "yield stress model '%s' is not supported.",
                                    yield_stress_model.c_str());
    }

    m_submodels["basal yield stress"] = m_basal_yield_stress_model;
  }
}

//! Allocate PISM's sub-models implementing some physical processes.
/*!
  This method is called after memory allocation but before filling any of
  IceModelVecs because all the physical parameters should be initialized before
  setting up the coupling or filling model-state variables.
 */
void IceModel::allocate_submodels() {

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
  atmosphere::Factory pa(m_grid);
  surface::Factory ps(m_grid);
  ocean::Factory po(m_grid);
  atmosphere::AtmosphereModel *atmosphere;

  if (m_surface == NULL) {

    m_log->message(2,
             "# Allocating a surface process model or coupler...\n");

    m_surface = new surface::InitializationHelper(m_grid, ps.create());

    atmosphere = pa.create();
    m_surface->attach_atmosphere_model(atmosphere);

    m_submodels["surface process model"] = m_surface;
  }

  if (m_ocean == NULL) {
    m_log->message(2,
             "# Allocating an ocean model or coupler...\n");

    m_ocean = new ocean::InitializationHelper(m_grid, po.create());

    m_submodels["ocean model"] = m_ocean;
  }
}

//! Miscellaneous initialization tasks plus tasks that need the fields that can come from regridding.
void IceModel::misc_setup() {

  m_log->message(3, "Finishing initialization...\n");
  InputOptions opts = process_input_options(m_ctx->com());

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
                                    "Please use '-o_format quilt' or re-build PISM with parallel NetCDF-4 or HDF5\n"
                                    "and use '-o_format netcdf4_parallel' or '-o_format hdf5' to proceed.",
                                    Mx * My * Mz * sizeof(double) / (1024 * 1024));
    }
  }

  m_output_vars = output_size_from_option("-o_size", "Sets the 'size' of an output file.",
                                          "medium");

  init_calving();
  init_diagnostics();
  init_snapshots();
  init_backups();
  init_timeseries();
  init_extras();
  init_viewers();

  // Make sure that we use the output.variable_order that works with NetCDF-4,
  // "quilt", and HDF5 parallel I/O. (For different reasons, but mainly because
  // it is faster.)
  std::string o_format = m_config->get_string("output.format");
  if ((o_format == "netcdf4_parallel" || o_format == "quilt" || o_format == "hdf5") &&
      m_config->get_string("output.variable_order") != "yxz") {
    throw RuntimeError(PISM_ERROR_LOCATION, "output formats netcdf4_parallel, quilt, and hdf5 require -o_order yxz.");
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
}

//! \brief Initialize calving mechanisms.
void IceModel::init_calving() {

  std::set<std::string> methods = set_split(m_config->get_string("calving.methods"), ',');

  if (methods.find("ocean_kill") != methods.end()) {

    if (m_ocean_kill_calving == NULL) {
      m_ocean_kill_calving = new calving::OceanKill(m_grid);
    }

    m_ocean_kill_calving->init();
    methods.erase("ocean_kill");

    m_submodels["ocean kill calving"] = m_ocean_kill_calving;
  }

  if (methods.find("thickness_calving") != methods.end()) {

    if (m_thickness_threshold_calving == NULL) {
      m_thickness_threshold_calving = new calving::CalvingAtThickness(m_grid);
    }

    m_thickness_threshold_calving->init();
    methods.erase("thickness_calving");

    m_submodels["thickness threshold calving"] = m_thickness_threshold_calving;
  }


  if (methods.find("eigen_calving") != methods.end()) {

    if (m_eigen_calving == NULL) {
      m_eigen_calving = new calving::EigenCalving(m_grid, m_stress_balance);
    }

    m_eigen_calving->init();
    methods.erase("eigen_calving");

    m_submodels["eigen calving"] = m_eigen_calving;
  }

  if (methods.find("vonmises_calving") != methods.end()) {

    if (m_vonmises_calving == NULL) {
      m_vonmises_calving = new calving::vonMisesCalving(m_grid, m_stress_balance);
    }

    m_vonmises_calving->init();
    methods.erase("vonmises_calving");

    m_submodels["von Mises calving"] = m_vonmises_calving;
  }

  if (methods.find("frontal_melt") != methods.end()) {

    if (m_frontal_melt == NULL) {
      m_frontal_melt = new FrontalMelt(m_grid, m_ocean);
    }

    m_frontal_melt->init();
    methods.erase("frontal_melt");

    m_submodels["frontal melt"] = m_frontal_melt;
  }

  if (methods.find("float_kill") != methods.end()) {
    if (m_float_kill_calving == NULL) {
      m_float_kill_calving = new calving::FloatKill(m_grid);
    }

    m_float_kill_calving->init();
    methods.erase("float_kill");

    m_submodels["float kill calving"] = m_float_kill_calving;
  }

  std::set<std::string>::iterator j = methods.begin();
  std::string unused;
  while (j != methods.end()) {
    unused += (*j + ",");
    ++j;
  }

  if (not unused.empty()) {
    m_log->message(2,
               "PISM ERROR: calving method(s) [%s] are unknown and are ignored.\n",
               unused.c_str());
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
    m_beddef = new bed::PBNull(m_grid);
  }
  else if (model == "iso") {
    m_beddef = new bed::PBPointwiseIsostasy(m_grid);
  }
  else if (model == "lc") {
    m_beddef = new bed::PBLingleClark(m_grid);
  }

  m_submodels["bed deformation"] = m_beddef;
}

} // end of namespace pism
