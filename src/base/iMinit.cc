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
  InputOptions opts = process_input_options(m_ctx->com());

  const bool use_input_file = opts.type == INIT_BOOTSTRAP or opts.type == INIT_RESTART;

  PIO input_file(m_grid->com, "guess_mode");

  if (use_input_file) {
    input_file.open(opts.filename, PISM_READONLY);
  }

  // Initialize 2D fields owned by IceModel (ice geometry, etc)
  {
    switch (opts.type) {
    case INIT_RESTART:
      restart_2d(input_file, opts.record);
      break;
    case INIT_BOOTSTRAP:
      bootstrap_2d(input_file);
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
  updateSurfaceElevationAndMask();

  // Now surface elevation is initialized, so we can initialize surface models (some use
  // elevation-based parameterizations of surface temperature and/or mass balance).
  m_surface->init();

  if (m_stress_balance) {
    m_stress_balance->init();

    if (m_config->get_boolean("geometry.update.use_basal_melt_rate")) {
      m_stress_balance->set_basal_melt_rate(m_basal_melt_rate);
    }
  }

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
    m_btu->init(opts);
  }

  // Initialize 3D (age and energy balance) parts of IceModel.
  {
    switch (opts.type) {
    case INIT_RESTART:
      restart_3d(input_file, opts.record);
      break;
    case INIT_BOOTSTRAP:
      bootstrap_3d();
      break;
    case INIT_OTHER:
    default:
      initialize_3d();
    }

    regrid(3);
  }

  // get projection information and compute cell areas
  {
    if (use_input_file) {
      get_projection_info(input_file);

      std::string proj4_string = m_output_global_attributes.get_string("proj4");
      if (not proj4_string.empty()) {
        m_log->message(2, "* Got projection parameters \"%s\" from \"%s\".\n",
                       proj4_string.c_str(), opts.filename.c_str());
      }

      std::string history = input_file.get_att_text("PISM_GLOBAL", "history");
      m_output_global_attributes.set_string("history",
                                            history + m_output_global_attributes.get_string("history"));
    }

    compute_cell_areas();
  }

  // miscellaneous steps
  {
    reset_counters();

    if (use_input_file) {
      std::string history = input_file.get_att_text("PISM_GLOBAL", "history");
      m_output_global_attributes.set_string("history",
                                            history + m_output_global_attributes.get_string("history"));

      initialize_cumulative_fluxes(input_file);
    } else {
      reset_cumulative_fluxes();
    }

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

      // skip "age", "enthalpy", and "Href" for now: we'll take care
      // of them a little later
      if (short_name == "enthalpy" ||
          short_name == "age"      ||
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

/*! @brief Initialize 3D fields managed by IceModel. */
/*!
 * This method should go away once we isolate "age" and "energy balance" sub-models.
 */
void IceModel::restart_3d(const PIO &input_file, unsigned int last_record) {

  // read the age field if present, otherwise set to zero
  if (m_config->get_boolean("age.enabled")) {
    bool age_exists = input_file.inq_var("age");

    if (age_exists) {
      m_ice_age.read(input_file, last_record);
    } else {
      m_log->message(2,
                     "PISM WARNING: input file '%s' does not have the 'age' variable.\n"
                     "  Setting it to zero...\n",
                     input_file.inq_filename().c_str());
      m_ice_age.set(0.0);
    }
  }

  // Initialize the enthalpy field by reading from a file or by using
  // temperature and liquid water fraction, or by using temperature
  // and assuming that the ice is cold.
  init_enthalpy(input_file, false, last_record);
}


void IceModel::initialize_2d() {
  throw RuntimeError("cannot initialize IceModel without an input file");
}


void IceModel::initialize_3d() {
  throw RuntimeError("cannot initialize IceModel without an input file");
}


void IceModel::reset_cumulative_fluxes() {
  // 2D
  m_climatic_mass_balance_cumulative.set(0.0);
  m_grounded_basal_flux_2D_cumulative.set(0.0);
  m_floating_basal_flux_2D_cumulative.set(0.0);
  m_nonneg_flux_2D_cumulative.set(0.0);
  // scalar
  m_cumulative_fluxes = FluxCounters();
}


void IceModel::initialize_cumulative_fluxes(const PIO &input_file) {
  // 2D
  {
    m_climatic_mass_balance_cumulative.regrid(input_file,  OPTIONAL, 0.0);
    m_grounded_basal_flux_2D_cumulative.regrid(input_file, OPTIONAL, 0.0);
    m_floating_basal_flux_2D_cumulative.regrid(input_file, OPTIONAL, 0.0);
    m_nonneg_flux_2D_cumulative.regrid(input_file,         OPTIONAL, 0.0);
  }

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
    sliding = new ZeroSliding(m_grid, EC);
  } else if (model == "prescribed_sliding" || model == "prescribed_sliding+sia") {
    sliding = new PrescribedSliding(m_grid, EC);
  } else if (model == "ssa" || model == "ssa+sia") {
    std::string method = m_config->get_string("stress_balance.ssa.method");

    if (method == "fem") {
      sliding = new SSAFEM(m_grid, EC);
    } else if (method == "fd") {
      sliding = new SSAFD(m_grid, EC);
    } else {
      throw RuntimeError::formatted("invalid ssa method: %s", method.c_str());
    }

  } else {
    throw RuntimeError::formatted("invalid stress balance model: %s", model.c_str());
  }

  SSB_Modifier *modifier = NULL;
  if (model == "none" || model == "ssa" || model == "prescribed_sliding") {
    modifier = new ConstantInColumn(m_grid, EC);
  } else if (model == "prescribed_sliding+sia" || model == "ssa+sia" || model == "sia") {
    modifier = new SIAFD(m_grid, EC);
  } else {
    throw RuntimeError::formatted("invalid stress balance model: %s", model.c_str());
  }

  // ~StressBalance() will de-allocate sliding and modifier.
  m_stress_balance = new StressBalance(m_grid, sliding, modifier);
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
  }
}

//! \brief Decide which bedrock thermal unit to use.
void IceModel::allocate_bedrock_thermal_unit() {

  if (m_btu != NULL) {
    return;
  }

  m_log->message(2, "# Allocating a bedrock thermal layer model...\n");

  m_btu = energy::BedThermalUnit::FromOptions(m_grid, m_ctx);
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
    throw RuntimeError::formatted("unknown value for configuration string 'hydrology.model':\n"
                                  "has value '%s'", hydrology_model.c_str());
  }
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
  if (m_config->get_boolean("energy.enabled") == true) {
    if (not m_config->get_boolean("energy.temperature_based")) {
      m_log->message(2,
                 "* Using the enthalpy-based energy balance model...\n");
    } else {
      m_log->message(2,
                 "* Using the temperature-based energy balance model...\n");
    }
  }

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
  }

  if (m_ocean == NULL) {
    m_log->message(2,
             "# Allocating an ocean model or coupler...\n");

    m_ocean = new ocean::InitializationHelper(m_grid, po.create());
  }
}


//! Allocates work vectors.
void IceModel::allocate_internal_objects() {
  const unsigned int WIDE_STENCIL = m_config->get_double("grid.max_stencil_width");

  // various internal quantities
  // 2d work vectors
  for (int j = 0; j < m_n_work2d; j++) {
    char namestr[30];
    snprintf(namestr, sizeof(namestr), "work_vector_%d", j);
    m_work2d[j].create(m_grid, namestr, WITH_GHOSTS, WIDE_STENCIL);
  }

  // 3d work vectors
  m_work3d.create(m_grid,"work_vector_3d",WITHOUT_GHOSTS);
  m_work3d.set_attrs("internal",
                    "e.g. new values of temperature or age or enthalpy during time step",
                    "", "");
}

/** 
 * Process the "mapping" variable and the "proj4" global attribute in the input file.
 */
void IceModel::get_projection_info(const PIO &input_file) {
  std::string proj4_string = input_file.get_att_text("PISM_GLOBAL", "proj4");
  if (not proj4_string.empty()) {
    m_output_global_attributes.set_string("proj4", proj4_string);
  }

  bool input_has_mapping = input_file.inq_var(m_mapping.get_name());
  if (input_has_mapping) {
    // Note: read_attributes clears attributes before reading
    io::read_attributes(input_file, m_mapping.get_name(), m_mapping);
    m_mapping.report_to_stdout(*m_log, 4);
  }

  std::string::size_type pos = proj4_string.find("+init=epsg:");
  if (pos == std::string::npos) {
    // if the PROJ.4 string does not contain "+init=epsg:", we're done
    return;
  }

  VariableMetadata epsg_mapping = epsg_to_cf(m_ctx->unit_system(), proj4_string);

  // Check that the EPSG code matches the projection information stored in the "mapping" variable
  // *unless* this variable has no attributes, in which case initialize it.
  if (input_has_mapping and
      ((not m_mapping.get_all_strings().empty()) or (not m_mapping.get_all_doubles().empty()))) {
    // Check if the "mapping" variable in the input file matches the EPSG code.
    // Check strings.
    VariableMetadata::StringAttrs strings = epsg_mapping.get_all_strings();
    VariableMetadata::StringAttrs::const_iterator j;
    for (j = strings.begin(); j != strings.end(); ++j) {
      if (not m_mapping.has_attribute(j->first)) {
        throw RuntimeError::formatted("input file '%s' has inconsistent metadata:\n"
                                      "%s requires %s = \"%s\",\n"
                                      "but the mapping variable has no %s.",
                                      input_file.inq_filename().c_str(),
                                      proj4_string.c_str(),
                                      j->first.c_str(), j->second.c_str(),
                                      j->first.c_str());
      }

      if (not (m_mapping.get_string(j->first) == j->second)) {
        throw RuntimeError::formatted("input file '%s' has inconsistent metadata:\n"
                                      "%s requires %s = \"%s\",\n"
                                      "but the mapping variable has %s = \"%s\".",
                                      input_file.inq_filename().c_str(),
                                      proj4_string.c_str(),
                                      j->first.c_str(), j->second.c_str(),
                                      j->first.c_str(),
                                      m_mapping.get_string(j->first).c_str());
      }
    }

    // Check doubles
    VariableMetadata::DoubleAttrs doubles = epsg_mapping.get_all_doubles();
    VariableMetadata::DoubleAttrs::const_iterator k;
    for (k = doubles.begin(); k != doubles.end(); ++k) {
      if (not m_mapping.has_attribute(k->first)) {
        throw RuntimeError::formatted("input file '%s' has inconsistent metadata:\n"
                                      "%s requires %s = %f,\n"
                                      "but the mapping variable has no %s.",
                                      input_file.inq_filename().c_str(),
                                      proj4_string.c_str(),
                                      k->first.c_str(), k->second[0],
                                      k->first.c_str());
      }

      if (fabs(m_mapping.get_double(k->first) - k->second[0]) > 1e-12) {
        throw RuntimeError::formatted("input file '%s' has inconsistent metadata:\n"
                                      "%s requires %s = %f,\n"
                                      "but the mapping variable has %s = %f.",
                                      input_file.inq_filename().c_str(),
                                      proj4_string.c_str(),
                                      k->first.c_str(), k->second[0],
                                      k->first.c_str(),
                                      m_mapping.get_double(k->first));
      }
    }
  } else {
    // Set "mapping" using the EPSG code.
    m_mapping = epsg_mapping;
  }
}

//! Miscellaneous initialization tasks plus tasks that need the fields that can come from regridding.
void IceModel::misc_setup() {

  m_log->message(3, "Finishing initialization...\n");

  // Check if we are initializing from a PISM output file:
  options::String input_file("-i", "Specifies a PISM input file");
  bool bootstrap = options::Bool("-bootstrap", "enable bootstrapping heuristics");

  if (input_file.is_set()) {
    PIO nc(m_grid->com, "guess_mode");

    nc.open(input_file, PISM_READONLY);

    std::string source = nc.get_att_text("PISM_GLOBAL", "source");
    nc.close();

    if (not bootstrap) {
      // If it's missing, print a warning
      if (source.empty()) {
        m_log->message(1,
                       "PISM WARNING: file '%s' does not have the 'source' global attribute.\n"
                       "     If '%s' is a PISM output file, please run the following to get rid of this warning:\n"
                       "     ncatted -a source,global,c,c,PISM %s\n",
                       input_file->c_str(), input_file->c_str(), input_file->c_str());
      } else if (source.find("PISM") == std::string::npos) {
        // If the 'source' attribute does not contain the string "PISM", then print
        // a message and stop:
        m_log->message(1,
                       "PISM WARNING: '%s' does not seem to be a PISM output file.\n"
                       "     If it is, please make sure that the 'source' global attribute contains the string \"PISM\".\n",
                       input_file->c_str());
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
      Mx_long = m_grid->Mx(),
      My_long = m_grid->My(),
      Mz_long = m_grid->Mz();
    std::string output_format = m_config->get_string("output.format");
    if (Mx_long * My_long * Mz_long * sizeof(double) > two_to_thirty_two - 4 and
        (output_format == "netcdf3" or output_format == "pnetcdf")) {
      throw RuntimeError::formatted("The computational grid is too big to fit in a NetCDF-3 file.\n"
                                    "Each 3D variable requires %lu Mb.\n"
                                    "Please use '-o_format quilt' or re-build PISM with parallel NetCDF-4 or HDF5\n"
                                    "and use '-o_format netcdf4_parallel' or '-o_format hdf5' to proceed.",
                                    Mx_long * My_long * Mz_long * sizeof(double) / (1024 * 1024));
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
    throw RuntimeError("output formats netcdf4_parallel, quilt, and hdf5 require -o_order yxz.");
  }

  // a report on whether PISM-PIK modifications of IceModel are in use
  {
    std::vector<std::string> pik_methods;
    if (m_config->get_boolean("geometry.part_grid.enabled")) {
      pik_methods.push_back("part_grid");
    }
    if (m_config->get_boolean("geometry.part_grid.redistribute_residual_volume")) {
      pik_methods.push_back("part_redist");
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

  std::istringstream arg(m_config->get_string("calving.methods"));
  std::string method_name;
  std::set<std::string> methods;

    while (getline(arg, method_name, ',')) {
      methods.insert(method_name);
    }

  if (methods.find("ocean_kill") != methods.end()) {

    if (m_ocean_kill_calving == NULL) {
      m_ocean_kill_calving = new calving::OceanKill(m_grid);
    }

    m_ocean_kill_calving->init();
    methods.erase("ocean_kill");
  }

  if (methods.find("thickness_calving") != methods.end()) {

    if (m_thickness_threshold_calving == NULL) {
      m_thickness_threshold_calving = new calving::CalvingAtThickness(m_grid);
    }

    m_thickness_threshold_calving->init();
    methods.erase("thickness_calving");
  }


  if (methods.find("eigen_calving") != methods.end()) {

    if (m_eigen_calving == NULL) {
      m_eigen_calving = new calving::EigenCalving(m_grid, m_stress_balance);
    }

    m_eigen_calving->init();
    methods.erase("eigen_calving");
  }

  if (methods.find("vonmises_calving") != methods.end()) {

    if (m_vonmises_calving == NULL) {
      m_vonmises_calving = new calving::vonMisesCalving(m_grid, m_stress_balance);
    }

    m_vonmises_calving->init();
    methods.erase("vonmises_calving");
  }

  if (methods.find("float_kill") != methods.end()) {
    if (m_float_kill_calving == NULL) {
      m_float_kill_calving = new calving::FloatKill(m_grid);
    }

    m_float_kill_calving->init();
    methods.erase("float_kill");
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
    return;
  }

  if (model == "iso") {
    m_beddef = new bed::PBPointwiseIsostasy(m_grid);
    return;
  }

  if (model == "lc") {
    m_beddef = new bed::PBLingleClark(m_grid);
    return;
  }
}

} // end of namespace pism
