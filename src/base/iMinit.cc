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

#include <petscdmda.h>
#include <cassert>
#include <algorithm>

#include "iceModel.hh"
#include "base/basalstrength/PISMConstantYieldStress.hh"
#include "base/basalstrength/PISMMohrCoulombYieldStress.hh"
#include "base/basalstrength/basal_resistance.hh"
#include "base/calving/PISMCalvingAtThickness.hh"
#include "base/calving/PISMEigenCalving.hh"
#include "base/calving/PISMFloatKill.hh"
#include "base/calving/PISMIcebergRemover.hh"
#include "base/calving/PISMOceanKill.hh"
#include "base/energy/bedrockThermalUnit.hh"
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
#include "coupler/surface/PSFactory.hh"
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
                  m_config->get_string("time_dimension_name"),
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

  reset_counters();

  // Initialize (or re-initialize) boundary models.
  init_couplers();

  // Check if we are initializing from a PISM output file:
  options::String input_file("-i", "Specifies the PISM input file");
  bool bootstrap = options::Bool("-bootstrap", "enable bootstrapping heuristics");

  if (input_file.is_set() and not bootstrap) {
    initFromFile(input_file);

    regrid(0);
    // Check consistency of geometry after initialization:
    updateSurfaceElevationAndMask();
  } else {
    set_vars_from_options();
  }

  // Initialize a bed deformation model (if needed); this should go
  // after the regrid(0) call but before other init() calls that need
  // bed elevation and uplift.
  if (m_beddef) {
    m_beddef->init();
    m_grid->variables().add(m_beddef->bed_elevation());
    m_grid->variables().add(m_beddef->uplift());
  }

  if (m_stress_balance) {
    m_stress_balance->init();

    if (m_config->get_boolean("include_bmr_in_continuity")) {
      m_stress_balance->set_basal_melt_rate(m_basal_melt_rate);
    }
  }

  if (btu) {
    bool bootstrapping_needed = false;
    btu->init(bootstrapping_needed);

    if (bootstrapping_needed) {
      // update surface and ocean models so that we can get the
      // temperature at the top of the bedrock
      m_log->message(2,
                 "getting surface B.C. from couplers...\n");
      init_step_couplers();

      get_bed_top_temp(m_bedtoptemp);

      btu->bootstrap();
    }
  }

  if (subglacial_hydrology) {
    subglacial_hydrology->init();
  }

  // basal_yield_stress_model->init() needs bwat so this must happen
  // after subglacial_hydrology->init()
  if (basal_yield_stress_model) {
    basal_yield_stress_model->init();
  }

  {
    if (input_file.is_set()) {
      m_log->message(2,
                 "* Trying to read cumulative climatic mass balance from '%s'...\n",
                 input_file->c_str());
      m_climatic_mass_balance_cumulative.regrid(input_file, OPTIONAL, 0.0);
    } else {
      m_climatic_mass_balance_cumulative.set(0.0);
    }
  }

  {
    if (input_file.is_set()) {
      m_log->message(2,
                 "* Trying to read cumulative grounded basal flux from '%s'...\n",
                 input_file->c_str());
      m_grounded_basal_flux_2D_cumulative.regrid(input_file, OPTIONAL, 0.0);
    } else {
      m_grounded_basal_flux_2D_cumulative.set(0.0);
    }
  }

  {
    if (input_file.is_set()) {
      m_log->message(2,
                 "* Trying to read cumulative floating basal flux from '%s'...\n",
                 input_file->c_str());
      m_floating_basal_flux_2D_cumulative.regrid(input_file, OPTIONAL, 0.0);
    } else {
      m_floating_basal_flux_2D_cumulative.set(0.0);
    }
  }

  {
    if (input_file.is_set()) {
      m_log->message(2,
                 "* Trying to read cumulative nonneg flux from '%s'...\n",
                 input_file->c_str());
      m_nonneg_flux_2D_cumulative.regrid(input_file, OPTIONAL, 0.0);
    } else {
      m_nonneg_flux_2D_cumulative.set(0.0);
    }
  }

  if (input_file.is_set()) {
    PIO nc(m_grid->com, "netcdf3");

    nc.open(input_file, PISM_READONLY);
    bool run_stats_exists = nc.inq_var("run_stats");
    if (run_stats_exists) {
      io::read_attributes(nc, run_stats.get_name(), run_stats);
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

  // get projection information and compute cell areas
  {
    if (input_file.is_set()) {
      PIO nc(m_grid->com, "guess_mode");
      nc.open(input_file, PISM_READONLY); // closed at the end of scope

      get_projection_info(nc);

      std::string proj4_string = m_output_global_attributes.get_string("proj4");
      if (not proj4_string.empty()) {
        m_log->message(2, "* Got projection parameters \"%s\" from \"%s\".\n",
                       proj4_string.c_str(), nc.inq_filename().c_str());
      }
    }

    compute_cell_areas();
  }

  // a report on whether PISM-PIK modifications of IceModel are in use
  std::vector<std::string> pik_methods;
  if (m_config->get_boolean("part_grid")) {
    pik_methods.push_back("part_grid");
  }
  if (m_config->get_boolean("part_redist")) {
    pik_methods.push_back("part_redist");
  }
  if (m_config->get_boolean("kill_icebergs")) {
    pik_methods.push_back("kill_icebergs");
  }

  if (not pik_methods.empty()) {
    m_log->message(2,
                   "* PISM-PIK mass/geometry methods are in use: %s\n",
                   join(pik_methods, ", ").c_str());
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

  m_log->message(3,
             "Setting initial values of model state variables...\n");

  options::String input_file("-i", "Specifies the input file");
  bool bootstrap = options::Bool("-bootstrap", "enable bootstrapping heuristics");

  if (bootstrap and input_file.is_set()) {
    bootstrapFromFile(input_file);
  } else {
    throw RuntimeError("No input file specified.");
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

  std::string model = m_config->get_string("stress_balance_model");

  ShallowStressBalance *sliding = NULL;
  if (model == "none" || model == "sia") {
    sliding = new ZeroSliding(m_grid, EC);
  } else if (model == "prescribed_sliding" || model == "prescribed_sliding+sia") {
    sliding = new PrescribedSliding(m_grid, EC);
  } else if (model == "ssa" || model == "ssa+sia") {
    std::string method = m_config->get_string("ssa_method");

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

  if (iceberg_remover != NULL) {
    return;
  }

  m_log->message(2,
             "# Allocating an iceberg remover (part of a calving model)...\n");

  if (m_config->get_boolean("kill_icebergs")) {

    // this will throw an exception on failure
    iceberg_remover = new calving::IcebergRemover(m_grid);

    // Iceberg Remover does not have a state, so it is OK to
    // initialize here.
    iceberg_remover->init();
  }
}

//! \brief Decide which bedrock thermal unit to use.
void IceModel::allocate_bedrock_thermal_unit() {

  if (btu != NULL) {
    return;
  }

  m_log->message(2,
             "# Allocating a bedrock thermal layer model...\n");

  btu = new energy::BedThermalUnit(m_grid);
}

//! \brief Decide which subglacial hydrology model to use.
void IceModel::allocate_subglacial_hydrology() {

  using namespace pism::hydrology;

  std::string hydrology_model = m_config->get_string("hydrology_model");

  if (subglacial_hydrology != NULL) { // indicates it has already been allocated
    return;
  }

  m_log->message(2,
             "# Allocating a subglacial hydrology model...\n");

  if (hydrology_model == "null") {
    subglacial_hydrology = new NullTransport(m_grid);
  } else if (hydrology_model == "routing") {
    subglacial_hydrology = new Routing(m_grid);
  } else if (hydrology_model == "distributed") {
    subglacial_hydrology = new Distributed(m_grid, m_stress_balance);
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

  m_log->message(2,
             "# Allocating a basal yield stress model...\n");

  std::string model = m_config->get_string("stress_balance_model");

  // only these two use the yield stress (so far):
  if (model == "ssa" || model == "ssa+sia") {
    std::string yield_stress_model = m_config->get_string("yield_stress_model");

    if (yield_stress_model == "constant") {
      basal_yield_stress_model = new ConstantYieldStress(m_grid);
    } else if (yield_stress_model == "mohr_coulomb") {
      basal_yield_stress_model = new MohrCoulombYieldStress(m_grid, subglacial_hydrology);
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
  if (m_config->get_boolean("do_energy") == true) {
    if (not m_config->get_boolean("do_cold_ice_methods")) {
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

    m_surface = ps.create();
    m_external_surface_model = false;

    atmosphere = pa.create();
    m_surface->attach_atmosphere_model(atmosphere);
  }

  if (m_ocean == NULL) {
    m_log->message(2,
             "# Allocating an ocean model or coupler...\n");

    m_ocean = po.create();
    m_external_ocean_model = false;
  }
}

//! Initializes atmosphere and ocean couplers.
void IceModel::init_couplers() {

  m_log->message(3,
             "Initializing boundary models...\n");

  assert(m_surface != NULL);
  m_surface->init();

  assert(m_ocean != NULL);
  m_ocean->init();
}


//! Some sub-models need fields provided by surface and ocean models
//! for initialization, so here we call update() to make sure that
//! surface and ocean models report a decent state
void IceModel::init_step_couplers() {

  const double
    now               = m_time->current(),
    one_year_from_now = m_time->increment_date(now, 1.0);

  // Take a one year long step if we can.
  MaxTimestep max_dt(one_year_from_now - now);

  assert(m_surface != NULL);
  max_dt = std::min(max_dt, m_surface->max_timestep(now));

  assert(m_ocean != NULL);
  max_dt = std::min(max_dt, m_ocean->max_timestep(now));

  // Do not take time-steps shorter than 1 second
  if (max_dt.value() < 1.0) {
    max_dt = MaxTimestep(1.0);
  }

  assert(max_dt.is_finite() == true);

  m_surface->update(now, max_dt.value());
  m_ocean->update(now, max_dt.value());
}


//! Allocates work vectors.
void IceModel::allocate_internal_objects() {
  const unsigned int WIDE_STENCIL = m_config->get_double("grid_max_stencil_width");

  // various internal quantities
  // 2d work vectors
  for (int j = 0; j < nWork2d; j++) {
    char namestr[30];
    snprintf(namestr, sizeof(namestr), "work_vector_%d", j);
    vWork2d[j].create(m_grid, namestr, WITH_GHOSTS, WIDE_STENCIL);
  }

  // 3d work vectors
  vWork3d.create(m_grid,"work_vector_3d",WITHOUT_GHOSTS);
  vWork3d.set_attrs("internal",
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

  bool input_has_mapping = input_file.inq_var(mapping.get_name());
  if (input_has_mapping) {
    // Note: read_attributes clears attributes before reading
    io::read_attributes(input_file, mapping.get_name(), mapping);
    mapping.report_to_stdout(*m_log, 4);
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
      ((not mapping.get_all_strings().empty()) or (not mapping.get_all_doubles().empty()))) {
    // Check if the "mapping" variable in the input file matches the EPSG code.
    // Check strings.
    VariableMetadata::StringAttrs strings = epsg_mapping.get_all_strings();
    VariableMetadata::StringAttrs::const_iterator j;
    for (j = strings.begin(); j != strings.end(); ++j) {
      if (not mapping.has_attribute(j->first)) {
        throw RuntimeError::formatted("input file '%s' has inconsistent metadata:\n"
                                      "%s requires %s = \"%s\",\n"
                                      "but the mapping variable has no %s.",
                                      input_file.inq_filename().c_str(),
                                      proj4_string.c_str(),
                                      j->first.c_str(), j->second.c_str(),
                                      j->first.c_str());
      }

      if (not (mapping.get_string(j->first) == j->second)) {
        throw RuntimeError::formatted("input file '%s' has inconsistent metadata:\n"
                                      "%s requires %s = \"%s\",\n"
                                      "but the mapping variable has %s = \"%s\".",
                                      input_file.inq_filename().c_str(),
                                      proj4_string.c_str(),
                                      j->first.c_str(), j->second.c_str(),
                                      j->first.c_str(),
                                      mapping.get_string(j->first).c_str());
      }
    }

    // Check doubles
    VariableMetadata::DoubleAttrs doubles = epsg_mapping.get_all_doubles();
    VariableMetadata::DoubleAttrs::const_iterator k;
    for (k = doubles.begin(); k != doubles.end(); ++k) {
      if (not mapping.has_attribute(k->first)) {
        throw RuntimeError::formatted("input file '%s' has inconsistent metadata:\n"
                                      "%s requires %s = %f,\n"
                                      "but the mapping variable has no %s.",
                                      input_file.inq_filename().c_str(),
                                      proj4_string.c_str(),
                                      k->first.c_str(), k->second[0],
                                      k->first.c_str());
      }

      if (fabs(mapping.get_double(k->first) - k->second[0]) > 1e-12) {
        throw RuntimeError::formatted("input file '%s' has inconsistent metadata:\n"
                                      "%s requires %s = %f,\n"
                                      "but the mapping variable has %s = %f.",
                                      input_file.inq_filename().c_str(),
                                      proj4_string.c_str(),
                                      k->first.c_str(), k->second[0],
                                      k->first.c_str(),
                                      mapping.get_double(k->first));
      }
    }
  } else {
    // Set "mapping" using the EPSG code.
    mapping = epsg_mapping;
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
    std::string output_format = m_config->get_string("output_format");
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

  // Quietly re-initialize couplers (they might have done one
  // time-step during initialization)
  {
    m_log->disable();
    init_couplers();
    m_log->enable();
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
  std::string o_format = m_config->get_string("output_format");
  if ((o_format == "netcdf4_parallel" || o_format == "quilt" || o_format == "hdf5") &&
      m_config->get_string("output_variable_order") != "yxz") {
    throw RuntimeError("output formats netcdf4_parallel, quilt, and hdf5 require -o_order yxz.");
  }
}

//! \brief Initialize calving mechanisms.
void IceModel::init_calving() {

  std::istringstream arg(m_config->get_string("calving_methods"));
  std::string method_name;
  std::set<std::string> methods;

    while (getline(arg, method_name, ',')) {
      methods.insert(method_name);
    }

  if (methods.find("ocean_kill") != methods.end()) {

    if (ocean_kill_calving == NULL) {
      ocean_kill_calving = new calving::OceanKill(m_grid);
    }

    ocean_kill_calving->init();
    methods.erase("ocean_kill");
  }

  if (methods.find("thickness_calving") != methods.end()) {

    if (thickness_threshold_calving == NULL) {
      thickness_threshold_calving = new calving::CalvingAtThickness(m_grid);
    }

    thickness_threshold_calving->init();
    methods.erase("thickness_calving");
  }


  if (methods.find("eigen_calving") != methods.end()) {

    if (eigen_calving == NULL) {
      eigen_calving = new calving::EigenCalving(m_grid, m_stress_balance);
    }

    eigen_calving->init();
    methods.erase("eigen_calving");
  }

  if (methods.find("float_kill") != methods.end()) {
    if (float_kill_calving == NULL) {
      float_kill_calving = new calving::FloatKill(m_grid);
    }

    float_kill_calving->init();
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
  std::string model = m_config->get_string("bed_deformation_model");

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
