// Copyright (C) 2004-2016 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <cmath>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <petscsys.h>

#include "iceModel.hh"

#include "base/basalstrength/PISMYieldStress.hh"
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
#include "base/util/IceGrid.hh"
#include "base/util/Mask.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/PISMDiagnostic.hh"
#include "base/util/error_handling.hh"
#include "base/util/pism_options.hh"
#include "coupler/PISMOcean.hh"
#include "coupler/PISMSurface.hh"
#include "earth/PISMBedDef.hh"
#include "enthalpyConverter.hh"
#include "pism_signal.h"
#include "base/util/PISMVars.hh"
#include "base/util/Profiling.hh"
#include "base/util/pism_utilities.hh"
#include "base/age/AgeModel.hh"
#include "base/energy/EnergyModel.hh"

namespace pism {

IceModel::FluxCounters::FluxCounters() {
  grounded_basal = 0.0;
  nonneg_rule    = 0.0;
  sub_shelf      = 0.0;
  surface        = 0.0;
  sum_divQ_SIA   = 0.0;
  sum_divQ_SSA   = 0.0;
  Href_to_H      = 0.0;
  H_to_Href      = 0.0;
}

FractureFields::FractureFields(IceGrid::ConstPtr grid) {
  Config::ConstPtr config = grid->ctx()->config();

  const unsigned int WIDE_STENCIL = config->get_double("grid.max_stencil_width");

  density.create(grid, "fracture_density", WITH_GHOSTS, WIDE_STENCIL);
  density.set_attrs("model_state", "fracture density in ice shelf", "", "");
  density.metadata().set_double("valid_max", 1.0);
  density.metadata().set_double("valid_min", 0.0);

  growth_rate.create(grid, "fracture_growth_rate", WITH_GHOSTS, WIDE_STENCIL);
  growth_rate.set_attrs("model_state", "fracture growth rate", "second-1", "");
  growth_rate.metadata().set_double("valid_min", 0.0);

  healing_rate.create(grid, "fracture_healing_rate", WITH_GHOSTS, WIDE_STENCIL);
  healing_rate.set_attrs("model_state", "fracture healing rate", "second-1", "");

  flow_enhancement.create(grid, "fracture_flow_enhancement", WITH_GHOSTS, WIDE_STENCIL);
  flow_enhancement.set_attrs("model_state", "fracture-induced flow enhancement", "", "");

  age.create(grid, "fracture_age", WITH_GHOSTS, WIDE_STENCIL);
  age.set_attrs("model_state", "age since fracturing", "years", "");

  toughness.create(grid, "fracture_toughness", WITH_GHOSTS, WIDE_STENCIL);
  toughness.set_attrs("model_state", "fracture toughness", "Pa", "");

  strain_rates.create(grid, "strain_rates", WITH_GHOSTS,
                      2, // stencil width, has to match or exceed the "offset" in eigenCalving
                      2);

  strain_rates.metadata(0).set_name("eigen1");
  strain_rates.set_attrs("internal",
                         "major principal component of horizontal strain-rate",
                         "second-1", "", 0);

  strain_rates.metadata(1).set_name("eigen2");
  strain_rates.set_attrs("internal",
                         "minor principal component of horizontal strain-rate",
                         "second-1", "", 1);

  deviatoric_stresses.create(grid, "sigma", WITH_GHOSTS,
                             2, // stencil width
                             3);

  deviatoric_stresses.metadata(0).set_name("sigma_xx");
  deviatoric_stresses.set_attrs("internal",
                                "deviatoric stress in x direction",
                                "Pa", "", 0);

  deviatoric_stresses.metadata(1).set_name("sigma_yy");
  deviatoric_stresses.set_attrs("internal",
                                "deviatoric stress in y direction",
                                "Pa", "", 1);

  deviatoric_stresses.metadata(2).set_name("sigma_xy");
  deviatoric_stresses.set_attrs("internal",
                                "deviatoric shear stress",
                                "Pa", "", 2);
}

IceModel::FluxFields::FluxFields(IceGrid::ConstPtr grid) {
  {
    climatic_mass_balance.create(grid, "climatic_mass_balance_cumulative",
                                 WITHOUT_GHOSTS);
    climatic_mass_balance.set_attrs("diagnostic",
                                    "cumulative surface mass balance",
                                    "kg m-2", "");
  }
  {
    basal_grounded.create(grid, "grounded_basal_flux_cumulative", WITHOUT_GHOSTS);
    basal_grounded.set_attrs("diagnostic",
                             "cumulative basal flux into the ice "
                             "in grounded areas (positive means ice gain)",
                             "kg m-2", "");
    basal_grounded.set_time_independent(false);
    basal_grounded.metadata().set_string("glaciological_units", "Gt m-2");
    basal_grounded.write_in_glaciological_units = true;
  }
  {
    basal_floating.create(grid, "floating_basal_flux_cumulative", WITHOUT_GHOSTS);
    basal_floating.set_attrs("diagnostic",
                             "cumulative basal flux into the ice "
                             "in floating areas (positive means ice gain)",
                             "kg m-2", "");
    basal_floating.set_time_independent(false);
    basal_floating.metadata().set_string("glaciological_units", "Gt m-2");
    basal_floating.write_in_glaciological_units = true;
  }
  {
    nonneg.create(grid, "nonneg_flux_cumulative", WITHOUT_GHOSTS);
    nonneg.set_attrs("diagnostic",
                     "cumulative nonnegative rule flux (positive means ice gain)",
                     "kg m-2", "");
    nonneg.set_time_independent(false);
    nonneg.metadata().set_string("glaciological_units", "Gt m-2");
    nonneg.write_in_glaciological_units = true;
  }
  {
    discharge.create(grid, "discharge_flux_cumulative", WITHOUT_GHOSTS);
    discharge.set_attrs("diagnostic",
                        "cumulative discharge (calving) flux (positive means ice loss)",
                        "kg", "");
    discharge.set_time_independent(false);
    discharge.metadata().set_string("glaciological_units", "Gt");
    discharge.write_in_glaciological_units = true;
  }
}

void IceModel::FluxFields::reset() {
  climatic_mass_balance.set(0.0);
  basal_grounded.set(0.0);
  basal_floating.set(0.0);
  nonneg.set(0.0);
  discharge.set(0.0);
}

void IceModel::FluxFields::regrid(const PIO &input_file) {
  climatic_mass_balance.regrid(input_file, OPTIONAL, 0.0);
  basal_grounded.regrid(input_file,        OPTIONAL, 0.0);
  basal_floating.regrid(input_file,        OPTIONAL, 0.0);
  nonneg.regrid(input_file,                OPTIONAL, 0.0);
  discharge.regrid(input_file,             OPTIONAL, 0.0);
}

IceModel::IceModel(IceGrid::Ptr g, Context::Ptr context)
  : m_grid(g),
    m_config(context->config()),
    m_ctx(context),
    m_sys(context->unit_system()),
    m_log(context->log()),
    m_time(context->time()),
    m_output_global_attributes("PISM_GLOBAL", m_sys),
    m_run_stats("run_stats", m_sys),
    m_cumulative_flux_fields(m_grid),
    m_extra_bounds("time_bounds", m_config->get_string("time.dimension_name"), m_sys),
    m_timestamp("timestamp", m_config->get_string("time.dimension_name"), m_sys) {

  m_extra_bounds.set_string("units", m_time->units_string());

  m_timestamp.set_string("units", "hours");
  m_timestamp.set_string("long_name", "wall-clock time since the beginning of the run");

  pism_signal = 0;
  signal(SIGTERM, pism_signal_handler);
  signal(SIGUSR1, pism_signal_handler);
  signal(SIGUSR2, pism_signal_handler);

  m_subglacial_hydrology = NULL;
  m_basal_yield_stress_model = NULL;

  m_stress_balance = NULL;

  m_surface = NULL;
  m_ocean   = NULL;
  m_beddef  = NULL;

  m_age_model = NULL;
  m_btu = NULL;
  m_energy_model = NULL;

  m_iceberg_remover             = NULL;
  m_ocean_kill_calving          = NULL;
  m_float_kill_calving          = NULL;
  m_thickness_threshold_calving = NULL;
  m_eigen_calving               = NULL;
  m_vonmises_calving            = NULL;
  m_frontal_melt                = NULL;

  m_output_global_attributes.set_string("Conventions", "CF-1.5");
  m_output_global_attributes.set_string("source", pism::version());

  // Do not save snapshots by default:
  m_save_snapshots = false;
  // Do not save time-series by default:
  m_save_ts        = false;
  m_save_extra     = false;

  m_fracture = NULL;

  reset_counters();

  // allocate temporary storage
  {
    const unsigned int WIDE_STENCIL = m_config->get_double("grid.max_stencil_width");

    // various internal quantities
    // 2d work vectors
    for (int j = 0; j < m_n_work2d; j++) {
      char namestr[30];
      snprintf(namestr, sizeof(namestr), "work_vector_%d", j);
      m_work2d[j].create(m_grid, namestr, WITH_GHOSTS, WIDE_STENCIL);
    }
  }
}

IceModel::FluxCounters IceModel::cumulative_fluxes() const {
  return m_cumulative_fluxes;
}

const IceModel::FluxFields& IceModel::cumulative_fluxes_2d() const {
  return m_cumulative_flux_fields;
}

const IceModelVec2S& IceModel::flux_divergence() const {
  return m_flux_divergence;
}

double IceModel::dt() const {
  return m_dt;
}

void IceModel::reset_counters() {
  dt_TempAge   = 0.0;

  m_dt              = 0.0;
  m_skip_countdown   = 0;

  m_timestep_hit_multiples_last_time = m_time->current();
}


IceModel::~IceModel() {

  delete m_age_model;

  delete m_stress_balance;

  delete m_ocean;

  delete m_surface;

  delete m_beddef;

  delete m_subglacial_hydrology;
  delete m_basal_yield_stress_model;
  delete m_btu;
  delete m_energy_model;

  delete m_fracture;

  delete m_iceberg_remover;
  delete m_ocean_kill_calving;
  delete m_float_kill_calving;
  delete m_thickness_threshold_calving;
  delete m_eigen_calving;
  delete m_vonmises_calving;
  delete m_frontal_melt;
}


//! Allocate all IceModelVecs defined in IceModel.
/*!
  This procedure allocates the memory used to store model state, diagnostic and
  work vectors and sets metadata.

  Default values should not be set here; please use set_vars_from_options().

  All the memory allocated here is freed by IceModelVecs' destructors.
*/
void IceModel::createVecs() {

  const unsigned int WIDE_STENCIL = m_config->get_double("grid.max_stencil_width");

  // The following code creates (and documents -- to some extent) the
  // variables. The main (and only) principle here is using standard names from
  // the CF conventions; see
  // http://cf-pcmdi.llnl.gov/documents/cf-standard-names

  // ice upper surface elevation
  m_ice_surface_elevation.create(m_grid, "usurf", WITH_GHOSTS, WIDE_STENCIL);
  m_ice_surface_elevation.set_attrs("diagnostic", "ice upper surface elevation",
                                  "m", "surface_altitude");
  m_grid->variables().add(m_ice_surface_elevation);

  // land ice thickness
  {
    m_ice_thickness.create(m_grid, "thk", WITH_GHOSTS, WIDE_STENCIL);
    m_ice_thickness.set_attrs("model_state", "land ice thickness",
                              "m", "land_ice_thickness");
    m_ice_thickness.metadata().set_double("valid_min", 0.0);
    m_grid->variables().add(m_ice_thickness);
  }

  {
    m_flux_divergence.create(m_grid, "flux_divergence", WITHOUT_GHOSTS);
    m_flux_divergence.set_attrs("diagnostic",
                                "flux divergence",
                                "m s-1", "");
    m_flux_divergence.metadata().set_string("glaciological_units", "m year-1");
    m_flux_divergence.write_in_glaciological_units = true;
  }

  if (m_config->get_boolean("geometry.grounded_cell_fraction")) {
    m_gl_mask.create(m_grid, "gl_mask", WITHOUT_GHOSTS);
    m_gl_mask.set_attrs("internal",
                      "fractional grounded/floating mask (floating=0, grounded=1)",
                      "", "");
    m_grid->variables().add(m_gl_mask);
  }

  // grounded_dragging_floating integer mask
  m_cell_type.create(m_grid, "mask", WITH_GHOSTS, WIDE_STENCIL);
  m_cell_type.set_attrs("diagnostic", "ice-type (ice-free/grounded/floating/ocean) integer mask",
                  "", "");
  std::vector<double> mask_values(4);
  mask_values[0] = MASK_ICE_FREE_BEDROCK;
  mask_values[1] = MASK_GROUNDED;
  mask_values[2] = MASK_FLOATING;
  mask_values[3] = MASK_ICE_FREE_OCEAN;
  m_cell_type.metadata().set_doubles("flag_values", mask_values);
  m_cell_type.metadata().set_string("flag_meanings",
                              "ice_free_bedrock grounded_ice floating_ice ice_free_ocean");
  m_grid->variables().add(m_cell_type);

  // yield stress for basal till (plastic or pseudo-plastic model)
  {
    m_basal_yield_stress.create(m_grid, "tauc", WITH_GHOSTS, WIDE_STENCIL);
    // PROPOSED standard_name = land_ice_basal_material_yield_stress
    m_basal_yield_stress.set_attrs("diagnostic",
                                 "yield stress for basal till (plastic or pseudo-plastic model)",
                                 "Pa", "");
    m_grid->variables().add(m_basal_yield_stress);
  }

  // basal melt rate
  m_basal_melt_rate.create(m_grid, "bmelt", WITHOUT_GHOSTS);
  m_basal_melt_rate.set_attrs("internal",
                              "ice basal melt rate from energy conservation and subshelf melt, in ice thickness per time",
                              "m s-1", "land_ice_basal_melt_rate");
  m_basal_melt_rate.metadata().set_string("glaciological_units", "m year-1");
  m_basal_melt_rate.write_in_glaciological_units = true;
  m_basal_melt_rate.metadata().set_string("comment", "positive basal melt rate corresponds to ice loss");
  m_grid->variables().add(m_basal_melt_rate);

  // longitude
  m_longitude.create(m_grid, "lon", WITH_GHOSTS);
  m_longitude.set_attrs("mapping", "longitude", "degree_east", "longitude");
  m_longitude.set_time_independent(true);
  m_longitude.metadata().set_string("coordinates", "");
  m_longitude.metadata().set_string("grid_mapping", "");
  m_longitude.metadata().set_double("valid_min", -180.0);
  m_longitude.metadata().set_double("valid_max",  180.0);
  m_grid->variables().add(m_longitude);

  // latitude
  m_latitude.create(m_grid, "lat", WITH_GHOSTS); // has ghosts so that we can compute cell areas
  m_latitude.set_attrs("mapping", "latitude", "degree_north", "latitude");
  m_latitude.set_time_independent(true);
  m_latitude.metadata().set_string("coordinates", "");
  m_latitude.metadata().set_string("grid_mapping", "");
  m_latitude.metadata().set_double("valid_min", -90.0);
  m_latitude.metadata().set_double("valid_max",  90.0);
  m_grid->variables().add(m_latitude);

  if (m_config->get_boolean("geometry.part_grid.enabled")) {
    // Href
    m_Href.create(m_grid, "Href", WITH_GHOSTS);
    m_Href.set_attrs("model_state", "temporary ice thickness at calving front boundary",
                    "m", "");
    m_grid->variables().add(m_Href);
  }

  if (m_config->get_boolean("stress_balance.ssa.dirichlet_bc")) {
    // bc_locations
    m_ssa_dirichlet_bc_mask.create(m_grid, "bc_mask", WITH_GHOSTS, WIDE_STENCIL);
    m_ssa_dirichlet_bc_mask.set_attrs("model_state", "Dirichlet boundary mask",
                      "", "");
    std::vector<double> bc_mask_values(2);
    bc_mask_values[0] = 0;
    bc_mask_values[1] = 1;
    m_ssa_dirichlet_bc_mask.metadata().set_doubles("flag_values", bc_mask_values);
    m_ssa_dirichlet_bc_mask.metadata().set_string("flag_meanings", "no_data bc_condition");
    m_grid->variables().add(m_ssa_dirichlet_bc_mask);

    double fill_value = units::convert(m_sys, m_config->get_double("output.fill_value"),
                                       "m year-1", "m second-1");
    // vel_bc
    m_ssa_dirichlet_bc_values.create(m_grid, "_ssa_bc", WITH_GHOSTS, WIDE_STENCIL); // u_ssa_bc and v_ssa_bc
    m_ssa_dirichlet_bc_values.set_attrs("model_state",
                     "X-component of the SSA velocity boundary conditions",
                     "m s-1", "", 0);
    m_ssa_dirichlet_bc_values.set_attrs("model_state",
                     "Y-component of the SSA velocity boundary conditions",
                     "m s-1", "", 1);
    for (int j = 0; j < 2; ++j) {
      m_ssa_dirichlet_bc_values.metadata(j).set_string("glaciological_units", "m year-1");
      m_ssa_dirichlet_bc_values.metadata(j).set_double("valid_min",  units::convert(m_sys, -1e6, "m year-1", "m second-1"));
      m_ssa_dirichlet_bc_values.metadata(j).set_double("valid_max",  units::convert(m_sys,  1e6, "m year-1", "m second-1"));
      m_ssa_dirichlet_bc_values.metadata(j).set_double("_FillValue", fill_value);
    }
    //just for diagnostics...
    m_grid->variables().add(m_ssa_dirichlet_bc_values);
  }


  // cell areas
  m_cell_area.create(m_grid, "cell_area", WITHOUT_GHOSTS);
  m_cell_area.set_attrs("diagnostic", "cell areas", "m2", "");
  m_cell_area.metadata().set_string("comment",
                                  "values are equal to dx*dy "
                                  "if projection parameters are not available; "
                                  "otherwise WGS84 ellipsoid is used");
  m_cell_area.set_time_independent(true);
  m_cell_area.metadata().set_string("glaciological_units", "km2");
  m_cell_area.write_in_glaciological_units = true;
  m_grid->variables().add(m_cell_area);

  if (m_config->get_boolean("fracture_density.enabled")) {
    m_fracture = new FractureFields(m_grid);

    m_grid->variables().add(m_fracture->toughness);
    m_grid->variables().add(m_fracture->age);
    m_grid->variables().add(m_fracture->flow_enhancement);
    m_grid->variables().add(m_fracture->healing_rate);
    m_grid->variables().add(m_fracture->growth_rate);
    m_grid->variables().add(m_fracture->density);
  }
}

//! The contents of the main PISM time-step.
/*!
During the time-step we perform the following actions:
 */
void IceModel::step(bool do_mass_continuity,
                    bool do_skip) {

  const Profiling &profiling = m_ctx->profiling();

  double current_time = m_time->current();

  //! \li call additionalAtStartTimestep() to let derived classes do more
  additionalAtStartTimestep();  // might set maxdt_temporary

  //! \li update the velocity field; in some cases the whole three-dimensional
  //! field is updated and in some cases just the vertically-averaged
  //! horizontal velocity is updated

  // always "update" ice velocity (possibly trivially); only update
  // SSA and only update velocities at depth if suggested by temp and age
  // stability criterion; note *lots* of communication is avoided by skipping
  // SSA (and temp/age)

  const bool updateAtDepth  = (m_skip_countdown == 0);

  //! \li update the yield stress for the plastic till model (if appropriate)
  if (updateAtDepth and m_basal_yield_stress_model) {
    profiling.begin("basal yield stress");
    m_basal_yield_stress_model->update();
    profiling.end("basal yield stress");
    m_basal_yield_stress.copy_from(m_basal_yield_stress_model->basal_material_yield_stress());
    m_stdout_flags += "y";
  } else {
    m_stdout_flags += "$";
  }

  // Update the fractional grounded/floating mask (used by the SSA
  // stress balance and the energy code)
  if (m_config->get_boolean("geometry.grounded_cell_fraction")) {
    enforce_consistency_of_geometry(); // update h and mask
    update_grounded_cell_fraction();
  }

  IceModelVec2S &melange_back_pressure = m_work2d[0];

  m_ocean->melange_back_pressure_fraction(melange_back_pressure);

  try {
    profiling.begin("stress balance");
    m_stress_balance->update(not updateAtDepth,
                             m_ocean->sea_level_elevation(),
                             melange_back_pressure);
    profiling.end("stress balance");
  } catch (RuntimeError &e) {
    options::String output_file("-o", "output file name",
                                "output.nc", options::DONT_ALLOW_EMPTY);

    std::string o_file = pism_filename_add_suffix(output_file,
                                                  "_stressbalance_failed", "");
    dumpToFile(o_file);

    e.add_context("performing a time step. (Note: Model state was saved to '%s'.)",
                  o_file.c_str());
    throw;
  }


  m_stdout_flags += m_stress_balance->stdout_report();

  m_stdout_flags += (updateAtDepth ? "v" : "V");

  //! \li determine the time step according to a variety of stability criteria;
  //!  see determineTimeStep()
  max_timestep(m_dt, m_skip_countdown);

  //! \li Update surface and ocean models.
  profiling.begin("surface");
  m_surface->update(current_time, m_dt);
  profiling.end("surface");

  profiling.begin("ocean");
  m_ocean->update(current_time, m_dt);
  profiling.end("ocean");

  dt_TempAge += m_dt;
  // IceModel::dt,dtTempAge are now set correctly according to
  // mass-continuity-eqn-diffusivity criteria, horizontal CFL criteria, and
  // other criteria from derived class additionalAtStartTimestep(), and from
  // "-skip" mechanism

  //! \li update the age of the ice (if appropriate)
  if (m_age_model != NULL and updateAtDepth) {
    AgeModelInputs inputs;
    inputs.ice_thickness = &m_ice_thickness;
    inputs.u3            = &m_stress_balance->velocity_u();
    inputs.v3            = &m_stress_balance->velocity_v();
    inputs.w3            = &m_stress_balance->velocity_w();

    profiling.begin("age");
    m_age_model->update(current_time, dt_TempAge, inputs);
    profiling.end("age");
    m_stdout_flags += "a";
  } else {
    m_stdout_flags += "$";
  }

  //! \li update the enthalpy (or temperature) field according to the conservation of
  //!  energy model based (especially) on the new velocity field; see
  //!  energyStep()
  if (updateAtDepth) { // do the energy step
    profiling.begin("energy");
    energyStep();
    profiling.end("energy");
    m_stdout_flags += "E";
  } else {
    m_stdout_flags += "$";
  }

  // Combine basal melt rate in grounded (computed during the energy
  // step) and floating (provided by an ocean model) areas.
  combine_basal_melt_rate();

  //! \li update the state variables in the subglacial hydrology model (typically
  //!  water thickness and sometimes pressure)
  profiling.begin("basal hydrology");
  m_subglacial_hydrology->update(current_time, m_dt);
  profiling.end("basal hydrology");

  //! \li update the fracture density field; see calculateFractureDensity()
  if (m_config->get_boolean("fracture_density.enabled")) {
    profiling.begin("fracture density");
    calculateFractureDensity();
    profiling.end("fracture density");
  }

  //! \li update the thickness of the ice according to the mass conservation
  //!  model; see massContExplicitStep()
  if (do_mass_continuity) {
    profiling.begin("mass transport");
    massContExplicitStep();
    enforce_consistency_of_geometry(); // update h and mask
    profiling.end("mass transport");

    // Note that there are three adaptive time-stepping criteria. Two of them
    // (using max. diffusion and 2D CFL) are limiting the mass-continuity
    // time-step and the third (3D CFL) limits the energy and age time-steps.

    // The mass-continuity time-step is usually smaller, and the skipping
    // mechanism lets us do several mass-continuity steps for each energy step.

    // When -no_mass is set, mass-continuity-related time-step restrictions are
    // disabled, making "skipping" unnecessary.

    // This is why the following two lines appear here and are executed only
    // if do_mass_continuity is true.
    if (do_skip and m_skip_countdown > 0) {
      m_skip_countdown--;
    }
    m_stdout_flags += "h";
  } else {
    m_stdout_flags += "$";
  }

  profiling.begin("calving");
  do_calving();
  profiling.end("calving");

  //! \li compute the bed deformation, which only depends on current thickness
  //! and bed elevation
  if (m_beddef) {
    const IceModelVec2S &bed_topography = m_beddef->bed_elevation();
    int topg_state_counter = bed_topography.get_state_counter();

    profiling.begin("bed deformation");
    m_beddef->update(current_time, m_dt);
    profiling.end("bed deformation");

    if (bed_topography.get_state_counter() != topg_state_counter) {
      m_stdout_flags += "b";
      enforce_consistency_of_geometry();
    } else {
      m_stdout_flags += " ";
    }
  }

  //! \li call additionalAtEndTimestep() to let derived classes do more
  additionalAtEndTimestep();

  // Done with the step; now adopt the new time.
  m_time->step(m_dt);

  if (updateAtDepth) {
    t_TempAge = m_time->current();
    dt_TempAge = 0.0;
  }

  // end the flag line
  m_stdout_flags += " " + m_adaptive_timestep_reason;
}


/**
 * Run the time-stepping loop from the current model time to `time`.
 *
 * This should be called by the coupler controlling PISM when it is
 * running alongside a GCM.
 *
 * @param run_end model time (in seconds) to run to
 *
 * @return 0 on success
 */
void IceModel::run_to(double run_end) {

  m_time->set_end(run_end);

  run();
}


/**
 * Run the time-stepping loop from the current time until the time
 * specified by the IceModel::grid::time object.
 *
 * This is the method used by PISM in the "standalone" mode.
 *
 * @return 0 on success
 */
void IceModel::run() {
  const Profiling &profiling = m_ctx->profiling();

  bool do_mass_conserve = m_config->get_boolean("geometry.update.enabled");
  bool do_energy = m_config->get_boolean("energy.enabled");
  bool do_skip = m_config->get_boolean("time_stepping.skip.enabled");

  int stepcount = m_config->get_boolean("time_stepping.count_steps") ? 0 : -1;

  enforce_consistency_of_geometry();

  // update diagnostics at the beginning of the run:
  write_timeseries();
  write_extras();

  m_log->message(2, "running forward ...\n");

  m_stdout_flags.erase(); // clear it out
  summaryPrintLine(true, do_energy, 0.0, 0.0, 0.0, 0.0, 0.0);
  m_adaptive_timestep_reason = '$'; // no reason for no timestep
  summary(do_energy);  // report starting state

  t_TempAge = m_time->current();
  dt_TempAge = 0.0;

  // main loop for time evolution
  // IceModel::step calls Time::step(dt), ensuring that this while loop
  // will terminate
  profiling.stage_begin("time-stepping loop");
  while (m_time->current() < m_time->end()) {

    m_stdout_flags.erase();  // clear it out

    step(do_mass_conserve, do_skip);

    // report a summary for major steps or the last one
    bool updateAtDepth = m_skip_countdown == 0;
    bool tempAgeStep = updateAtDepth and (do_energy or m_age_model != NULL);

    const bool show_step = tempAgeStep or m_adaptive_timestep_reason == "end of the run";
    summary(show_step);

    // writing these fields here ensures that we do it after the last time-step
    profiling.begin("I/O during run");
    write_snapshot();
    write_timeseries();
    write_extras();
    write_backup();
    profiling.end("I/O during run");

    update_viewers();

    if (stepcount >= 0) {
      stepcount++;
    }
    if (endOfTimeStepHook() != 0) {
      break;
    }
  } // end of the time-stepping loop

  profiling.stage_end("time-stepping loop");

  options::Integer pause_time("-pause", "Pause after the run, seconds", 0);
  if (pause_time > 0) {
    m_log->message(2, "pausing for %d secs ...\n", pause_time.value());
    PetscErrorCode ierr = PetscSleep(pause_time); PISM_CHK(ierr, "PetscSleep");
  }

  if (stepcount >= 0) {
    m_log->message(1,
               "count_time_steps:  run() took %d steps\n"
               "average dt = %.6f years\n",
               stepcount,
               units::convert(m_sys, m_time->end() - m_time->start(), "seconds", "years")/(double)stepcount);
  }
}

//! Manage the initialization of the IceModel object.
/*!
Please see the documenting comments of the functions called below to find
explanations of their intended uses.
 */
void IceModel::init() {
  const Profiling &profiling = m_ctx->profiling();

  profiling.begin("initialization");

  //! The IceModel initialization sequence is this:

  //! 1) Initialize model time:
  time_setup();

  //! 2) Process the options:
  setFromOptions();

  //! 3) Memory allocation:
  createVecs();

  //! 4) Allocate PISM components modeling some physical processes.
  allocate_submodels();

  //! 6) Initialize coupler models and fill the model state variables
  //! (from a PISM output file, from a bootstrapping file using some
  //! modeling choices or using formulas). Calls IceModel::regrid()
  model_state_setup();

  //! 7) Report grid parameters:
  m_grid->report_parameters();

  //! 8) Miscellaneous stuff: set up the bed deformation model, initialize the
  //! basal till model, initialize snapshots. This has to happen *after*
  //! regridding.
  misc_setup();

  //! The following flow-chart illustrates the process.
  //!
  //! \dotfile initialization-sequence.dot "IceModel initialization sequence"

  // Get the start time in seconds and ensure that it is consistent
  // across all processors.
  m_start_time = GlobalMax(m_grid->com, GetTime());

  profiling.end("initialization");
}

const stressbalance::StressBalance* IceModel::stress_balance() const {
  return this->m_stress_balance;
}

const ocean::OceanModel* IceModel::ocean_model() const {
  return this->m_ocean;
}

const bed::BedDef* IceModel::bed_model() const {
  return m_beddef;
}

const energy::BedThermalUnit* IceModel::bedrock_thermal_model() const {
  return m_btu;
}

const IceModelVec2S& IceModel::ice_thickness() const {
  return m_ice_thickness;
}

const IceModelVec2S & IceModel::cell_area() const {
  return m_cell_area;
}

const IceModelVec2CellType & IceModel::cell_type() const {
  return m_cell_type;
}

const IceModelVec2S& IceModel::ice_surface_elevation() const {
  return m_ice_surface_elevation;
}

const energy::EnergyModel* IceModel::energy_balance_model() const {
  return m_energy_model;
}

} // end of namespace pism
