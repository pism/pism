// Copyright (C) 2004-2021 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <algorithm>
#include <petscsys.h>

#include "pism/pism_config.hh"

#include "IceModel.hh"

#include "pism/basalstrength/YieldStress.hh"
#include "pism/basalstrength/basal_resistance.hh"
#include "pism/frontretreat/util/IcebergRemover.hh"
#include "pism/frontretreat/calving/CalvingAtThickness.hh"
#include "pism/frontretreat/calving/EigenCalving.hh"
#include "pism/frontretreat/calving/FloatKill.hh"
#include "pism/frontretreat/calving/HayhurstCalving.hh"
#include "pism/frontretreat/calving/vonMisesCalving.hh"
#include "pism/energy/BedThermalUnit.hh"
#include "pism/hydrology/Hydrology.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Diagnostic.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/pism_options.hh"
#include "pism/coupler/SeaLevel.hh"
#include "pism/coupler/OceanModel.hh"
#include "pism/coupler/SurfaceModel.hh"
#include "pism/earth/BedDef.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/util/pism_signal.h"
#include "pism/util/Vars.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/age/AgeModel.hh"
#include "pism/energy/EnergyModel.hh"
#include "pism/util/io/File.hh"
#include "pism/util/iceModelVec2T.hh"
#include "pism/fracturedensity/FractureDensity.hh"
#include "pism/coupler/util/options.hh" // ForcingOptions
#include "pism/coupler/util/ScalarForcing.hh"

#include <mpi.h>
extern "C"{
#include "cdipio.h"
#include "cdi.h"
#include "yaxt.h"
}

namespace pism {

IceModel::IceModel(IceGrid::Ptr grid, std::shared_ptr<Context> context)
  : m_grid(grid),
    m_opened(false),
    m_config(context->config()),
    m_ctx(context),
    m_sys(context->unit_system()),
    m_log(context->log()),
    m_time(context->time()),
    m_wide_stencil(m_config->get_number("grid.max_stencil_width")),
    m_output_global_attributes("PISM_GLOBAL", m_sys),
    m_run_stats("run_stats", m_sys),
    m_geometry(m_grid),
    m_new_bed_elevation(true),
    m_basal_yield_stress(m_grid, "tauc", WITH_GHOSTS, m_wide_stencil),
    m_basal_melt_rate(m_grid, "bmelt", WITHOUT_GHOSTS),
    m_bedtoptemp(m_grid, "bedtoptemp", WITHOUT_GHOSTS),
    m_ssa_dirichlet_bc_mask(m_grid, "bc_mask", WITH_GHOSTS, m_wide_stencil),
    m_ssa_dirichlet_bc_values(m_grid, "_ssa_bc", WITH_GHOSTS, m_wide_stencil), // u_ssa_bc and v_ssa_bc
    m_thickness_change(grid),
    m_ts_times(new std::vector<double>()),
    m_extra_bounds("time_bounds", m_sys),
    m_timestamp("timestamp", m_sys) {

  // time-independent info
  {
    m_run_stats.set_string("source", std::string("PISM ") + pism::revision);
    m_run_stats.set_string("long_name", "Run statistics");
  }

  m_extra_bounds.set_string("units", m_time->units_string());

  m_timestamp.set_string("units", "hours");
  m_timestamp.set_string("long_name", "wall-clock time since the beginning of the run");

  pism_signal = 0;
  signal(SIGTERM, pism_signal_handler);
  signal(SIGUSR1, pism_signal_handler);
  signal(SIGUSR2, pism_signal_handler);

  m_surface = nullptr;
  m_ocean   = nullptr;
  m_beddef  = nullptr;

  m_btu = nullptr;
  m_energy_model = nullptr;

  m_output_global_attributes.set_string("Conventions", "CF-1.6");
  m_output_global_attributes.set_string("source", pism::version());

  // Do not save snapshots by default:
  m_save_snapshots = false;
  // Do not save time-series by default:
  m_save_extra     = false;

  m_fracture = nullptr;

  reset_counters();

  // allocate temporary storage
  {
    // 2d work vectors
    for (int j = 0; j < m_n_work2d; j++) {
      std::shared_ptr<IceModelVec2S> ptr(new IceModelVec2S(m_grid,
                                                           pism::printf("work_vector_%d", j),
                                                           WITH_GHOSTS, m_wide_stencil));
      m_work2d.push_back(ptr);
    }
  }

  auto surface_input_file = m_config->get_string("hydrology.surface_input.file");
  if (not surface_input_file.empty()) {
    ForcingOptions surface_input(*m_ctx, "hydrology.surface_input");
    int buffer_size = m_config->get_number("input.forcing.buffer_size");
    int evaluations_per_year = m_config->get_number("input.forcing.evaluations_per_year");

    File file(m_grid->com, surface_input.filename, PISM_NETCDF3, PISM_READONLY);

    m_surface_input_for_hydrology = IceModelVec2T::ForcingField(m_grid,
                                                                file,
                                                                "water_input_rate",
                                                                "", // no standard name
                                                                buffer_size,
                                                                evaluations_per_year,
                                                                surface_input.period);
    m_surface_input_for_hydrology->set_attrs("diagnostic",
                                             "water input rate for the subglacial hydrology model",
                                             "kg m-2 s-1", "kg m-2 year-1", "", 0);
    m_surface_input_for_hydrology->metadata().set_number("valid_min", 0.0);
  }
}

double IceModel::dt() const {
  return m_dt;
}

void IceModel::reset_counters() {
  dt_TempAge       = 0.0;
  m_dt             = 0.0;
  m_skip_countdown = 0;

  m_timestep_hit_multiples_last_time = m_time->current();
}


IceModel::~IceModel() {

  delete m_beddef;

  delete m_btu;
  delete m_energy_model;
}


//! Allocate all IceModelVecs defined in IceModel.
/*!
  This procedure allocates the memory used to store model state, diagnostic and
  work vectors and sets metadata.

  Default values should not be set here; please use set_vars_from_options().

  All the memory allocated here is freed by IceModelVecs' destructors.
*/
void IceModel::allocate_storage() {

  // FIXME: this should do for now, but we should pass a const reference to Geometry to sub-models
  // as a function argument.
  m_grid->variables().add(m_geometry.ice_surface_elevation);
  m_grid->variables().add(m_geometry.ice_thickness);
  m_grid->variables().add(m_geometry.cell_type);
  m_grid->variables().add(m_geometry.sea_level_elevation);
  m_grid->variables().add(m_geometry.longitude);
  m_grid->variables().add(m_geometry.latitude);

  if (m_config->get_flag("geometry.grounded_cell_fraction")) {
    m_grid->variables().add(m_geometry.cell_grounded_fraction);
  }

  if (m_config->get_flag("geometry.part_grid.enabled")) {
    m_grid->variables().add(m_geometry.ice_area_specific_volume);
  }

  // yield stress for basal till (plastic or pseudo-plastic model)
  {
    // PROPOSED standard_name = land_ice_basal_material_yield_stress
    m_basal_yield_stress.set_attrs("diagnostic",
                                 "yield stress for basal till (plastic or pseudo-plastic model)",
                                   "Pa", "Pa", "", 0);
    m_grid->variables().add(m_basal_yield_stress);
  }

  {
    m_bedtoptemp.set_attrs("diagnostic",
                           "temperature at the top surface of the bedrock thermal layer",
                           "Kelvin", "Kelvin", "", 0);
  }

  // basal melt rate
  m_basal_melt_rate.set_attrs("internal",
                              "ice basal melt rate from energy conservation and subshelf melt, in ice thickness per time",
                              "m s-1", "m year-1", "land_ice_basal_melt_rate", 0);
  m_basal_melt_rate.metadata().set_string("comment", "positive basal melt rate corresponds to ice loss");
  m_grid->variables().add(m_basal_melt_rate);

  // SSA Dirichlet B.C. locations and values
  //
  // The mask m_ssa_dirichlet_bc_mask is also used to prescribe locations of ice thickness Dirichlet
  // B.C. (FIXME)
  {
    m_ssa_dirichlet_bc_mask.set_attrs("model_state", "Dirichlet boundary mask",
                                      "", "", "", 0);
    m_ssa_dirichlet_bc_mask.metadata().set_numbers("flag_values", {0, 1});
    m_ssa_dirichlet_bc_mask.metadata().set_string("flag_meanings", "no_data bc_condition");
    m_ssa_dirichlet_bc_mask.metadata().set_output_type(PISM_INT);
    m_ssa_dirichlet_bc_mask.set_time_independent(true);

    // FIXME: this is used by the inverse modeling code. Do NOT get
    // this field from m_grid->variables() elsewhere in the code!
    m_grid->variables().add(m_ssa_dirichlet_bc_mask);

    m_ssa_dirichlet_bc_mask.set(0.0);
  }
  // SSA Dirichlet B.C. values
  {
    double fill_value = units::convert(m_sys, m_config->get_number("output.fill_value"),
                                       "m year-1", "m second-1");
    double valid_range = units::convert(m_sys, 1e6, "m year-1", "m second-1");
    // vel_bc
    m_ssa_dirichlet_bc_values.set_attrs("model_state",
                                        "X-component of the SSA velocity boundary conditions",
                                        "m s-1", "m year-1", "", 0);
    m_ssa_dirichlet_bc_values.set_attrs("model_state",
                                        "Y-component of the SSA velocity boundary conditions",
                                        "m s-1", "m year-1", "", 1);
    for (int j = 0; j < 2; ++j) {
      m_ssa_dirichlet_bc_values.metadata(j).set_numbers("valid_range", {-valid_range, valid_range});
      m_ssa_dirichlet_bc_values.metadata(j).set_number("_FillValue", fill_value);
    }

    // FIXME: this is used by the inverse modeling code. Do NOT get
    // this field from m_grid->variables() elsewhere in the code!
    m_grid->variables().add(m_ssa_dirichlet_bc_values);
  }

  // Add some variables to the list of "model state" fields.
  m_model_state.insert(&m_ssa_dirichlet_bc_mask);
  m_model_state.insert(&m_ssa_dirichlet_bc_values);

  m_model_state.insert(&m_geometry.latitude);
  m_model_state.insert(&m_geometry.longitude);
  m_model_state.insert(&m_geometry.ice_thickness);
  m_model_state.insert(&m_geometry.ice_area_specific_volume);
}

//! Update the surface elevation and the flow-type mask when the geometry has changed.
/*!
  This procedure should be called whenever necessary to maintain consistency of geometry.

  For instance, it should be called when either ice thickness or bed elevation change.
  In particular we always want \f$h = H + b\f$ to apply at grounded points, and, on the
  other hand, we want the mask to reflect that the ice is floating if the flotation
  criterion applies at a point.

  If `flag == REMOVE_ICEBERGS`, also calls the code which removes icebergs, to avoid
  stress balance solver problems caused by ice that is not attached to the grounded ice
  sheet.
*/
void IceModel::enforce_consistency_of_geometry(ConsistencyFlag flag) {

  m_geometry.bed_elevation.copy_from(m_beddef->bed_elevation());
  m_geometry.sea_level_elevation.copy_from(m_sea_level->elevation());

  if (m_iceberg_remover and flag == REMOVE_ICEBERGS) {
    // The iceberg remover has to use the same mask as the stress balance code, hence the
    // stress-balance-related threshold here.
    m_geometry.ensure_consistency(m_config->get_number("stress_balance.ice_free_thickness_standard"));

    m_iceberg_remover->update(m_ssa_dirichlet_bc_mask,
                              m_geometry.cell_type,
                              m_geometry.ice_thickness);
    // The call above modifies ice thickness and updates the mask accordingly, but we re-compute the
    // mask (we need to use a different threshold).
  }

  // This will ensure that ice area specific volume is zero if ice thickness is greater
  // than zero, then compute new surface elevation and mask.
  m_geometry.ensure_consistency(m_config->get_number("geometry.ice_free_thickness_standard"));

  if (flag == REMOVE_ICEBERGS) {
    // clean up partially-filled cells that are not next to ice
    IceModelVec::AccessList list{&m_geometry.ice_area_specific_volume,
                                 &m_geometry.cell_type};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_geometry.ice_area_specific_volume(i, j) > 0.0 and
          not m_geometry.cell_type.next_to_ice(i, j)) {
        m_geometry.ice_area_specific_volume(i, j) = 0.0;
      }
    }
  }
}

stressbalance::Inputs IceModel::stress_balance_inputs() {
  stressbalance::Inputs result;
  if (m_config->get_flag("geometry.update.use_basal_melt_rate")) {
    result.basal_melt_rate = &m_basal_melt_rate;
  }

  result.basal_yield_stress = &m_basal_yield_stress;
  result.geometry           = &m_geometry;
  result.new_bed_elevation  = m_new_bed_elevation;
  result.enthalpy           = &m_energy_model->enthalpy();
  result.age                = m_age_model ? &m_age_model->age() : nullptr;

  result.water_column_pressure = &m_ocean->average_water_column_pressure();

  if (m_config->get_flag("stress_balance.ssa.dirichlet_bc")) {
    result.bc_mask   = &m_ssa_dirichlet_bc_mask;
    result.bc_values = &m_ssa_dirichlet_bc_values;
  }

  if (m_config->get_flag("fracture_density.enabled")) {
    result.fracture_density = &m_fracture->density();
  }

  return result;
}

energy::Inputs IceModel::energy_model_inputs() {
  energy::Inputs result;

  result.basal_frictional_heating = &m_stress_balance->basal_frictional_heating();
  result.basal_heat_flux          = &m_btu->flux_through_top_surface(); // bedrock thermal layer
  result.cell_type                = &m_geometry.cell_type;              // geometry
  result.ice_thickness            = &m_geometry.ice_thickness;          // geometry
  result.shelf_base_temp          = &m_ocean->shelf_base_temperature(); // ocean model
  result.till_water_thickness     = &m_subglacial_hydrology->till_water_thickness();
  result.surface_liquid_fraction  = &m_surface->liquid_water_fraction(); // surface model
  result.surface_temp             = &m_surface->temperature();           // surface model

  result.volumetric_heating_rate  = &m_stress_balance->volumetric_strain_heating();
  result.u3                       = &m_stress_balance->velocity_u();
  result.v3                       = &m_stress_balance->velocity_v();
  result.w3                       = &m_stress_balance->velocity_w();

  result.check();             // make sure all data members were set

  return result;
}

YieldStressInputs IceModel::yield_stress_inputs() {
  YieldStressInputs result;

  result.geometry                   = &m_geometry;
  result.till_water_thickness       = &m_subglacial_hydrology->till_water_thickness();
  result.subglacial_water_thickness = &m_subglacial_hydrology->subglacial_water_thickness();

  return result;
}

//! The contents of the main PISM time-step.
/*!
During the time-step we perform the following actions:
 */
void IceModel::step(bool do_mass_continuity,
                    bool do_skip) {

  const Profiling &profiling = m_ctx->profiling();

  double current_time = m_time->current();

  //! \li call pre_step_hook() to let derived classes do more
  pre_step_hook();

  //! \li update the velocity field; in some cases the whole three-dimensional
  //! field is updated and in some cases just the vertically-averaged
  //! horizontal velocity is updated

  // always "update" ice velocity (possibly trivially); only update
  // SSA and only update velocities at depth if suggested by temp and age
  // stability criterion; note *lots* of communication is avoided by skipping
  // SSA (and temp/age)

  const bool updateAtDepth  = (m_skip_countdown == 0);

  // Combine basal melt rate in grounded (computed during the energy
  // step) and floating (provided by an ocean model) areas.
  //
  // Basal melt rate may be used by a stress balance model to compute vertical velocity of
  // ice.
  {
    combine_basal_melt_rate(m_geometry,
                            m_ocean->shelf_base_mass_flux(),
                            m_energy_model->basal_melt_rate(),
                            m_basal_melt_rate);
  }

  try {
    profiling.begin("stress_balance");
    m_stress_balance->update(stress_balance_inputs(), updateAtDepth);
    profiling.end("stress_balance");
  } catch (RuntimeError &e) {
    std::string output_file = m_config->get_string("output.file_name");

    if (output_file.empty()) {
      m_log->message(2, "WARNING: output.file_name is empty. Using unnamed.nc instead.\n");
      output_file = "unnamed.nc";
    }

    std::string o_file = filename_add_suffix(output_file,
                                             "_stressbalance_failed", "");
    File file(m_grid->com, o_file,
              string_to_backend(m_config->get_string("output.format")),
              PISM_READWRITE_MOVE,
              m_ctx->pio_iosys_id());

    update_run_stats();
    write_metadata(file, WRITE_MAPPING, PREPEND_HISTORY);

    save_variables(file, INCLUDE_MODEL_STATE, output_variables("medium"),
                   m_time->current());

    e.add_context("performing a time step. (Note: Model state was saved to '%s'.)",
                  o_file.c_str());
    throw;
  }


  m_stdout_flags += m_stress_balance->stdout_report();

  m_stdout_flags += (updateAtDepth ? "v" : "V");

  //! \li determine the time step according to a variety of stability criteria
  auto dt_info = max_timestep(m_skip_countdown);
  m_dt                       = dt_info.dt;
  m_adaptive_timestep_reason = dt_info.reason;
  m_skip_countdown           = dt_info.skip_counter;

  //! \li update the yield stress for the plastic till model (if appropriate)
  if (m_basal_yield_stress_model) {
    profiling.begin("basal_yield_stress");
    m_basal_yield_stress_model->update(yield_stress_inputs(), current_time, m_dt);
    profiling.end("basal_yield_stress");
    m_basal_yield_stress.copy_from(m_basal_yield_stress_model->basal_material_yield_stress());
    m_stdout_flags += "y";
  } else {
    m_stdout_flags += "$";
  }

  dt_TempAge += m_dt;

  //! \li update the age of the ice (if appropriate)
  if (m_age_model and updateAtDepth) {
    AgeModelInputs inputs;
    inputs.ice_thickness = &m_geometry.ice_thickness;
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
  //!  energy_step()
  if (updateAtDepth) { // do the energy step
    profiling.begin("energy");
    energy_step();
    profiling.end("energy");
    m_stdout_flags += "E";
  } else {
    m_stdout_flags += "$";
  }

  //! \li update the fracture density field; see update_fracture_density()
  if (m_config->get_flag("fracture_density.enabled")) {
    profiling.begin("fracture_density");
    update_fracture_density();
    profiling.end("fracture_density");
  }

  //! \li update the thickness of the ice according to the mass conservation model and calving
  //! parameterizations

  // FIXME: thickness B.C. mask should be separate
  IceModelVec2Int &thickness_bc_mask = m_ssa_dirichlet_bc_mask;

  if (do_mass_continuity) {
    profiling.begin("mass_transport");
    {
      // Note that there are three adaptive time-stepping criteria. Two of them (using max.
      // diffusion and 2D CFL) are limiting the mass-continuity time-step and the third (3D
      // CFL) limits the energy and age time-steps.

      // The mass-continuity time-step is usually smaller, and the skipping mechanism lets us
      // do several mass-continuity steps for each energy step.

      // When -no_mass is set, mass-continuity-related time-step restrictions are disabled,
      // making "skipping" unnecessary.

      // This is why the following two lines appear here and are executed only if
      // do_mass_continuity is true.
      if (do_skip and m_skip_countdown > 0) {
        m_skip_countdown--;
      }

      m_geometry_evolution->flow_step(m_geometry,
                                      m_dt,
                                      m_stress_balance->advective_velocity(),
                                      m_stress_balance->diffusive_flux(),
                                      thickness_bc_mask);

      m_geometry_evolution->apply_flux_divergence(m_geometry);

      enforce_consistency_of_geometry(DONT_REMOVE_ICEBERGS);
    }
    profiling.end("mass_transport");

    // calving, frontal melt, and discharge accounting
    profiling.begin("front_retreat");
    front_retreat_step();
    profiling.end("front_retreat");

    m_stdout_flags += "h";
  } else {
    m_stdout_flags += "$";
  }

  profiling.begin("sea_level");
  m_sea_level->update(m_geometry, current_time, m_dt);
  profiling.end("sea_level");

  profiling.begin("ocean");
  m_ocean->update(m_geometry, current_time, m_dt);
  profiling.end("ocean");

  // The sea level elevation might have changed, so we need to update the mask, etc. Note
  // that THIS MAY PRODUCE ICEBERGS, but we assume that the surface model does not care.
  enforce_consistency_of_geometry(DONT_REMOVE_ICEBERGS);

  //! \li Update surface and ocean models.
  profiling.begin("surface");
  m_surface->update(m_geometry, current_time, m_dt);
  profiling.end("surface");


  if (do_mass_continuity) {
    // compute and apply effective surface and basal mass balance

    m_geometry_evolution->source_term_step(m_geometry, m_dt,
                                           thickness_bc_mask,
                                           m_surface->mass_flux(),
                                           m_basal_melt_rate);
    m_geometry_evolution->apply_mass_fluxes(m_geometry);

    // add removed icebergs to discharge due to calving
    {
      IceModelVec2S
        &old_H    = *m_work2d[0],
        &old_Href = *m_work2d[1];

      {
        old_H.copy_from(m_geometry.ice_thickness);
        old_Href.copy_from(m_geometry.ice_area_specific_volume);
      }

      // the last call has to remove icebergs
      enforce_consistency_of_geometry(REMOVE_ICEBERGS);

      bool add_values = true;
      compute_geometry_change(m_geometry.ice_thickness,
                              m_geometry.ice_area_specific_volume,
                              old_H, old_Href,
                              add_values,
                              m_thickness_change.calving);
    }
  }

  //! \li update the state variables in the subglacial hydrology model (typically
  //!  water thickness and sometimes pressure)
  profiling.begin("basal_hydrology");
  hydrology_step();
  profiling.end("basal_hydrology");

  //! \li compute the bed deformation, which depends on current thickness, bed elevation,
  //! and sea level
  if (m_beddef) {
    int topg_state_counter = m_beddef->bed_elevation().state_counter();

    profiling.begin("bed_deformation");
    m_beddef->update(m_geometry.ice_thickness,
                     m_geometry.sea_level_elevation,
                     current_time, m_dt);
    profiling.end("bed_deformation");

    if (m_beddef->bed_elevation().state_counter() != topg_state_counter) {
      m_new_bed_elevation = true;
    } else {
      m_new_bed_elevation = false;
    }
  } else {
    m_new_bed_elevation = false;
  }

  if (m_new_bed_elevation) {
    enforce_consistency_of_geometry(DONT_REMOVE_ICEBERGS);
    m_stdout_flags += "b";
  } else {
    m_stdout_flags += " ";
  }

  //! \li call post_step_hook() to let derived classes do more
  post_step_hook();

  // Done with the step; now adopt the new time.
  m_time->step(m_dt);

  if (updateAtDepth) {
    t_TempAge  = m_time->current();
    dt_TempAge = 0.0;
  }

  // Check if the ice thickness exceeded the height of the computational box and stop if it did.
  const bool thickness_too_high = check_maximum_ice_thickness(m_geometry.ice_thickness);

  if (thickness_too_high) {
    std::string output_file = m_config->get_string("output.file_name");

    if (output_file.empty()) {
      m_log->message(2, "WARNING: output.file_name is empty. Using unnamed.nc instead.");
      output_file = "unnamed.nc";
    }

    std::string o_file = filename_add_suffix(output_file,
                                             "_max_thickness", "");
    File file(m_grid->com,
              o_file,
              string_to_backend(m_config->get_string("output.format")),
              PISM_READWRITE_MOVE,
              m_ctx->pio_iosys_id());

    update_run_stats();
    write_metadata(file, WRITE_MAPPING, PREPEND_HISTORY);

    save_variables(file, INCLUDE_MODEL_STATE, output_variables("small"),
                   m_time->current());

    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "Ice thickness exceeds the height of the computational box (%7.4f m).\n"
                                  "The model state was saved to '%s'. To continue this simulation,\n"
                                  "run with\n"
                                  "-i %s -bootstrap -regrid_file %s -allow_extrapolation -Lz N [other options]\n"
                                  "where N > %7.4f.",
                                  m_grid->Lz(), o_file.c_str(), o_file.c_str(), o_file.c_str(), m_grid->Lz());
  }

  // end the flag line
  m_stdout_flags += " " + m_adaptive_timestep_reason;
}

/*!
 * Note: don't forget to update IceRegionalModel::hydrology_step() if necessary.
 */
void IceModel::hydrology_step() {
  hydrology::Inputs inputs;

  IceModelVec2S &sliding_speed = *m_work2d[0];
  sliding_speed.set_to_magnitude(m_stress_balance->advective_velocity());

  inputs.no_model_mask      = nullptr;
  inputs.geometry           = &m_geometry;
  inputs.surface_input_rate = nullptr;
  inputs.basal_melt_rate    = &m_basal_melt_rate;
  inputs.ice_sliding_speed  = &sliding_speed;

  if (m_surface_input_for_hydrology) {
    m_surface_input_for_hydrology->update(m_time->current(), m_dt);
    m_surface_input_for_hydrology->average(m_time->current(), m_dt);
    inputs.surface_input_rate = m_surface_input_for_hydrology.get();
  } else if (m_config->get_flag("hydrology.surface_input_from_runoff")) {
    // convert [kg m-2] to [kg m-2 s-1]
    IceModelVec2S &surface_input_rate = *m_work2d[1];
    surface_input_rate.copy_from(m_surface->runoff());
    surface_input_rate.scale(1.0 / m_dt);
    inputs.surface_input_rate = &surface_input_rate;
  }

  m_subglacial_hydrology->update(m_time->current(), m_dt, inputs);
}

//! Virtual.  Does nothing in `IceModel`.  Derived classes can do more computation in each time step.
void IceModel::pre_step_hook() {
  // empty
}

//! Virtual.  Does nothing in `IceModel`.  Derived classes can do more computation in each time step.
void IceModel::post_step_hook() {
  // empty
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
 */
void IceModel::run() {
  const Profiling &profiling = m_ctx->profiling();

  bool do_mass_conserve = m_config->get_flag("geometry.update.enabled");
  bool do_energy = m_config->get_flag("energy.enabled");
  bool do_skip = m_config->get_flag("time_stepping.skip.enabled");

  int stepcount = m_config->get_flag("time_stepping.count_steps") ? 0 : -1;

  // de-allocate diagnostics that are not needed
  prune_diagnostics();

  // Enforce consistency *and* remove icebergs. During time-stepping we remove icebergs at
  // the end of the time step, so we need to ensure that ice geometry is "OK" before the
  // first step.
  enforce_consistency_of_geometry(REMOVE_ICEBERGS);

  // Update spatially-variable diagnostics at the beginning of the run.
  write_extras();

  // Update scalar time series to remember the state at the beginning of the run.
  // This is needed to compute rates of change of the ice mass, volume, etc.
  {
    const double time = m_time->current();
    for (auto d : m_ts_diagnostics) {
      d.second->update(time, time);
    }
  }

  m_log->message(2, "running forward ...\n");

  m_stdout_flags.erase(); // clear it out
  print_summary_line(true, do_energy, 0.0, 0.0, 0.0, 0.0, 0.0);
  print_summary(do_energy);  // report starting state

  t_TempAge = m_time->current();
  dt_TempAge = 0.0;

  // main loop for time evolution
  // IceModel::step calls Time::step(dt), ensuring that this while loop
  // will terminate
  profiling.stage_begin("time-stepping loop");
  while (m_time->current() < m_time->end()) {

    m_stdout_flags.erase();  // clear it out

    step(do_mass_conserve, do_skip);

    update_diagnostics(m_dt);

    // report a summary for major steps or the last one
    bool updateAtDepth = m_skip_countdown == 0;
    bool tempAgeStep   = updateAtDepth and (m_age_model or do_energy);

    double one_second = 1.0;
    const bool show_step = tempAgeStep or fabs(m_time->current() - m_time->end()) < one_second;
    print_summary(show_step);

    // update viewers before writing extras because writing extras resets diagnostics
    update_viewers();

    // writing these fields here ensures that we do it after the last time-step
    profiling.begin("io");
    if (not m_opened) {
      open_files();
      m_opened = true;
    }
    write_snapshot();
    write_extras();
    write_backup();
    if (m_sthwritten) {
      pioWriteTimestep();
      m_sthwritten = false;
    }
    profiling.end("io");

    if (stepcount >= 0) {
      stepcount++;
    }
    if (process_signals() != 0) {
      break;
    }
  } // end of the time-stepping loop

  profiling.stage_end("time-stepping loop");

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
  // Get the start time in seconds and ensure that it is consistent
  // across all processors.
  m_start_time = GlobalMax(m_grid->com, get_time());

  const Profiling &profiling = m_ctx->profiling();

  profiling.begin("initialization");

  //! The IceModel initialization sequence is this:

  //! 1) Initialize model time:
  time_setup();

  //! 2) Process the options:
  process_options();

  //! 3) Memory allocation:
  allocate_storage();

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

  profiling.end("initialization");
}

const Geometry& IceModel::geometry() const {
  return m_geometry;
}

const GeometryEvolution& IceModel::geometry_evolution() const {
  return *m_geometry_evolution;
}

const stressbalance::StressBalance* IceModel::stress_balance() const {
  return this->m_stress_balance.get();
}

const ocean::OceanModel* IceModel::ocean_model() const {
  return m_ocean.get();
}

const bed::BedDef* IceModel::bed_model() const {
  return m_beddef;
}

const energy::BedThermalUnit* IceModel::bedrock_thermal_model() const {
  return m_btu;
}

const energy::EnergyModel* IceModel::energy_balance_model() const {
  return m_energy_model;
}

/*!
 * Return thickness change due to calving (over the last time step).
 */
const IceModelVec2S& IceModel::calving() const {
  return m_thickness_change.calving;
}

/*!
 * Return thickness change due to frontal melt (over the last time step).
 */
const IceModelVec2S& IceModel::frontal_melt() const {
  return m_thickness_change.frontal_melt;
}

/*!
 * Return thickness change due to forced retreat (over the last time step).
 */
const IceModelVec2S& IceModel::forced_retreat() const {
  return m_thickness_change.forced_retreat;
}

void warn_about_missing(const Logger &log,
                        const std::set<std::string> &vars,
                        const std::string &type,
                        const std::set<std::string> &available,
                        bool stop) {
  std::vector<std::string> missing;
  for (auto v : vars) {
    if (available.find(v) == available.end()) {
      missing.push_back(v);
    }
  }

  if (not missing.empty()) {
    size_t N = missing.size();
    const char
      *ending = N > 1 ? "s" : "",
      *verb   = N > 1 ? "are" : "is";
    if (stop) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "%s variable%s %s %s not available!\n"
                                    "Available variables:\n- %s",
                                    type.c_str(),
                                    ending,
                                    join(missing, ",").c_str(),
                                    verb,
                                    set_join(available, ",\n- ").c_str());
    } else {
      log.message(2,
                  "\nWARNING: %s variable%s %s %s not available!\n\n",
                  type.c_str(),
                  ending,
                  join(missing, ",").c_str(),
                  verb);
    }
  }
}

/*!
 * De-allocate diagnostics that were not requested.
 *
 * Checks viewers, -extra_vars, -backup, -save_vars, and regular output.
 *
 * FIXME: I need to make sure that these reporting mechanisms are active. It is possible that
 * variables are on a list, but that list is not actually used.
 */
void IceModel::prune_diagnostics() {

  // get the list of available diagnostics
  std::set<std::string> available;
  for (auto d : m_diagnostics) {
    available.insert(d.first);
  }

  auto m_extra_stop = m_config->get_flag("output.extra.stop_missing");
  warn_about_missing(*m_log, m_output_vars,   "output",     available, false);
  warn_about_missing(*m_log, m_snapshot_vars, "snapshot",   available, false);
  warn_about_missing(*m_log, m_backup_vars,   "backup",     available, false);
  warn_about_missing(*m_log, m_extra_vars,    "diagnostic", available, m_extra_stop);

  // get the list of requested diagnostics
  auto requested = set_split(m_config->get_string("output.runtime.viewer.variables"), ',');
  requested = combine(requested, m_output_vars);
  requested = combine(requested, m_snapshot_vars);
  requested = combine(requested, m_extra_vars);
  requested = combine(requested, m_backup_vars);

  // de-allocate diagnostics that were not requested
  for (auto v : available) {
    if (requested.find(v) == requested.end()) {
      m_diagnostics.erase(v);
    }
  }

  // scalar time series
  std::vector<std::string> missing;
  if (not m_ts_filename.empty() and m_ts_vars.empty()) {
    // use all diagnostics
  } else {
    TSDiagnosticList diagnostics;
    for (auto v : m_ts_vars) {
      if (m_ts_diagnostics.find(v) != m_ts_diagnostics.end()) {
        diagnostics[v] = m_ts_diagnostics[v];
      } else {
        missing.push_back(v);
      }
    }
    // replace m_ts_diagnostics with requested diagnostics, de-allocating the rest
    m_ts_diagnostics = diagnostics;
  }

  if (not missing.empty()) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "requested scalar diagnostics %s are not available",
                                  join(missing, ",").c_str());
  }
}

/*!
 * Update diagnostics.
 *
 * This usually involves accumulating data needed to computed time-averaged quantities.
 *
 * Call this after prune_diagnostics() to avoid unnecessary work.
 */
void IceModel::update_diagnostics(double dt) {
  for (auto d : m_diagnostics) {
    d.second->update(dt);
  }

  const double time = m_time->current();
  for (auto d : m_ts_diagnostics) {
    d.second->update(time - dt, time);
  }
}

/*!
 * Reset accumulators in diagnostics that compute time-averaged quantities.
 */
void IceModel::reset_diagnostics() {
  for (auto d : m_diagnostics) {
    d.second->reset();
  }
}

IceModel::ThicknessChanges::ThicknessChanges(IceGrid::ConstPtr grid)
  : calving(grid, "thickness_change_due_to_calving", WITHOUT_GHOSTS),
    frontal_melt(grid, "thickness_change_due_to_frontal_melt", WITHOUT_GHOSTS),
    forced_retreat(grid, "thickness_change_due_to_forced_retreat", WITHOUT_GHOSTS) {
  // empty
}

} // end of namespace pism
