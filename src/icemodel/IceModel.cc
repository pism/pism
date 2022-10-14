// Copyright (C) 2004-2023 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "pism/icemodel/IceModel.hh"

#include "pism/basalstrength/YieldStress.hh"
#include "pism/frontretreat/util/IcebergRemover.hh"
#include "pism/energy/BedThermalUnit.hh"
#include "pism/hydrology/Hydrology.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/util/Grid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Diagnostic.hh"
#include "pism/util/error_handling.hh"
#include "pism/coupler/SeaLevel.hh"
#include "pism/coupler/OceanModel.hh"
#include "pism/coupler/SurfaceModel.hh"
#include "pism/earth/BedDef.hh"
#include "pism/util/pism_signal.h"
#include "pism/util/Vars.hh"
#include "pism/util/Profiling.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/age/AgeModel.hh"
#include "pism/age/Isochrones.hh"
#include "pism/energy/EnergyModel.hh"
#include "pism/util/io/File.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/fracturedensity/FractureDensity.hh"
#include "pism/coupler/util/options.hh" // ForcingOptions
#include "pism/coupler/ocean/PyOceanModel.hh"

namespace pism {

IceModel::IceModel(std::shared_ptr<Grid> grid, const std::shared_ptr<Context> &context)
    : m_grid(grid),
      m_config(context->config()),
      m_ctx(context),
      m_sys(context->unit_system()),
      m_log(context->log()),
      m_time(context->time()),
      m_wide_stencil(static_cast<int>(m_config->get_number("grid.max_stencil_width"))),
      m_output_global_attributes("PISM_GLOBAL", m_sys),
      m_geometry(m_grid),
      m_new_bed_elevation(true),
      m_basal_yield_stress(m_grid, "tauc"),
      m_basal_melt_rate(m_grid, "bmelt"),
      m_bedtoptemp(m_grid, "bedtoptemp"),
      m_velocity_bc_mask(m_grid, "vel_bc_mask"),
      m_velocity_bc_values(m_grid, "_bc"), // u_bc and v_bc
      m_ice_thickness_bc_mask(grid, "thk_bc_mask"),
      m_step_counter(0),
      m_thickness_change(grid),
      m_ts_times(new std::vector<double>()),
      m_extra_bounds("time_bounds", m_sys),
      m_timestamp("timestamp", m_sys) {

  m_velocity_bc_mask.set_interpolation_type(NEAREST);
  m_ice_thickness_bc_mask.set_interpolation_type(NEAREST);

  m_extra_bounds["units"] = m_time->units_string();

  m_timestamp["units"] = "hours";
  m_timestamp["long_name"] = "wall-clock time since the beginning of the run";

  pism_signal = 0;
  signal(SIGTERM, pism_signal_handler);
  signal(SIGUSR1, pism_signal_handler);
  signal(SIGUSR2, pism_signal_handler);

  m_surface = nullptr;
  m_ocean   = nullptr;
  m_beddef  = nullptr;

  m_btu = nullptr;
  m_energy_model = nullptr;

  m_output_global_attributes["Conventions"] = "CF-1.6";
  m_output_global_attributes["source"] = pism::version();

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
      m_work2d.push_back(std::make_shared<array::Scalar2>(m_grid,
                                                             pism::printf("work_vector_%d", j)));
    }
  }

  auto surface_input_file = m_config->get_string("hydrology.surface_input.file");
  if (not surface_input_file.empty()) {
    ForcingOptions surface_input(*m_ctx, "hydrology.surface_input");
    int buffer_size = static_cast<int>(m_config->get_number("input.forcing.buffer_size"));

    File file(m_grid->com, surface_input.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    m_surface_input_for_hydrology =
        std::make_shared<array::Forcing>(m_grid, file, "water_input_rate",
                                         "", // no standard name
                                         buffer_size, surface_input.periodic);
    m_surface_input_for_hydrology->metadata(0)
        .long_name("water input rate for the subglacial hydrology model")
        .units("kg m-2 s-1")
        .output_units("kg m-2 year-1");
    m_surface_input_for_hydrology->metadata()["valid_min"] = { 0.0 };
  }
}

double IceModel::dt() const {
  return m_dt;
}

void IceModel::reset_counters() {
  dt_TempAge       = 0.0;
  m_dt             = 0.0;
  m_skip_countdown = 0;

  {
    int year_increment = static_cast<int>(m_config->get_number("time_stepping.hit_multiples"));

    if (year_increment > 0) {
      // start year:
      auto year = m_time->units().date(m_time->current(), m_time->calendar()).year;

      int last_multiple = year - year % year_increment;
      // correct last_multiple if 'year' is negative
      // and not a multiple of year_increment:
      last_multiple -= year_increment * static_cast<int>(year % year_increment < 0);

      units::DateTime last_date{ last_multiple, 1, 1, 0, 0, 0.0 };

      m_timestep_hit_multiples_last_time = m_time->units().time(last_date, m_time->calendar());
    } else {
      m_timestep_hit_multiples_last_time = m_time->current();
    }
  }
}


IceModel::~IceModel() {
  // empty; defined here to be able to use more forward-declared classes in IceModel.hh
}


//! Allocate all Arrays defined in IceModel.
/*!
  This procedure allocates the memory used to store model state, diagnostic and
  work vectors and sets metadata.

  Default values should not be set here; please use set_vars_from_options().

  All the memory allocated here is freed by Arrays' destructors.
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
    m_basal_yield_stress.metadata(0)
        .long_name("yield stress for basal till (plastic or pseudo-plastic model)")
        .units("Pa");
    m_grid->variables().add(m_basal_yield_stress);
  }

  {
    m_bedtoptemp.metadata(0)
        .long_name("temperature at the top surface of the bedrock thermal layer")
        .units("Kelvin");
  }

  // basal melt rate
  m_basal_melt_rate.metadata(0)
      .long_name(
          "ice basal melt rate from energy conservation and subshelf melt, in ice thickness per time")
      .units("m s-1")
      .output_units("m year-1")
      .standard_name("land_ice_basal_melt_rate");
  m_basal_melt_rate.metadata()["comment"] = "positive basal melt rate corresponds to ice loss";
  m_grid->variables().add(m_basal_melt_rate);

  // Sliding velocity (usually SSA) Dirichlet B.C. locations and values
  {
    m_velocity_bc_mask.metadata(0)
        .long_name("Mask prescribing Dirichlet boundary locations for the sliding velocity")
        .set_output_type(io::PISM_INT)
        .set_time_independent(true);
    m_velocity_bc_mask.metadata()["flag_values"]   = { 0, 1 };
    m_velocity_bc_mask.metadata()["flag_meanings"] = "no_data boundary_condition";

    m_velocity_bc_mask.set(0.0);
  }
  // SSA Dirichlet B.C. values
  {
    double fill_value       = m_config->get_number("output.fill_value");
    const double huge_value = 1e6;
    // vel_bc
    m_velocity_bc_values.metadata(0).long_name(
        "X-component of the SSA velocity boundary conditions");
    m_velocity_bc_values.metadata(1).long_name(
        "Y-component of the SSA velocity boundary conditions");
    for (int j : { 0, 1 }) {
      m_velocity_bc_values.metadata(j)["valid_range"] = { -huge_value, huge_value };
      m_velocity_bc_values.metadata(j)["_FillValue"]  = { fill_value };
      m_velocity_bc_values.metadata(j).units("m s-1");
    }
  }

  // Ice thickness BC mask
  {
    m_ice_thickness_bc_mask.metadata(0)
        .long_name("Mask specifying locations where ice thickness is held constant")
        .set_time_independent(true)
        .set_output_type(io::PISM_INT);
    m_ice_thickness_bc_mask.metadata()["flag_values"] = {0, 1};
    m_ice_thickness_bc_mask.metadata()["flag_meanings"] = "no_data boundary_condition";

    m_ice_thickness_bc_mask.set(0.0);
  }

  // Add some variables to the list of "model state" fields.
  m_model_state.insert(&m_velocity_bc_mask);
  m_model_state.insert(&m_velocity_bc_values);

  m_model_state.insert(&m_ice_thickness_bc_mask);

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

    m_iceberg_remover->update(m_ice_thickness_bc_mask,
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
    array::AccessScope list{&m_geometry.ice_area_specific_volume,
                                 &m_geometry.cell_type};

    for (auto p = m_grid->points(); p; p.next()) {
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
    result.bc_mask   = &m_velocity_bc_mask;
    result.bc_values = &m_velocity_bc_values;
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

std::string IceModel::save_state_on_error(const std::string &suffix,
                                          const std::set<std::string> &additional_variables) {
  std::string output_file = m_config->get_string("output.file");

  if (output_file.empty()) {
    m_log->message(2, "WARNING: output.file is empty. Using unnamed.nc instead.");
    output_file = "unnamed.nc";
  }

  output_file = filename_add_suffix(output_file, suffix, "");

  File file(m_grid->com,
            output_file,
            string_to_backend(m_config->get_string("output.format")),
            io::PISM_READWRITE_MOVE,
            m_ctx->pio_iosys_id());

  run_stats();

  write_metadata(file, WRITE_MAPPING, PREPEND_HISTORY);

  auto variables = output_variables("small");
  for (const auto &v : additional_variables) {
    variables.insert(v);
  }

  save_variables(file, INCLUDE_MODEL_STATE, variables, m_time->current());

  return output_file;
}

//! The contents of the main PISM time-step.
/*!
During the time-step we perform the following actions:
 */
void IceModel::step(bool do_mass_continuity,
                    bool do_skip) {

  m_step_counter++;

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
    std::string output_file = save_state_on_error("_stressbalance_failed", {});

    e.add_context("performing a time step. (Note: Model state was saved to '%s'.)",
                  output_file.c_str());
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

  if (do_mass_continuity) {
    // reset the conservation error field:
    m_geometry_evolution->reset();

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

      try {
        m_geometry_evolution->flow_step(m_geometry,
                                        m_dt,
                                        m_stress_balance->advective_velocity(),
                                        m_stress_balance->diffusive_flux(),
                                        m_ice_thickness_bc_mask);
      } catch (RuntimeError &e) {
        std::string output_file = save_state_on_error("_mass_transport_failed",
                                                      {"flux_staggered", "flux_divergence"});

        e.add_context("performing a mass transport time step (dt=%f s). (Note: Model state was saved to '%s'.)",
                      m_dt, output_file.c_str());
        throw;
      }

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
                                           m_ice_thickness_bc_mask,
                                           m_surface->mass_flux(),
                                           m_basal_melt_rate);
    m_geometry_evolution->apply_mass_fluxes(m_geometry);

    // add removed icebergs to discharge due to calving
    {
      auto &old_H    = *m_work2d[0];
      auto &old_Href = *m_work2d[1];

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

  if (m_isochrones) {
    m_isochrones->update(current_time, m_dt,
                         m_stress_balance->velocity_u(),
                         m_stress_balance->velocity_v(),
                         m_geometry.ice_thickness,
                         m_geometry_evolution->top_surface_mass_balance(),
                         m_geometry_evolution->bottom_surface_mass_balance());
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

    m_new_bed_elevation = m_beddef->bed_elevation().state_counter() != topg_state_counter;
  } else {
    m_new_bed_elevation = false;
  }

  if (m_config->get_flag("time_stepping.assume_bed_elevation_changed")) {
    m_new_bed_elevation = true;
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
  if (max(m_geometry.ice_thickness) > m_grid->Lz()) {
    auto o_file = save_state_on_error("_max_thickness", {});

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

  array::Scalar &sliding_speed = *m_work2d[0];
  compute_magnitude(m_stress_balance->advective_velocity(), sliding_speed);

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
    array::Scalar &surface_input_rate = *m_work2d[1];
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
IceModelTerminationReason IceModel::run_to(double run_end) {

  m_time->set_end(run_end);

  return run();
}


/**
 * Run the time-stepping loop from the current time until the time
 * specified by the IceModel::grid::time object.
 *
 * This is the method used by PISM in the "standalone" mode.
 */
IceModelTerminationReason IceModel::run() {
  const Profiling &profiling = m_ctx->profiling();

  bool do_mass_conserve = m_config->get_flag("geometry.update.enabled");
  bool do_energy = m_config->get_flag("energy.enabled");
  bool do_skip = m_config->get_flag("time_stepping.skip.enabled");

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
    for (const auto &d : m_ts_diagnostics) {
      d.second->update(time, time);
    }
  }

  m_log->message(2, "running forward ...\n");

  m_stdout_flags.erase(); // clear it out
  print_summary_line(true, do_energy, 0.0, 0.0, 0.0, 0.0, 0.0);
  print_summary(do_energy);  // report starting state

  t_TempAge = m_time->current();
  dt_TempAge = 0.0;

  IceModelTerminationReason termination_reason = PISM_DONE;
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

    double time_resolution = m_config->get_number("time_stepping.resolution");
    const bool show_step = tempAgeStep or fabs(m_time->current() - m_time->end()) < time_resolution;
    print_summary(show_step);

    // update viewers before writing extras because writing extras resets diagnostics
    update_viewers();

    // writing these fields here ensures that we do it after the last time-step
    profiling.begin("io");
    write_snapshot();
    write_extras();
    bool stop_after_chekpoint = write_checkpoint();
    profiling.end("io");

    if (stop_after_chekpoint) {
      termination_reason = PISM_CHEKPOINT;
      break;
    }

    if (process_signals() != 0) {
      termination_reason = PISM_SIGNAL;
      break;
    }
  } // end of the time-stepping loop
  profiling.stage_end("time-stepping loop");

  return termination_reason;
}

//! Manage the initialization of the IceModel object.
/*!
Please see the documenting comments of the functions called below to find
explanations of their intended uses.
 */
void IceModel::init() {
  // Get the start time in seconds and ensure that it is consistent
  // across all processors.
  m_start_time = get_time(m_grid->com);

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

const energy::BedThermalUnit* IceModel::bedrock_thermal_model() const {
  return m_btu.get();
}

const energy::EnergyModel* IceModel::energy_balance_model() const {
  return m_energy_model.get();
}

const YieldStress* IceModel::basal_yield_stress_model() const {
  return m_basal_yield_stress_model.get();
}

const bed::BedDef* IceModel::bed_deformation_model() const {
  return m_beddef.get();
}

/*!
 * Return thickness change due to calving (over the last time step).
 */
const array::Scalar& IceModel::calving() const {
  return m_thickness_change.calving;
}

/*!
 * Return thickness change due to frontal melt (over the last time step).
 */
const array::Scalar& IceModel::frontal_melt() const {
  return m_thickness_change.frontal_melt;
}

/*!
 * Return thickness change due to forced retreat (over the last time step).
 */
const array::Scalar& IceModel::forced_retreat() const {
  return m_thickness_change.forced_retreat;
}

void warn_about_missing(const Logger &log,
                        const std::set<std::string> &vars,
                        const std::string &type,
                        const std::set<std::string> &available,
                        bool stop) {
  std::vector<std::string> missing;
  for (const auto &v : vars) {
    if (available.find(v) == available.end()) {
      missing.push_back(v);
    }
  }

  if (not missing.empty()) {
    size_t N = missing.size();
    const char *ending = N > 1 ? "s" : "";
    const char *verb   = N > 1 ? "are" : "is";
    if (stop) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "%s variable%s %s %s not available!\n"
                                    "Available variables:\n- %s",
                                    type.c_str(),
                                    ending,
                                    join(missing, ",").c_str(),
                                    verb,
                                    set_join(available, ",\n- ").c_str());
    }

    log.message(2,
                "\nWARNING: %s variable%s %s %s not available!\n\n",
                type.c_str(),
                ending,
                join(missing, ",").c_str(),
                verb);
  }
}

/*!
 * De-allocate diagnostics that were not requested.
 *
 * Checks viewers, -spatial_vars, -checkpoint, -save_vars, and regular output.
 *
 * FIXME: I need to make sure that these reporting mechanisms are active. It is possible that
 * variables are on a list, but that list is not actually used.
 */
void IceModel::prune_diagnostics() {

  // get the list of available diagnostics
  std::set<std::string> available;
  for (const auto &d : m_diagnostics) {
    available.insert(d.first);
  }

  auto m_extra_stop = m_config->get_flag("output.extra.stop_missing");
  warn_about_missing(*m_log, m_output_vars,     "output",     available, false);
  warn_about_missing(*m_log, m_snapshot_vars,   "snapshot",   available, false);
  warn_about_missing(*m_log, m_checkpoint_vars, "checkpoint", available, false);
  warn_about_missing(*m_log, m_extra_vars,      "diagnostic", available, m_extra_stop);

  // get the list of requested diagnostics
  auto requested = set_split(m_config->get_string("output.runtime.viewer.variables"), ',');
  requested = combine(requested, m_output_vars);
  requested = combine(requested, m_snapshot_vars);
  requested = combine(requested, m_extra_vars);
  requested = combine(requested, m_checkpoint_vars);

  // de-allocate diagnostics that were not requested
  for (const auto &v : available) {
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
    for (const auto &v : m_ts_vars) {
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
  for (const auto &d : m_diagnostics) {
    d.second->update(dt);
  }

  const double time = m_time->current();
  for (const auto &d : m_ts_diagnostics) {
    d.second->update(time - dt, time);
  }
}

/*!
 * Reset accumulators in diagnostics that compute time-averaged quantities.
 */
void IceModel::reset_diagnostics() {
  for (auto &d : m_diagnostics) {
    d.second->reset();
  }
}

IceModel::ThicknessChanges::ThicknessChanges(const std::shared_ptr<const Grid> &grid)
  : calving(grid, "thickness_change_due_to_calving"),
    frontal_melt(grid, "thickness_change_due_to_frontal_melt"),
    forced_retreat(grid, "thickness_change_due_to_forced_retreat") {
  // empty
}

void IceModel::set_python_ocean_model(std::shared_ptr<ocean::PyOceanModel> model) {
  m_ocean = std::make_shared<ocean::PyOceanModelAdapter>(m_grid, model);
  m_submodels["ocean model"] = m_ocean.get();
}

} // end of namespace pism
