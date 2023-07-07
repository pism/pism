/* Copyright (C) 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "pism/regional/IceRegionalModel.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/regional/RegionalYieldStress.hh"
#include "pism/util/io/File.hh"
#include "pism/coupler/SurfaceModel.hh"
#include "pism/regional/EnthalpyModel_Regional.hh"
#include "pism/energy/CHSystem.hh"
#include "pism/energy/BedThermalUnit.hh"
#include "pism/energy/utilities.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/hydrology/Hydrology.hh"

namespace pism {

IceRegionalModel::IceRegionalModel(std::shared_ptr<Grid> g, std::shared_ptr<Context> c)
    : IceModel(g, c),
      m_no_model_mask(m_grid, "no_model_mask"),
      m_usurf_stored(m_grid, "usurfstore"),
      m_thk_stored(m_grid, "thkstore") {
  m_no_model_mask.set_interpolation_type(NEAREST);

  if (m_config->get_flag("energy.ch_warming.enabled")) {
    m_ch_warming_flux.reset(
        new array::Array3D(m_grid, "ch_warming_flux", array::WITHOUT_GHOSTS, m_grid->z()));
  }
}


void IceRegionalModel::allocate_storage() {

  IceModel::allocate_storage();

  m_log->message(2, "  creating IceRegionalModel vecs ...\n");

  // stencil width of 2 needed by SIAFD_Regional::compute_surface_gradient()
  m_no_model_mask.metadata(0)
      .long_name("mask: zeros (modeling domain) and ones (no-model buffer near grid edges)")
      .set_time_independent(true)
      .set_output_type(io::PISM_INT); // no units and no standard name
  m_no_model_mask.metadata()["flag_values"]   = { 0, 1 };
  m_no_model_mask.metadata()["flag_meanings"] = "normal special_treatment";

  m_no_model_mask.set(0);

  // stencil width of 2 needed for differentiation because GHOSTS=1
  m_usurf_stored.metadata(0)
      .long_name(
          "saved surface elevation for use to keep surface gradient constant in no_model strip")
      .units("m"); //  no standard name

  // stencil width of 1 needed for differentiation
  m_thk_stored.metadata(0)
      .long_name("saved ice thickness for use to keep driving stress constant in no_model strip")
      .units("m"); //  no standard name

  m_model_state.insert(&m_thk_stored);
  m_model_state.insert(&m_usurf_stored);
  m_model_state.insert(&m_no_model_mask);
}

void IceRegionalModel::model_state_setup() {

  // initialize the model state (including special fields)
  IceModel::model_state_setup();

  InputOptions input = process_input_options(m_ctx->com(), m_config);

  // Initialize stored ice thickness and surface elevation. This goes here and not in
  // bootstrap_2d because bed topography is not initialized at the time bootstrap_2d is
  // called.
  if (input.type == INIT_BOOTSTRAP) {
    if (m_config->get_flag("regional.zero_gradient")) {
      m_usurf_stored.set(0.0);
      m_thk_stored.set(0.0);
    } else {
      m_usurf_stored.copy_from(m_geometry.ice_surface_elevation);
      m_thk_stored.copy_from(m_geometry.ice_thickness);
    }
  }

  m_geometry_evolution->set_no_model_mask(m_no_model_mask);

  if (m_ch_system) {
    const bool use_input_file = input.type == INIT_BOOTSTRAP or input.type == INIT_RESTART;

    std::unique_ptr<File> input_file;

    if (use_input_file) {
      input_file.reset(new File(m_grid->com, input.filename, io::PISM_GUESS, io::PISM_READONLY));
    }

    switch (input.type) {
    case INIT_RESTART:
      {
        m_ch_system->restart(*input_file, input.record);
        break;
      }
    case INIT_BOOTSTRAP:
      {

        m_ch_system->bootstrap(*input_file,
                               m_geometry.ice_thickness,
                               m_surface->temperature(),
                               m_surface->mass_flux(),
                               m_btu->flux_through_top_surface());
        break;
      }
    case INIT_OTHER:
    default:
      {
        m_basal_melt_rate.set(m_config->get_number("bootstrapping.defaults.bmelt"));

        m_ch_system->initialize(m_basal_melt_rate,
                                m_geometry.ice_thickness,
                                m_surface->temperature(),
                                m_surface->mass_flux(),
                                m_btu->flux_through_top_surface());

      }
    }
  }
}

void IceRegionalModel::allocate_geometry_evolution() {
  if (m_geometry_evolution) {
    return;
  }

  m_log->message(2, "# Allocating the geometry evolution model (regional version)...\n");

  m_geometry_evolution.reset(new RegionalGeometryEvolution(m_grid));

  m_submodels["geometry_evolution"] = m_geometry_evolution.get();
}

void IceRegionalModel::allocate_energy_model() {

  if (m_energy_model != NULL) {
    return;
  }

  m_log->message(2, "# Allocating an energy balance model...\n");

  if (m_config->get_flag("energy.enabled")) {
    if (m_config->get_flag("energy.temperature_based")) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "pismr -regional does not support the '-energy cold' mode.");
    }

    m_energy_model =
        std::make_shared<energy::EnthalpyModel_Regional>(m_grid, m_stress_balance.get());
  } else {
    m_energy_model = std::make_shared<energy::DummyEnergyModel>(m_grid, m_stress_balance.get());
  }

  m_submodels["energy balance model"] = m_energy_model.get();

  if (m_config->get_flag("energy.ch_warming.enabled") and
      not m_ch_system) {

    m_log->message(2, "# Allocating the cryo-hydrologic warming model...\n");

    m_ch_system.reset(new energy::CHSystem(m_grid, m_stress_balance.get()));
    m_submodels["cryo-hydrologic warming"] = m_ch_system.get();
  }
}

void IceRegionalModel::allocate_stressbalance() {

  if (m_stress_balance) {
    return;
  }

  bool regional = true;
  m_stress_balance = stressbalance::create(m_config->get_string("stress_balance.model"),
                                           m_grid, regional);

  m_submodels["stress balance"] = m_stress_balance.get();
}


void IceRegionalModel::allocate_basal_yield_stress() {

  if (m_basal_yield_stress_model) {
    return;
  }

  IceModel::allocate_basal_yield_stress();

  if (m_basal_yield_stress_model) {
    // IceModel allocated a basal yield stress model. This means that we are using a
    // stress balance model that uses it and we need to add regional modifications.
    m_basal_yield_stress_model.reset(new RegionalYieldStress(m_basal_yield_stress_model));

    m_submodels["basal yield stress"] = m_basal_yield_stress_model.get();
  }
}

//! Bootstrap a "regional" model.
/*!
 * Need to initialize all the variables managed by IceModel, as well as
 * - no_model_mask
 * - usurf_stored
 * - thk_stored
 */
void IceRegionalModel::bootstrap_2d(const File &input_file) {

  IceModel::bootstrap_2d(input_file);

  // no_model_mask
  {
    if (input_file.find_variable(m_no_model_mask.metadata().get_name())) {
      m_no_model_mask.regrid(input_file, io::CRITICAL);
    } else {
      // set using the no_model_strip parameter
      double strip_width = m_config->get_number("regional.no_model_strip", "meters");
      set_no_model_strip(*m_grid, strip_width, m_no_model_mask);
    }

    // m_no_model_mask was added to m_model_state, so
    // IceModel::regrid() will take care of regridding it, if necessary.
  }

  if (m_config->get_flag("stress_balance.ssa.dirichlet_bc")) {
    array::AccessScope list
      {&m_no_model_mask, &m_velocity_bc_mask, &m_ice_thickness_bc_mask};

    for (auto p = m_grid->points(); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_no_model_mask(i, j) > 0.5) {
        m_velocity_bc_mask(i, j)      = 1;
        m_ice_thickness_bc_mask(i, j) = 1;
      }
    }
  }
}

stressbalance::Inputs IceRegionalModel::stress_balance_inputs() {
  stressbalance::Inputs result = IceModel::stress_balance_inputs();

  result.no_model_mask              = &m_no_model_mask;
  result.no_model_ice_thickness     = &m_thk_stored;
  result.no_model_surface_elevation = &m_usurf_stored;

  return result;
}

energy::Inputs IceRegionalModel::energy_model_inputs() {
  energy::Inputs result = IceModel::energy_model_inputs();

  result.no_model_mask = &m_no_model_mask;

  return result;
}

void IceRegionalModel::energy_step() {

  if (m_ch_system) {
    bedrock_thermal_model_step();

    energy::Inputs inputs = energy_model_inputs();
    const array::Array3D *strain_heating = inputs.volumetric_heating_rate;
    inputs.volumetric_heating_rate = m_ch_warming_flux.get();

    energy::cryo_hydrologic_warming_flux(m_config->get_number("constants.ice.thermal_conductivity"),
                                         m_config->get_number("energy.ch_warming.average_channel_spacing"),
                                         m_geometry.ice_thickness,
                                         m_energy_model->enthalpy(),
                                         m_ch_system->enthalpy(),
                                         *m_ch_warming_flux);

    // Convert to the loss of energy by the CH system:
    m_ch_warming_flux->scale(-1.0);

    m_ch_system->update(t_TempAge, dt_TempAge, inputs);

    // Add CH warming flux to the strain heating term:
    m_ch_warming_flux->scale(-1.0);
    m_ch_warming_flux->add(1.0, *strain_heating);

    m_energy_model->update(t_TempAge, dt_TempAge, inputs);

    m_stdout_flags = m_energy_model->stdout_flags() + m_stdout_flags;
  } else {
    IceModel::energy_step();
  }
}

YieldStressInputs IceRegionalModel::yield_stress_inputs() {
  YieldStressInputs result = IceModel::yield_stress_inputs();

  result.no_model_mask = &m_no_model_mask;

  return result;
}

const energy::CHSystem* IceRegionalModel::cryo_hydrologic_system() const {
  return m_ch_system.get();
}

/*! @brief Report temperature of the cryo-hydrologic system */
class CHTemperature : public Diag<IceRegionalModel>
{
public:
  CHTemperature(const IceRegionalModel *m)
    : Diag<IceRegionalModel>(m) {

    m_vars = {SpatialVariableMetadata(m_sys, "ch_temp", m_grid->z())};

    set_attrs("temperature of the cryo-hydrologic system", "",
              "Kelvin", "Kelvin", 0);
  }

protected:
  std::shared_ptr<array::Array> compute_impl() const {

    std::shared_ptr<array::Array3D> result(new array::Array3D(m_grid, "ch_temp", array::WITHOUT_GHOSTS, m_grid->z()));

    energy::compute_temperature(model->cryo_hydrologic_system()->enthalpy(),
                                model->geometry().ice_thickness,
                                *result);
    result->metadata(0) = m_vars[0];

    return result;
  }
};

/*! @brief Report liquid water fraction in the cryo-hydrologic system */
class CHLiquidWaterFraction : public Diag<IceRegionalModel>
{
public:
  CHLiquidWaterFraction(const IceRegionalModel *m)
    : Diag<IceRegionalModel>(m) {

    m_vars = {SpatialVariableMetadata(m_sys, "ch_liqfrac", m_grid->z())};

    set_attrs("liquid water fraction in the cryo-hydrologic system", "",
              "1", "1", 0);
  }

protected:
  std::shared_ptr<array::Array> compute_impl() const {

    std::shared_ptr<array::Array3D> result(new array::Array3D(m_grid, "ch_liqfrac", array::WITHOUT_GHOSTS, m_grid->z()));

    energy::compute_liquid_water_fraction(model->cryo_hydrologic_system()->enthalpy(),
                                          model->geometry().ice_thickness,
                                          *result);
    result->metadata(0) = m_vars[0];
    return result;
  }
};


/*! @brief Report rate of cryo-hydrologic warming */
class CHHeatFlux : public Diag<IceRegionalModel>
{
public:
  CHHeatFlux(const IceRegionalModel *m)
    : Diag<IceRegionalModel>(m) {

    m_vars = {SpatialVariableMetadata(m_sys, "ch_heat_flux", m_grid->z())};

    set_attrs("rate of cryo-hydrologic warming", "",
              "W m-3", "W m-3", 0);
  }

protected:
  std::shared_ptr<array::Array> compute_impl() const {

    std::shared_ptr<array::Array3D> result(new array::Array3D(m_grid, "ch_heat_flux", array::WITHOUT_GHOSTS, m_grid->z()));
    result->metadata(0) = m_vars[0];

    energy::cryo_hydrologic_warming_flux(m_config->get_number("constants.ice.thermal_conductivity"),
                                         m_config->get_number("energy.ch_warming.average_channel_spacing"),
                                         model->geometry().ice_thickness,
                                         model->energy_balance_model()->enthalpy(),
                                         model->cryo_hydrologic_system()->enthalpy(),
                                         *result);
    return result;
  }
};

void IceRegionalModel::init_diagnostics() {
  IceModel::init_diagnostics();

  if (m_ch_system) {
    m_diagnostics["ch_temp"]      = Diagnostic::Ptr(new CHTemperature(this));
    m_diagnostics["ch_liqfrac"]   = Diagnostic::Ptr(new CHLiquidWaterFraction(this));
    m_diagnostics["ch_heat_flux"] = Diagnostic::Ptr(new CHHeatFlux(this));
  }
}

void IceRegionalModel::hydrology_step() {
  hydrology::Inputs inputs;

  array::Scalar &sliding_speed = *m_work2d[0];
  compute_magnitude(m_stress_balance->advective_velocity(), sliding_speed);

  inputs.no_model_mask      = &m_no_model_mask;
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


} // end of namespace pism
