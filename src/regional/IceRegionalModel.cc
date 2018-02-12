/* Copyright (C) 2015, 2016, 2017, 2018 PISM Authors
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

#include "IceRegionalModel.hh"
#include "pism/util/Vars.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/stressbalance/ShallowStressBalance.hh"
#include "SSAFD_Regional.hh"
#include "pism/stressbalance/SSB_Modifier.hh"
#include "SIAFD_Regional.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/basalstrength/ConstantYieldStress.hh"
#include "RegionalDefaultYieldStress.hh"
#include "pism/util/io/PIO.hh"
#include "pism/coupler/OceanModel.hh"
#include "EnthalpyModel_Regional.hh"

namespace pism {

//! \brief Set no_model_mask variable to have value 1 in strip of width 'strip'
//! m around edge of computational domain, and value 0 otherwise.
static void set_no_model_strip(const IceGrid &grid, double width, IceModelVec2Int &result) {

  if (width <= 0.0) {
    return;
  }

  IceModelVec::AccessList list(result);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (in_null_strip(grid, i, j, width)) {
      result(i, j) = 1;
    } else {
      result(i, j) = 0;
    }
  }

  result.update_ghosts();
}


IceRegionalModel::IceRegionalModel(IceGrid::Ptr g, Context::Ptr c)
  : IceModel(g, c) {
  // empty
}


void IceRegionalModel::allocate_storage() {

  IceModel::allocate_storage();

  m_log->message(2, 
                 "  creating IceRegionalModel vecs ...\n");

  // stencil width of 2 needed by SIAFD_Regional::compute_surface_gradient()
  m_no_model_mask.create(m_grid, "no_model_mask", WITH_GHOSTS, 2);
  m_no_model_mask.set_attrs("model_state",
                            "mask: zeros (modeling domain) and ones"
                            " (no-model buffer near grid edges)",
                            "", ""); // no units and no standard name
  m_no_model_mask.metadata().set_doubles("flag_values", {0, 1});
  m_no_model_mask.metadata().set_string("flag_meanings", "normal special_treatment");
  m_no_model_mask.set_time_independent(true);
  m_no_model_mask.metadata().set_output_type(PISM_BYTE);
  m_no_model_mask.set(0);

  // stencil width of 2 needed for differentiation because GHOSTS=1
  m_usurf_stored.create(m_grid, "usurfstore", WITH_GHOSTS, 2);
  m_usurf_stored.set_attrs("model_state",
                           "saved surface elevation for use to keep surface gradient constant"
                           " in no_model strip",
                           "m",
                           ""); //  no standard name

  // stencil width of 1 needed for differentiation
  m_thk_stored.create(m_grid, "thkstore", WITH_GHOSTS, 1);
  m_thk_stored.set_attrs("model_state",
                         "saved ice thickness for use to keep driving stress constant"
                         " in no_model strip",
                         "m",
                         ""); //  no standard name

  m_model_state.insert(&m_thk_stored);
  m_model_state.insert(&m_usurf_stored);
  m_model_state.insert(&m_no_model_mask);
}

void IceRegionalModel::model_state_setup() {

  // initialize the model state (including special fields)
  IceModel::model_state_setup();

  InputOptions input = process_input_options(m_ctx->com());

  // Initialize stored ice thickness and surface elevation. This goes here and not in
  // bootstrap_2d because bed topography is not initialized at the time bootstrap_2d is
  // called.
  if (input.type == INIT_BOOTSTRAP) {

    if (m_config->get_boolean("regional.zero_gradient")) {
      m_usurf_stored.set(0.0);
      m_thk_stored.set(0.0);
    } else {
      GeometryCalculator gc(*m_config);
      // Use values set by IceModel::bootstrap_2d() to initialize usurf_stored and thk_stored.
      gc.compute_surface(m_ocean->sea_level_elevation(),
                         m_geometry.bed_elevation,
                         m_geometry.ice_thickness,
                         m_usurf_stored);

      m_thk_stored.copy_from(m_geometry.ice_thickness);
    }
  }

  // initialize the no_model_mask
  m_geometry_evolution->set_no_model_mask(m_no_model_mask);
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

  if (m_config->get_boolean("energy.enabled")) {
    if (m_config->get_boolean("energy.temperature_based")) {
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "pismr -regional does not support the '-energy cold' mode.");
    } else {
      m_energy_model = new energy::EnthalpyModel_Regional(m_grid, m_stress_balance.get());
    }
  } else {
    m_energy_model = new energy::DummyEnergyModel(m_grid, m_stress_balance.get());
  }

  m_submodels["energy balance model"] = m_energy_model;
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

  if (m_basal_yield_stress_model != NULL) {
    return;
  }

  std::string model = m_config->get_string("stress_balance.model");

  // only these two use the yield stress (so far):
  if (model == "ssa" || model == "ssa+sia") {
    std::string yield_stress_model = m_config->get_string("basal_yield_stress.model");

    if (yield_stress_model == "constant") {
      m_basal_yield_stress_model = new ConstantYieldStress(m_grid);
    } else if (yield_stress_model == "mohr_coulomb") {
      m_basal_yield_stress_model = new RegionalDefaultYieldStress(m_grid, m_subglacial_hydrology);
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "yield stress model '%s' is not supported.",
                                    yield_stress_model.c_str());
    }
    m_submodels["basal yield stress"] = m_basal_yield_stress_model;
  }
}

//! Bootstrap a "regional" model.
/*!
 * Need to initialize all the variables managed by IceModel, as well as
 * - no_model_mask
 * - usurf_stored
 * - thk_stored
 */
void IceRegionalModel::bootstrap_2d(const PIO &input_file) {

  IceModel::bootstrap_2d(input_file);

  // no_model_mask
  {
    // set using the no_model_strip parameter
    double strip_width = m_config->get_double("regional.no_model_strip", "meters");
    set_no_model_strip(*m_grid, strip_width, m_no_model_mask);

    // m_no_model_mask was added to m_model_state, so
    // IceModel::regrid() will take care of regridding it, if necessary.
  }

  if (m_config->get_boolean("stress_balance.ssa.dirichlet_bc")) {
    IceModelVec::AccessList list{&m_no_model_mask, &m_ssa_dirichlet_bc_mask};

    for (Points p(*m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      if (m_no_model_mask(i, j) > 0.5) {
        m_ssa_dirichlet_bc_mask(i, j) = 1;
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

YieldStressInputs IceRegionalModel::yield_stress_inputs() {
  YieldStressInputs result = IceModel::yield_stress_inputs();

  result.no_model_mask = &m_no_model_mask;

  return result;
}

} // end of namespace pism
