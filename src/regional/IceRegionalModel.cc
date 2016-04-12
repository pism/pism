/* Copyright (C) 2015, 2016 PISM Authors
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
#include "base/util/PISMVars.hh"
#include "base/util/pism_options.hh"
#include "base/enthalpyConverter.hh"
#include "base/stressbalance/ShallowStressBalance.hh"
#include "SSAFD_Regional.hh"
#include "base/stressbalance/SSB_Modifier.hh"
#include "SIAFD_Regional.hh"
#include "base/stressbalance/PISMStressBalance.hh"
#include "base/basalstrength/PISMConstantYieldStress.hh"
#include "RegionalDefaultYieldStress.hh"
#include "base/util/io/PIO.hh"

namespace pism {

//! \brief Set no_model_mask variable to have value 1 in strip of width 'strip'
//! m around edge of computational domain, and value 0 otherwise.
static void set_no_model_strip(const IceGrid &grid, double width, IceModelVec2Int &result) {

  IceModelVec::AccessList list(result);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (in_null_strip(grid, i, j, width) == true) {
      result(i, j) = 1;
    } else {
      result(i, j) = 0;
    }
  }

  result.metadata().set_string("pism_intent", "model_state");

  result.update_ghosts();
}


IceRegionalModel::IceRegionalModel(IceGrid::Ptr g, Context::Ptr c)
  : IceModel(g, c) {
  // empty
}


void IceRegionalModel::createVecs() {

  IceModel::createVecs();

  m_log->message(2, 
             "  creating IceRegionalModel vecs ...\n");

  // stencil width of 2 needed for surfaceGradientSIA() action
  m_no_model_mask.create(m_grid, "no_model_mask", WITH_GHOSTS, 2);
  m_no_model_mask.set_attrs("model_state", // ensures that it gets written at the end of the run
                          "mask: zeros (modeling domain) and ones (no-model buffer near grid edges)",
                          "", ""); // no units and no standard name
  double NMMASK_NORMAL   = 0.0,
         NMMASK_ZERO_OUT = 1.0;
  std::vector<double> mask_values(2);
  mask_values[0] = NMMASK_NORMAL;
  mask_values[1] = NMMASK_ZERO_OUT;
  m_no_model_mask.metadata().set_doubles("flag_values", mask_values);
  m_no_model_mask.metadata().set_string("flag_meanings", "normal special_treatment");
  m_no_model_mask.set_time_independent(true);
  m_no_model_mask.set(NMMASK_NORMAL);
  m_grid->variables().add(m_no_model_mask);

  // stencil width of 2 needed for differentiation because GHOSTS=1
  m_usurf_stored.create(m_grid, "usurfstore", WITH_GHOSTS, 2);
  m_usurf_stored.set_attrs("model_state", // ensures that it gets written at the end of the run
                       "saved surface elevation for use to keep surface gradient constant in no_model strip",
                       "m",
                       ""); //  no standard name
  m_grid->variables().add(m_usurf_stored);

  // stencil width of 1 needed for differentiation
  m_thk_stored.create(m_grid, "thkstore", WITH_GHOSTS, 1);
  m_thk_stored.set_attrs("model_state", // ensures that it gets written at the end of the run
                     "saved ice thickness for use to keep driving stress constant in no_model strip",
                     "m",
                     ""); //  no standard name
  m_grid->variables().add(m_thk_stored);

  // Note that the name of this variable (bmr_stored) does not matter: it is
  // *never* read or written. We make a copy of bmelt instead.
  m_bmr_stored.create(m_grid, "bmr_stored", WITH_GHOSTS, 2);
  m_bmr_stored.set_attrs("internal",
                       "time-independent basal melt rate in the no-model-strip",
                       "m s-1", "");

  if (m_config->get_boolean("ssa_dirichlet_bc")) {
    // remove the bc_mask variable from the dictionary
    m_grid->variables().remove("bc_mask");

    m_grid->variables().add(m_no_model_mask, "bc_mask");
  }
}

void IceRegionalModel::model_state_setup() {

  IceModel::model_state_setup();

  // Now save the basal melt rate at the beginning of the run.
  m_bmr_stored.copy_from(m_basal_melt_rate);

  bool zgwnm = options::Bool("-zero_grad_where_no_model",
                             "set zero surface gradient in no model strip");
  if (zgwnm) {
    m_thk_stored.set(0.0);
    m_usurf_stored.set(0.0);
  }

  options::Real strip_km("-no_model_strip", 
                        "width in km of strip near boundary in which modeling is turned off",
                        0.0);

  if (strip_km.is_set()) {
    m_log->message(2, 
                   "* Option -no_model_strip read... setting boundary strip width to %.2f km\n",
                   strip_km.value());
    set_no_model_strip(*m_grid,
                       units::convert(m_sys, strip_km, "km", "m"),
                       m_no_model_mask);
  } else {
    double strip = m_config->get_double("regional_no_model_strip", "m");
    m_log->message(2,
                   "* Option -no_model_strip is not set... setting boundary strip width to %.2f km\n",
                   units::convert(m_sys, strip, "m", "km"));
    set_no_model_strip(*m_grid, strip, m_no_model_mask);
  }
}

void IceRegionalModel::allocate_stressbalance() {

  using namespace pism::stressbalance;

  if (m_stress_balance != NULL) {
    return;
  }

  EnthalpyConverter::Ptr EC = m_ctx->enthalpy_converter();

  std::string model = m_config->get_string("stress_balance_model");

  ShallowStressBalance *sliding = NULL;
  if (model == "none" || model == "sia") {
    sliding = new ZeroSliding(m_grid, EC);
  } else if (model == "prescribed_sliding" || model == "prescribed_sliding+sia") {
    sliding = new PrescribedSliding(m_grid, EC);
  } else if (model == "ssa" || model == "ssa+sia") {
    sliding = new SSAFD_Regional(m_grid, EC);
  } else {
    throw RuntimeError::formatted("invalid stress balance model: %s", model.c_str());
  }

  SSB_Modifier *modifier = NULL;
  if (model == "none" || model == "ssa" || model == "prescribed_sliding") {
    modifier = new ConstantInColumn(m_grid, EC);
  } else if (model == "prescribed_sliding+sia" ||
             model == "ssa+sia" ||
             model == "sia") {
    modifier = new SIAFD_Regional(m_grid, EC);
  } else {
    throw RuntimeError::formatted("invalid stress balance model: %s", model.c_str());
  }

  // ~StressBalance() will de-allocate sliding and modifier.
  m_stress_balance = new StressBalance(m_grid, sliding, modifier);
}


void IceRegionalModel::allocate_basal_yield_stress() {

  if (basal_yield_stress_model != NULL) {
    return;
  }

  std::string model = m_config->get_string("stress_balance_model");

  // only these two use the yield stress (so far):
  if (model == "ssa" || model == "ssa+sia") {
    std::string yield_stress_model = m_config->get_string("yield_stress_model");

    if (yield_stress_model == "constant") {
      basal_yield_stress_model = new ConstantYieldStress(m_grid);
    } else if (yield_stress_model == "mohr_coulomb") {
      basal_yield_stress_model = new RegionalDefaultYieldStress(m_grid, subglacial_hydrology);
    } else {
      throw RuntimeError::formatted("yield stress model '%s' is not supported.",
                                    yield_stress_model.c_str());
    }
  }
}


void IceRegionalModel::bootstrap_2d(const std::string &filename) {

  IceModel::bootstrap_2d(filename);

  // read usurfstore from usurf, then restore its name
  m_usurf_stored.metadata().set_name("usurf");
  m_usurf_stored.regrid(filename, OPTIONAL, 0.0);
  m_usurf_stored.metadata().set_name("usurfstore");

  // read thkstore from thk, then restore its name
  m_thk_stored.metadata().set_name("thk");
  m_thk_stored.regrid(filename, OPTIONAL, 0.0);
  m_thk_stored.metadata().set_name("thkstore");
}


void IceRegionalModel::initFromFile(const std::string &filename) {
  PIO nc(m_grid->com, "guess_mode");

  bool no_model_strip_set = options::Bool("-no_model_strip", "No-model strip, in km");

  if (no_model_strip_set) {
    m_no_model_mask.metadata().set_string("pism_intent", "internal");
  }

  m_log->message(2, 
             "* Initializing IceRegionalModel from NetCDF file '%s'...\n",
             filename.c_str());

  // Allow re-starting from a file that does not contain u_ssa_bc and v_ssa_bc.
  // The user is probably using -regrid_file to bring in SSA B.C. data.
  if (m_config->get_boolean("ssa_dirichlet_bc")) {
    bool u_ssa_exists, v_ssa_exists;

    nc.open(filename, PISM_READONLY);
    u_ssa_exists = nc.inq_var("u_ssa_bc");
    v_ssa_exists = nc.inq_var("v_ssa_bc");
    nc.close();

    if (! (u_ssa_exists && v_ssa_exists)) {
      m_ssa_dirichlet_bc_values.metadata().set_string("pism_intent", "internal");
      m_log->message(2, 
                 "PISM WARNING: u_ssa_bc and/or v_ssa_bc not found in %s. Setting them to zero.\n"
                 "              This may be overridden by the -regrid_file option.\n",
                 filename.c_str());

      m_ssa_dirichlet_bc_values.set(0.0);
    }
  }

  bool zgwnm = options::Bool("-zero_grad_where_no_model",
                             "zero surface gradient in no model strip");
  if (zgwnm) {
    m_thk_stored.metadata().set_string("pism_intent", "internal");
    m_usurf_stored.metadata().set_string("pism_intent", "internal");
  }

  IceModel::initFromFile(filename);

  if (m_config->get_boolean("ssa_dirichlet_bc")) {
      m_ssa_dirichlet_bc_values.metadata().set_string("pism_intent", "model_state");
  }

  if (zgwnm) {
    m_thk_stored.metadata().set_string("pism_intent", "model_state");
    m_usurf_stored.metadata().set_string("pism_intent", "model_state");
  }
}


void IceRegionalModel::set_vars_from_options() {

  // base class reads the -bootstrap option and does the bootstrapping:
  IceModel::set_vars_from_options();

  if (m_config->get_boolean("do_cold_ice_methods")) {
    throw RuntimeError("pismo does not support the 'cold' mode.");
  }
}

void IceRegionalModel::massContExplicitStep() {

  // This ensures that no_model_mask is available in
  // IceRegionalModel::cell_interface_fluxes() below.
  IceModelVec::AccessList list(m_no_model_mask);

  IceModel::massContExplicitStep();
}

void IceRegionalModel::cell_interface_fluxes(bool dirichlet_bc,
                                             int i, int j,
                                             StarStencil<Vector2> input_velocity,
                                             StarStencil<double> input_flux,
                                             StarStencil<double> &output_velocity,
                                             StarStencil<double> &output_flux) {

  IceModel::cell_interface_fluxes(dirichlet_bc, i, j,
                                  input_velocity,
                                  input_flux,
                                  output_velocity,
                                  output_flux);

  StarStencil<int> nmm = m_no_model_mask.int_star(i,j);
  Direction dirs[4] = {North, East, South, West};

  for (int n = 0; n < 4; ++n) {
    Direction direction = dirs[n];

      if ((nmm.ij == 1) || (nmm.ij == 0 && nmm[direction] == 1)) {
      output_velocity[direction] = 0.0;
      output_flux[direction] = 0.0;
    }
  }
  //
}

void IceRegionalModel::enthalpyAndDrainageStep(unsigned int *vertSacrCount,
                                               double *liquifiedVol,
                                               unsigned int *bulgeCount) {

  IceModel::enthalpyAndDrainageStep(vertSacrCount, liquifiedVol, bulgeCount);

  // note that the call above sets vWork3d; ghosts are comminucated later (in
  // IceModel::energyStep()).
  IceModelVec::AccessList list;
  list.add(m_no_model_mask);
  list.add(vWork3d);
  list.add(m_ice_enthalpy);

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_no_model_mask(i, j) < 0.5) {
      continue;
    }

    double *new_enthalpy = vWork3d.get_column(i, j);
    double *old_enthalpy = m_ice_enthalpy.get_column(i, j);

    for (unsigned int k = 0; k < m_grid->Mz(); ++k) {
      new_enthalpy[k] = old_enthalpy[k];
    }
  }

  // set basal_melt_rate; ghosts are comminucated later (in IceModel::energyStep()).
  list.add(m_basal_melt_rate);
  list.add(m_bmr_stored);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (m_no_model_mask(i, j) < 0.5) {
      continue;
    }

    m_basal_melt_rate(i, j) = m_bmr_stored(i, j);
  }
}

} // end of namespace pism
