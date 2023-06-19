/* Copyright (C) 2015, 2016, 2017, 2018, 2019, 2021, 2022, 2023 PISM Authors
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

#include "YieldStress.hh"

#include "pism/util/ConfigInterface.hh"
#include "pism/util/Logger.hh"

namespace pism {

YieldStressInputs::YieldStressInputs() {
  geometry                   = nullptr;
  no_model_mask              = nullptr;
  till_water_thickness       = nullptr;
  subglacial_water_thickness = nullptr;
}

YieldStress::YieldStress(std::shared_ptr<const Grid> g)
  : Component(g),
  m_basal_yield_stress(m_grid, "tauc") {

  // PROPOSED standard_name = land_ice_basal_material_yield_stress
  m_basal_yield_stress.set_attrs("model_state",
                                 "yield stress for basal till (plastic or pseudo-plastic model)",
                                 "Pa", "Pa", "", 0);
}

/*!
 * Restart a yield stress model from an input file.
 */
void YieldStress::restart(const File &input_file, int record) {
  m_log->message(2, "* Initializing the %s...\n", name().c_str());

  this->restart_impl(input_file, record);
}

/*!
 * Bootstrap a yield stress model using incomplete inputs.
 */
void YieldStress::bootstrap(const File &input_file, const YieldStressInputs &inputs) {
  m_log->message(2, "Initializing the %s...\n", name().c_str());

  this->bootstrap_impl(input_file, inputs);
}

/*!
 * Initialize a yield stress model using inputs from other models and configuration
 * parameters.
 */
void YieldStress::init(const YieldStressInputs &inputs) {
  m_log->message(2, "Initializing the %s...\n", name().c_str());

  this->init_impl(inputs);
}

/*!
 * Update a yield stress model.
 */
void YieldStress::update(const YieldStressInputs &inputs, double t, double dt) {
  this->update_impl(inputs, t, dt);
}

const array::Scalar& YieldStress::basal_material_yield_stress() {
  return m_basal_yield_stress;
}

/*!
 * Define model state variables.
 *
 * All yield stress models have to write basal yield stress to output files and read it
 * from and input file during initialization because yield stress may be used by PISM's
 * stress balance model. The stress balance code has to be executed early during an update
 * of the model because its output (ice velocity) is used to compute the maximum allowed
 * time step.
 *
 * Now that PISM's yield stress models are time-dependent YieldStress::update() will be
 * called *after* the maximum time step is found. This means that during the first time
 * step basal_material_yield_stress() gets called before update().
 */
void YieldStress::define_model_state_impl(const File &output) const {
  m_basal_yield_stress.define(output, io::PISM_DOUBLE);
}

void YieldStress::write_model_state_impl(const File &output) const {
  m_basal_yield_stress.write(output);
}

DiagnosticList YieldStress::diagnostics_impl() const {
  return {{"tauc", Diagnostic::wrap(m_basal_yield_stress)}};
}

std::string YieldStress::name() const {
  return m_name;
}

} // end of namespace pism
