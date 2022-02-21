/* Copyright (C) 2015, 2016, 2017, 2018, 2019, 2022 PISM Authors
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

#include "RegionalYieldStress.hh"
#include "pism/util/pism_utilities.hh" // pism::combine()
#include "pism/util/MaxTimestep.hh"
#include "pism/util/IceModelVec2Int.hh"

namespace pism {

RegionalYieldStress::RegionalYieldStress(std::shared_ptr<YieldStress> input)
  : YieldStress(input->grid()), m_input(input) {

  m_high_tauc = m_config->get_number("regional.no_model_yield_stress", "Pa");

  m_name = "regional " + m_input->name();
}

/*!
 * Set `basal_yield_stress` to `tauc` in areas indicated using `mask`.
 */
static void set_no_model_yield_stress(double tauc,
                                      const IceModelVec2Int &mask,
                                      IceModelVec2S &basal_yield_stress) {
  auto grid = mask.grid();

  IceModelVec::AccessList list{&mask, &basal_yield_stress};

  for (Points p(*grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    if (mask(i, j) > 0.5) {
      basal_yield_stress(i, j) = tauc;
    }
  }
}

void RegionalYieldStress::restart_impl(const File &input_file, int record) {
  m_input->restart(input_file, record);

  // m_basal_yield_stress is a part of the model state for all yield stress models, so the
  // call above should read it in.
  m_basal_yield_stress.copy_from(m_input->basal_material_yield_stress());

  IceModelVec2Int no_model_mask(m_grid, "no_model_mask", WITHOUT_GHOSTS);
  no_model_mask.set_attrs("model_state",
                          "mask: zeros (modeling domain) and ones"
                          " (no-model buffer near grid edges)",
                          "", "", "", 0); // no units and no standard name
  // We are re-starting a simulation, so the input file has to contain no_model_mask.
  no_model_mask.read(input_file, record);
  // However, the used can set "-regrid_vars no_model_mask,...", so we have to try this,
  // too.
  regrid(name(), no_model_mask);

  set_no_model_yield_stress(m_high_tauc, no_model_mask, m_basal_yield_stress);
}

void RegionalYieldStress::bootstrap_impl(const File &input_file, const YieldStressInputs &inputs) {
  m_input->bootstrap(input_file, inputs);

  m_basal_yield_stress.copy_from(m_input->basal_material_yield_stress());

  set_no_model_yield_stress(m_high_tauc, *inputs.no_model_mask, m_basal_yield_stress);
}

void RegionalYieldStress::init_impl(const YieldStressInputs &inputs) {
  m_input->init(inputs);

  m_basal_yield_stress.copy_from(m_input->basal_material_yield_stress());

  set_no_model_yield_stress(m_high_tauc, *inputs.no_model_mask, m_basal_yield_stress);
}

void RegionalYieldStress::update_impl(const YieldStressInputs &inputs,
                                      double t, double dt) {
  m_input->update(inputs, t, dt);

  m_basal_yield_stress.copy_from(m_input->basal_material_yield_stress());

  set_no_model_yield_stress(m_high_tauc, *inputs.no_model_mask, m_basal_yield_stress);
}

void RegionalYieldStress::define_model_state_impl(const File &output) const {
  m_input->define_model_state(output);

  // define tauc (this is likely to be a no-op because m_input should have defined it by
  // now)
  m_basal_yield_stress.define(output);
}

void RegionalYieldStress::write_model_state_impl(const File &output) const {
  m_input->write_model_state(output);
  // Write basal yield stress that includes the modification containing high yield stress
  // in "no model" areas, overwriting the field written by m_input.
  m_basal_yield_stress.write(output);
}

DiagnosticList RegionalYieldStress::diagnostics_impl() const {
  // Override the tauc diagnostic with the one that includes the regional modification
  return combine({{"tauc", Diagnostic::wrap(m_basal_yield_stress)}},
                 m_input->diagnostics());
}

MaxTimestep RegionalYieldStress::max_timestep_impl(double t) const {
  auto dt = m_input->max_timestep(t);

  if (dt.finite()) {
    return MaxTimestep(dt.value(), name());
  }

  return MaxTimestep(name());
}

TSDiagnosticList RegionalYieldStress::ts_diagnostics_impl() const {
  return m_input->ts_diagnostics();
}

} // end of namespace pism
