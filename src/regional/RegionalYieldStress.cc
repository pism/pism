/* Copyright (C) 2015, 2016, 2017, 2018, 2019 PISM Authors
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
#include "pism/util/pism_utilities.hh"

namespace pism {

RegionalYieldStress::RegionalYieldStress(std::shared_ptr<YieldStress> input)
  : YieldStress(input->grid()), m_input(input) {

  m_high_tauc = m_config->get_number("regional.no_model_yield_stress", "Pa");

  m_name = "regional " + m_input->name();
}

RegionalYieldStress::~RegionalYieldStress() {
  // empty
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

  // Read in tauc from the input file (this field would have been written by this class
  // and so it should contain the modification in "no model" areas).
  m_basal_yield_stress.read(input_file, record);
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

} // end of namespace pism
