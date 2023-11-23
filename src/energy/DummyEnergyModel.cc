/* Copyright (C) 2016, 2017, 2023 PISM Authors
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

#include "pism/energy/EnthalpyModel.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace energy {

DummyEnergyModel::DummyEnergyModel(
    std::shared_ptr<const Grid> grid,
    std::shared_ptr<const stressbalance::StressBalance> stress_balance)
    : EnthalpyModel(grid, stress_balance) {
  // empty
}

void DummyEnergyModel::restart_impl(const File &input_file, int record) {
  EnthalpyModel::restart_impl(input_file, record);

  m_log->message(2,
                 "NOTE: this \"energy balance\" model holds enthalpy and basal melt rate constant in time.\n");
}

void DummyEnergyModel::bootstrap_impl(const File &input_file,
                                   const array::Scalar &ice_thickness,
                                   const array::Scalar &surface_temperature,
                                   const array::Scalar &climatic_mass_balance,
                                   const array::Scalar &basal_heat_flux) {
  EnthalpyModel::bootstrap_impl(input_file,
                                ice_thickness, surface_temperature,
                                climatic_mass_balance, basal_heat_flux);
  m_log->message(2,
                 "NOTE: this \"energy balance\" model holds enthalpy and basal melt rate constant in time.\n");
}

void DummyEnergyModel::update_impl(double t, double dt, const Inputs &inputs) {
  (void) t;
  (void) dt;
  (void) inputs;
}

MaxTimestep DummyEnergyModel::max_timestep_impl(double t) const {
  // silence a compiler warning
  (void) t;

  // no time step restriction
  return MaxTimestep("dummy energy model");
}


} // end of namespace energy
} // end of namespace pism
