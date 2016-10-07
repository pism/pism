/* Copyright (C) 2016 PISM Authors
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

#include "EnergyModel.hh"

namespace pism {
namespace energy {

TemperatureModel::TemperatureModel(IceGrid::ConstPtr grid,
                                   stressbalance::StressBalance *stress_balance)
  : EnergyModel(grid, stress_balance) {

  m_ice_temperature.create(m_grid, "temp", WITH_GHOSTS);
  m_ice_temperature.set_attrs("model_state",
                              "ice temperature", "K", "land_ice_temperature");
  m_ice_temperature.metadata().set_double("valid_min", 0.0);
}

void TemperatureModel::init_impl(const InputOptions &opts) {

}

void TemperatureModel::update_impl(double t, double dt, const EnergyModelInputs &inputs) {

}

void TemperatureModel::define_model_state_impl(const PIO &output) const {
  m_ice_temperature.define(output);
  m_basal_melt_rate.define(output);
}

void TemperatureModel::write_model_state_impl(const PIO &output) const {
  m_ice_temperature.write(output);
  m_basal_melt_rate.write(output);
}

} // end of namespace energy
} // end of namespace pism
