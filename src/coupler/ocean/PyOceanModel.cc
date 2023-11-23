/* Copyright (C) 2023 PISM Authors
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

#include "pism/coupler/ocean/PyOceanModel.hh"

#include "pism/util/MaxTimestep.hh"

#include <exception>
#include <memory>

namespace pism {
namespace ocean {

void PyOceanModel::allocate(std::shared_ptr<const Grid> grid) {
  shelf_base_mass_flux   = OceanModel::allocate_shelf_base_mass_flux(grid);
  shelf_base_temperature = OceanModel::allocate_shelf_base_temperature(grid);
  water_column_pressure  = OceanModel::allocate_water_column_pressure(grid);
}

MaxTimestep PyOceanModel::max_timestep(double t) const {
  (void) t;
  return {};
}

void PyOceanModel::init(const Geometry &geometry) {
  (void) geometry;
  throw RuntimeError(PISM_ERROR_LOCATION, "PyOceanModel.init(geometry) is not implemented");
}

void PyOceanModel::update(const Geometry &geometry, double t, double dt) {
  (void) geometry;
  (void) t;
  (void) dt;
  throw RuntimeError(PISM_ERROR_LOCATION, "PyOceanModel.update(geometry, t, dt) is not implemented");
}

void PyOceanModel::define_model_state(const File &output) const {
  (void) output;
  // empty
}

void PyOceanModel::write_model_state(const File &output) const {
  (void) output;
  // empty
}

PyOceanModelAdapter::PyOceanModelAdapter(std::shared_ptr<const Grid> grid,
                                         std::shared_ptr<PyOceanModel> implementation)
    : CompleteOceanModel(grid), m_impl(implementation) {

  m_impl->shelf_base_mass_flux = m_shelf_base_mass_flux;
  m_impl->shelf_base_temperature = m_shelf_base_temperature;
  m_impl->water_column_pressure = m_water_column_pressure;
}

MaxTimestep PyOceanModelAdapter::max_timestep_impl(double t) const {
  return m_impl->max_timestep(t);
}

void PyOceanModelAdapter::update_impl(const Geometry &geometry, double t, double dt) {
  m_impl->update(geometry, t, dt);
}

void PyOceanModelAdapter::init_impl(const Geometry &geometry) {
  m_impl->init(geometry);
}


void PyOceanModelAdapter::define_model_state_impl(const File &output) const {
  m_impl->define_model_state(output);
}

void PyOceanModelAdapter::write_model_state_impl(const File &output) const {
  m_impl->write_model_state(output);
}

} // namespace ocean
} // namespace pism
