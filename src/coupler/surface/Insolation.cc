/* Copyright (C) 2026 PISM Authors
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
#include <memory>

#include "pism/coupler/surface/DEBMSimplePointwise.hh"

#include "pism/coupler/surface/Insolation.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/util/Time.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/io/IO_Flags.hh"

namespace pism {
namespace surface {

Insolation::Insolation(std::shared_ptr<const Grid> grid) : Component(grid) {
}

void Insolation::update(double t, double dt, const array::Scalar1 &surface_elevation) {
  update_impl(t, dt, surface_elevation);
}

void Insolation::init(const InputOptions &options) {
  init_impl(options);
}

void Insolation::surface_insolation_energy_series(int i, int j,
                                                  const std::vector<OrbitalParameters> &orbital,
                                                  double dt, std::vector<double> &result) const {
  surface_insolation_energy_series_impl(i, j, orbital, dt, result);
}

DEBMSimpleInsolation::DEBMSimpleInsolation(std::shared_ptr<const Grid> grid,
                                           DEBMSimplePointwise &pointwise_model)
  : Insolation(grid),
    m_pointwise_model(pointwise_model){
  // empty
}

void DEBMSimpleInsolation::init_impl(const InputOptions &/*options*/) {
}

void DEBMSimpleInsolation::update_impl(double /*t*/, double /*dt*/,
                                       const array::Scalar1 &surface_elevation) {
  m_surface_elevation = &surface_elevation;
}

void DEBMSimpleInsolation::begin_access() const {
  m_grid->latitude().begin_access();
  m_surface_elevation->begin_access();
}

void DEBMSimpleInsolation::end_access() const {
  m_surface_elevation->end_access();
  m_grid->latitude().end_access();
}

void DEBMSimpleInsolation::surface_insolation_energy_series_impl(
    int i, int j, const std::vector<OrbitalParameters> &orbital, double dt,
    std::vector<double> &result) const {

  double surface_elevation = (*m_surface_elevation)(i, j);
  double transmissivity    = m_pointwise_model.atmosphere_transmissivity(surface_elevation);

  const auto &latitude = m_grid->latitude();

  for (unsigned k = 0; k < result.size(); ++k) {
    double toa_insolation_energy = m_pointwise_model.toa_insolation_energy(
        orbital[k].solar_declination, orbital[k].distance_factor, latitude(i, j), dt);
    result[k] = transmissivity * toa_insolation_energy;
  }
}

} // end of namespace surface
} // end of namespace pism
