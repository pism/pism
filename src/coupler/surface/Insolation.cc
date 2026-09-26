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
#include "pism/coupler/surface/TerrainInsolation.hh"
#include "pism/coupler/surface/DEBMSimple.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/util/Time.hh"

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

DEBMTerrainShadedInsolation::DEBMTerrainShadedInsolation(
    std::shared_ptr<const Grid> grid, std::function<double(double)> atmosphere_transmissivity)
    : Insolation(grid),
      m_model(new TerrainInsolation(grid, atmosphere_transmissivity)),
      m_time_name(time().variable_name() + "_insolation"),
      m_t_last(time().current()) {

  m_update_interval = m_config->get_number("bed_deformation.update_interval", "seconds");
  m_t_eps           = m_config->get_number("time_stepping.resolution", "seconds");
}

void DEBMTerrainShadedInsolation::init_impl(const InputOptions &options) {
}

void DEBMTerrainShadedInsolation::update_impl(double t, double dt,
                                              const array::Scalar1 &surface_elevation) {

  m_model->update_shading(surface_elevation);

  m_model->update_daily_insolation(t + 0.5 * dt, surface_elevation);
}

void DEBMTerrainShadedInsolation::begin_access() const {
  m_model->insolation().begin_access();
}

void DEBMTerrainShadedInsolation::end_access() const {
  m_model->insolation().end_access();
}

void DEBMTerrainShadedInsolation::surface_insolation_energy_series_impl(
    int i, int j, const std::vector<OrbitalParameters> &orbital, double dt,
    std::vector<double> &result) const {
  // the daily-mean insolation rate (W m^-2) at this cell was computed in
  // update_insolation_input
  double rate = m_model->insolation()(i, j);
  for (size_t k = 0; k < orbital.size(); ++k) {
    result[k] = rate * dt;
  }
}

MaxTimestep DEBMTerrainShadedInsolation::max_timestep_impl(double t,
                                                           const CFLData * /*cfl_data*/) const {

  if (t < m_t_last) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "time %f is less than the previous time %f",
                                  t, m_t_last);
  }

  if (m_update_interval == 0.0) {
    return {};
  }

  // Find the smallest time of the form m_t_last + k * m_update_interval that is greater
  // than t
  double k = std::ceil((t - m_t_last) / m_update_interval);

  double t_next = m_t_last + k * m_update_interval;

  double dt_max = m_update_interval;
  if (t < t_next) {
    dt_max = t_next - t;
  }

  if (dt_max < m_t_eps) {
    dt_max = m_update_interval;
  }

  return MaxTimestep(dt_max, "bed_def");
}

std::set<VariableMetadata> DEBMTerrainShadedInsolation::state_impl() const {
  auto variables = array::metadata({ &m_model->horizon_altitude() });

  VariableMetadata T(m_time_name, m_sys);
  T.long_name("time of the last update of the terrain-shaded insolation model")
      .units(time().units())
      .set_time_dependent(true);
  T["calendar"] = time().calendar();

  variables.insert(T);

  return variables;
}

void DEBMTerrainShadedInsolation::write_state_impl(const OutputFile &output) const {
  m_model->horizon_altitude().write(output);

  auto t_length = output.time_dimension_length();
  auto t_start = t_length > 0 ? t_length - 1 : 0;

  output.write_timeseries(m_time_name, { t_start }, { 1 }, { m_t_last });
}

DiagnosticList DEBMTerrainShadedInsolation::spatial_diagnostics_impl() const {
  DiagnosticList result = { { "insolation", Diagnostic::wrap(m_model->insolation()) }};

  return result;
}


} // end of namespace surface
} // end of namespace pism
