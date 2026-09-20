// Copyright (C) 2026 Andy Aschwanden
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include <algorithm>

#include "pism/coupler/ocean/Plume.hh"
#include "pism/coupler/ocean/PlumePhysics.hh"
#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Config.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Logger.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/util/Time.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"

namespace pism {
namespace ocean {

Plume::Plume(std::shared_ptr<const Grid> grid)
  : PlumeModel(grid),
    m_ambient_temperature(grid, "plume_temperature"),
    m_ambient_salinity(grid, "plume_salinity") {

  ForcingOptions opt(*m_grid->ctx(), "ocean.plume");

  {
    auto buffer_size = static_cast<int>(m_config->get_number("input.forcing.buffer_size"));

    File file(m_grid->com, opt.filename, io::PISM_NETCDF3, io::PISM_READONLY);

    m_theta_ocean = std::make_shared<array::Forcing>(m_grid,
                                                     file,
                                                     "theta_ocean",
                                                     "", // no standard name
                                                     buffer_size,
                                                     opt.periodic,
                                                     LINEAR);

    m_salinity_ocean = std::make_shared<array::Forcing>(m_grid,
                                                        file,
                                                        "salinity_ocean",
                                                        "", // no standard name
                                                        buffer_size,
                                                        opt.periodic,
                                                        LINEAR);
  }

  m_thermal_forcing = m_config->get_flag("ocean.plume.temperature_as_thermal_forcing");

  // Thermal forcing is a temperature difference, but it is read as an absolute
  // temperature so that a file in degrees Celsius picks up the same kelvin offset PICO
  // relies on; compute_ambient_temperature() undoes it.
  m_theta_ocean->metadata(0)
      .long_name(m_thermal_forcing ? "thermal forcing of the adjacent ocean"
                                   : "potential temperature of the adjacent ocean")
      .units("kelvin");

  m_salinity_ocean->metadata(0)
      .long_name("salinity of the adjacent ocean")
      .units("g/kg");

  m_ambient_temperature.metadata(0)
      .long_name("ambient ocean temperature driving the plume")
      .units("kelvin");
  m_ambient_temperature.metadata()["_FillValue"] = { 0.0 };
  m_ambient_temperature.set(0.0);

  m_ambient_salinity.metadata(0)
      .long_name("ambient ocean salinity driving the plume")
      .units("g/kg");
  m_ambient_salinity.metadata()["_FillValue"] = { 0.0 };
  m_ambient_salinity.set(0.0);
}

void Plume::init_impl(const Geometry &geometry) {

  m_log->message(2,
                 "* Initializing the buoyant plume ocean model,\n"
                 "  reading ambient ocean temperature and salinity from a file...\n");
  m_log->message(2, "  Note: the plume model requires stress balance computation to be enabled.\n");

  if (m_thermal_forcing) {
    m_log->message(2,
                   "  Interpreting 'theta_ocean' as thermal forcing (temperature above the\n"
                   "  freezing point at the grounding-line depth).\n");
  }

  ForcingOptions opt(*m_grid->ctx(), "ocean.plume");

  m_theta_ocean->init(opt.filename, opt.periodic);

  // read ocean salinity from a file if present, otherwise use a constant
  {
    File input(m_grid->com, opt.filename, io::PISM_GUESS, io::PISM_READONLY);

    auto variable_name = m_salinity_ocean->metadata().get_name();

    if (input.variable_exists(variable_name)) {
      m_salinity_ocean->init(opt.filename, opt.periodic);
    } else {
      double salinity = m_config->get_number("constants.sea_water.salinity", "g / kg");

      m_salinity_ocean = array::Forcing::Constant(m_grid, variable_name, salinity);

      m_log->message(2, "  Variable '%s' not found; using constant salinity: %f (g / kg).\n",
                     variable_name.c_str(), salinity);
    }
  }

  // read time-independent data right away:
  if (m_theta_ocean->buffer_size() == 1 and m_salinity_ocean->buffer_size() == 1) {
    m_theta_ocean->update(time().current(), 0.0);
    m_salinity_ocean->update(time().current(), 0.0);
  }

  double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                        *m_water_column_pressure);
}

void Plume::update_impl(const Inputs &inputs, double t, double dt) {

  m_theta_ocean->update(t, dt);
  m_salinity_ocean->update(t, dt);

  m_theta_ocean->average(t, dt);
  m_salinity_ocean->average(t, dt);

  m_ambient_salinity.copy_from(*m_salinity_ocean);

  PlumePhysics physics(*m_config);

  if (inputs.stress_balance == nullptr) {
    // No ice velocity to transport the grounding-line elevation with, e.g. during the
    // "fake" ocean time step of bootstrapping: no plume, no melt.
    m_log->message(3,
                   "WARNING: the plume model requires the stress balance to transport the\n"
                   "         grounding-line elevation. Stress balance not available: setting\n"
                   "         the sub-shelf melt rate to zero.\n");

    compute_shelf_base_temperature(physics, *inputs.geometry, *m_shelf_base_temperature);
    m_basal_melt_rate.set(0.0);
    m_shelf_base_mass_flux->set(0.0);

    double
      ice_density   = m_config->get_number("constants.ice.density"),
      water_density = m_config->get_number("constants.sea_water.density"),
      g             = m_config->get_number("constants.standard_gravity");

    compute_average_water_column_pressure(*inputs.geometry, ice_density, water_density, g,
                                          *m_water_column_pressure);
    return;
  }

  if (inputs.hydrology == nullptr) {
    m_log->message(3,
                   "WARNING: the plume model requires hydrology routing for subglacial discharge.\n");
  }

  m_log->message(3, "  Plume: Computing plume-based melt rates...\n");

  update_geometry(inputs);

  compute_ambient_temperature(physics, *inputs.geometry, m_ambient_temperature);

  compute_shelf_base_temperature(physics, *inputs.geometry, *m_shelf_base_temperature);

  update_melt_rate(inputs, physics, m_ambient_temperature, m_ambient_salinity);
}

//! Ambient temperature `T_a` (kelvin) on floating cells.
/*!
 * With thermal forcing as input, `T_a = T_f(S_a, z_gl) + TF` so that the plume's driving
 * temperature difference `T_a - T_f(S_a, z_gl)` is exactly the forcing read from the
 * file. `z_gl` is clamped to the local shelf base as in compute_melt_rate().
 */
void Plume::compute_ambient_temperature(const PlumePhysics &physics, const Geometry &geometry,
                                        array::Scalar &result) const {

  const double T0 = m_config->get_number("constants.fresh_water.melting_point_temperature");

  const auto &cell_type = geometry.cell_type;
  const auto &surf      = geometry.ice_surface_elevation;
  const auto &thk       = geometry.ice_thickness;

  const array::Scalar &theta    = *m_theta_ocean;
  const array::Scalar &salinity = *m_salinity_ocean;

  array::AccessScope scope{&cell_type, &surf, &thk, &m_grounding_line_elevation,
                           &theta, &salinity, &result};

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    if (not cell_type.floating_ice(i, j)) {
      result(i, j) = 0.0;
      continue;
    }

    if (m_thermal_forcing) {
      const double z_b  = surf(i, j) - thk(i, j);
      const double z_gl = std::min(m_grounding_line_elevation(i, j), z_b);

      result(i, j) = physics.characteristic_freezing_point(salinity(i, j), z_gl) + (theta(i, j) - T0);
    } else {
      result(i, j) = theta(i, j);
    }
  }
}

//! Sub-shelf ice temperature: the freezing point at the shelf base under floating ice,
//! the fresh-water melting point elsewhere.
void Plume::compute_shelf_base_temperature(const PlumePhysics &physics, const Geometry &geometry,
                                           array::Scalar &result) const {

  const double T0 = m_config->get_number("constants.fresh_water.melting_point_temperature");

  const auto &cell_type = geometry.cell_type;
  const auto &surf      = geometry.ice_surface_elevation;
  const auto &thk       = geometry.ice_thickness;

  const array::Scalar &salinity = *m_salinity_ocean;

  array::AccessScope scope{&cell_type, &surf, &thk, &salinity, &result};

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    if (cell_type.floating_ice(i, j)) {
      const double z_b = surf(i, j) - thk(i, j);
      result(i, j) = physics.characteristic_freezing_point(salinity(i, j), z_b);
    } else {
      result(i, j) = T0;
    }
  }
}

MaxTimestep Plume::max_timestep_impl(double t, const CFLData *cfl_data) const {
  (void) cfl_data;

  auto dt_theta    = m_theta_ocean->max_timestep(t);
  auto dt_salinity = m_salinity_ocean->max_timestep(t);

  if (dt_theta.finite() and dt_salinity.finite()) {
    return { std::min(dt_theta.value(), dt_salinity.value()), "ocean plume" };
  }
  if (dt_theta.finite()) {
    return { dt_theta.value(), "ocean plume" };
  }
  if (dt_salinity.finite()) {
    return { dt_salinity.value(), "ocean plume" };
  }

  return MaxTimestep("ocean plume");
}

DiagnosticList Plume::spatial_diagnostics_impl() const {
  return plume_diagnostics(m_ambient_temperature, m_ambient_salinity);
}

} // end of namespace ocean
} // end of namespace pism
