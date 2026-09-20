// Copyright (C) 2012-2019, 2021, 2022, 2023, 2024, 2025, 2026 Constantine Khrulev, Ricarda Winkelmann, Ronja Reese, Torsten
// Albrecht, Matthias Mengel, and Andy Aschwanden
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
//
// Please cite this model as:
// 1.
// Antarctic sub-shelf melt rates via PICO
// R. Reese, T. Albrecht, M. Mengel, X. Asay-Davis and R. Winkelmann
// The Cryosphere, 12, 1969-1985, (2018)
// DOI: 10.5194/tc-12-1969-2018
//
// 2.
// A box model of circulation and melting in ice shelf caverns
// D. Olbers & H. Hellmer
// Ocean Dynamics (2010), Volume 60, Issue 1, pp 141–153
// DOI: 10.1007/s10236-009-0252-z
//
// 3.
// PICOP, a new ocean melt parameterization under ice shelves
// combining PICO and a plume model.
// T. Pelle, M. Morlighem, J.H. Bondzio
// The Cryosphere, 13, 1043-49, (2019)
// DOI: 10.5194/tc-13-1043-2019


#include "pism/coupler/ocean/Picop.hh"
#include "pism/coupler/ocean/PlumePhysics.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/Config.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Logger.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {

namespace ocean {

Picop::Picop(std::shared_ptr<const Grid> grid)
  : PlumeModel(grid),
    m_pico(std::make_shared<Pico>(grid)),
    m_theta_ocean(m_pico->get_temperature()),
    m_salinity_ocean(m_pico->get_salinity()) {
  // empty
}

void Picop::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_pico->init(geometry);
  m_log->message(2, "* Initializing the Plume extension of PICO (PICOP) for the ocean ...\n");
  m_log->message(2, "  Note: PICOP requires stress balance computation to be enabled.\n");

  double
    ice_density   = m_config->get_number("constants.ice.density"),
    water_density = m_config->get_number("constants.sea_water.density"),
    g             = m_config->get_number("constants.standard_gravity");

  compute_average_water_column_pressure(geometry, ice_density, water_density, g,
                                        *m_water_column_pressure);
}

std::set<VariableMetadata> Picop::state_impl() const {
  return m_pico->state();
}

void Picop::write_state_impl(const OutputFile &output) const {
  m_pico->write_state(output);
}

void Picop::update_impl(const Inputs &inputs, double t, double dt) {

  m_pico->update(inputs, t, dt);

  if (inputs.stress_balance == nullptr) {
    // Use outputs from PICO if the stress balance is not available
    m_log->message(3,
                   "WARNING: PICOP requires stress balance for plume transport calculations.\n"
                   "         Stress balance not available - falling back to PICO melt rates.\n"
                   "         To use PICOP, enable stress balance computation.\n");
    m_shelf_base_temperature->copy_from(m_pico->shelf_base_temperature());
    m_shelf_base_mass_flux->copy_from(m_pico->shelf_base_mass_flux());
    return;
  }

  if (inputs.hydrology == nullptr) {
    // Use outputs from PICOP if the stress balance is not available
    m_log->message(3,
                   "WARNING: PICOP requires hydrology routing for plume transport calculations.\n"
                   );
  }

  m_log->message(3, "  PICOP: Computing plume-based melt rates...\n");

  PlumePhysics picop_physics(*m_config);

  update_geometry(inputs);

  m_shelf_base_temperature->copy_from(m_pico->shelf_base_temperature());

  update_melt_rate(inputs, picop_physics, m_theta_ocean, m_salinity_ocean);
}

MaxTimestep Picop::max_timestep_impl(double t, const CFLData *cfl_data) const {

  auto pico_dt_max = m_pico->max_timestep(t, cfl_data);
  if (pico_dt_max.finite()) {
    return { pico_dt_max.value(), "ocean picop" };
  }

  return MaxTimestep("ocean picop");
}

DiagnosticList Picop::spatial_diagnostics_impl() const {
  return plume_diagnostics(m_theta_ocean, m_salinity_ocean);
}

} // end of namespace ocean
} // end of namespace pism
