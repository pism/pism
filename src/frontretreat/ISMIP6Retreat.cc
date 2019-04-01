/* Copyright (C) 2019 PISM Authors
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

#include "ISMIP6Retreat.hh"
#include "pism/coupler/util/options.hh"

namespace pism {

ISMIP6Retreat::ISMIP6Retreat(IceGrid::ConstPtr grid)
  : Component(grid) {
  ForcingOptions opt(*m_grid->ctx(), "geometry.front_retreat.ismip6");
  {
    unsigned int buffer_size = m_config->get_double("climate_forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");
    bool periodic = opt.period > 0;

    PIO file(m_grid->com, "netcdf3", opt.filename, PISM_READONLY);

    m_retreat_mask = IceModelVec2T::ForcingField(m_grid,
                                                 file,
                                                 "land_ice_area_fraction_retreat",
                                                 "", // no standard name
                                                 buffer_size,
                                                 evaluations_per_year,
                                                 periodic);
  }
}

ISMIP6Retreat::~ISMIP6Retreat() {
  // empty
}

void ISMIP6Retreat::init() {

  ForcingOptions opt(*m_grid->ctx(), "geometry.front_retreat.ismip6");

  m_log->message(2,
                 "* Initializing the ISMIP6 front retreat parameterization\n"
                 "  using a time-dependent ice extent mask read from '%s'...",
                 opt.filename.c_str());

  m_retreat_mask->init(opt.filename, opt.period, opt.reference_time);
}

void ISMIP6Retreat::update(double t,
                           double dt,
                           IceModelVec2S& ice_thickness,
                           IceModelVec2S& ice_area_specific_volume) {
  m_retreat_mask->update(t, dt);
  m_retreat_mask->average(t, dt);

  IceModelVec::AccessList list{m_retreat_mask.get(), &ice_thickness, &ice_area_specific_volume};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    double f = (*m_retreat_mask)(i, j);

    if (f <= 0.0) {
      ice_area_specific_volume(i, j) = 0.0;
      ice_thickness(i, j)            = 0.0;
    } else if (f < 1.0) {
      ice_area_specific_volume(i, j) = ice_thickness(i, j) * f;
      ice_thickness(i, j)            = 0.0;
    } else {
      // M == 1.0: do nothing
    }
  }
}

MaxTimestep ISMIP6Retreat::max_timestep_impl(double t) const {
  auto dt = m_retreat_mask->max_timestep(t);
  return MaxTimestep(dt.value(), "ISMIP6 retreat");
}

} // end of namespace pism
