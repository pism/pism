/* Copyright (C) 2019, 2021, 2022 PISM Authors
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

#include "PrescribedRetreat.hh"
#include "pism/coupler/util/options.hh"

namespace pism {

PrescribedRetreat::PrescribedRetreat(IceGrid::ConstPtr grid)
  : Component(grid) {
  ForcingOptions opt(*m_grid->ctx(), "geometry.front_retreat.prescribed");
  {
    unsigned int buffer_size = m_config->get_number("input.forcing.buffer_size");

    File file(m_grid->com, opt.filename, PISM_NETCDF3, PISM_READONLY);

    m_retreat_mask = std::make_shared<array::Forcing>(m_grid,
                                                      file,
                                                      "land_ice_area_fraction_retreat",
                                                      "", // no standard name
                                                      buffer_size,
                                                      opt.periodic);
    m_retreat_mask->set_attrs("forcing", "maximum ice extent mask",
                              "1", "1", "", 0);
  }
}

void PrescribedRetreat::init() {

  ForcingOptions opt(*m_grid->ctx(), "geometry.front_retreat.prescribed");

  m_log->message(2,
                 "* Initializing the prescribed front retreat mechanism\n"
                 "  using a time-dependent ice extent mask '%s' in '%s'...\n",
                 m_retreat_mask->get_name().c_str(),
                 opt.filename.c_str());

  m_retreat_mask->init(opt.filename, opt.periodic);
}

void PrescribedRetreat::update(double t,
                               double dt,
                               array::Scalar& ice_thickness,
                               array::Scalar& ice_area_specific_volume) {
  m_retreat_mask->update(t, dt);
  m_retreat_mask->average(t, dt);

  array::AccessScope list{m_retreat_mask.get(), &ice_thickness, &ice_area_specific_volume};

  for (auto p = m_grid->points(); p; p.next()) {
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

MaxTimestep PrescribedRetreat::max_timestep_impl(double t) const {
  auto dt = m_retreat_mask->max_timestep(t);

  if (dt.finite()) {
    return MaxTimestep(dt.value(), "prescribed ice retreat");
  } else {
    return MaxTimestep("prescribed ice retreat");
  }
}

} // end of namespace pism
