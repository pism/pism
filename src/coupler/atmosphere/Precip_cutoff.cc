// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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

#include "Precip_cutoff.hh"

#include "pism/util/ConfigInterface.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace atmosphere {

Precip_cutoff::Precip_cutoff(IceGrid::ConstPtr grid, std::shared_ptr<AtmosphereModel> in)
  : AtmosphereModel(grid, in),
    m_usurf(grid, "surface_elevation", WITHOUT_GHOSTS),
    m_mask(grid, "precip_mask", WITHOUT_GHOSTS) {

  m_precipitation = allocate_precipitation(grid);

  m_option = "atmosphere.precip_cutoff";

  m_use_cutoff_height = m_config->get_flag(m_option + ".use_cutoff_height");
  m_cutoff_height = m_config->get_number(m_option + ".cutoff_height");
}

Precip_cutoff::~Precip_cutoff() {
  // empty
}

void Precip_cutoff::init_impl(const Geometry &geometry) {
  m_input_model->init(geometry);

  m_log->message(2, "* Initializing precipitation cut-off modifier...\n");

  auto filename = m_config->get_string(m_option + ".cutoff_mask_file");

  if (not filename.empty()) {
    m_mask.regrid(filename, OPTIONAL, 0);
  } else {
    m_mask.set(0);
  }

}

void Precip_cutoff::init_timeseries_impl(const std::vector<double> &ts) const {
  AtmosphereModel::init_timeseries_impl(ts);
}

void Precip_cutoff::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  m_precipitation->copy_from(m_input_model->mean_precipitation());

  m_usurf.copy_from(geometry.ice_surface_elevation);

  IceModelVec::AccessList list{&m_mask, &m_usurf,
                               m_precipitation.get()};

  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    const bool cutoff = (m_mask(i, j) == 1) or ( m_use_cutoff_height and (m_usurf(i, j) > m_cutoff_height) );
    if (cutoff) {
      (*m_precipitation)(i, j) = 0.0;
    }
  }
}

const IceModelVec2S& Precip_cutoff::mean_precipitation_impl() const {
  return *m_precipitation;
}

void Precip_cutoff::precip_time_series_impl(int i, int j, std::vector<double> &result) const {
  m_input_model->precip_time_series(i, j, result);

  const bool cutoff = (m_mask(i, j) == 1) or ( m_use_cutoff_height and (m_usurf(i, j) > m_cutoff_height) );
  if (cutoff) {
    for (unsigned int k = 0; k < result.size(); ++k) {
      result[k] = 0.0;
    }
  }
}

} // end of namespace atmosphere
} // end of namespace pism
