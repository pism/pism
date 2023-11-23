// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2023 PISM Authors
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

#include "pism/coupler/atmosphere/OrographicPrecipitation.hh"

#include "pism/coupler/atmosphere/OrographicPrecipitationSerial.hh"
#include "pism/coupler/util/options.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Grid.hh"

namespace pism {
namespace atmosphere {

OrographicPrecipitation::OrographicPrecipitation(std::shared_ptr<const Grid> grid,
                                                 std::shared_ptr<AtmosphereModel> in)
    : AtmosphereModel(grid, in) {

  m_precipitation = allocate_precipitation(grid);

  m_work0 = m_precipitation->allocate_proc0_copy();

  const int
    Mx = m_grid->Mx(),
    My = m_grid->My(),
    Z  = m_config->get_number("atmosphere.orographic_precipitation.grid_size_factor"),
    Nx = m_grid->periodicity() & grid::X_PERIODIC ? Mx : Z * (Mx - 1) + 1,
    Ny = m_grid->periodicity() & grid::Y_PERIODIC ? My : Z * (My - 1) + 1;

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) {
      m_serial_model.reset(new OrographicPrecipitationSerial(*m_config,
                                                             Mx, My,
                                                             m_grid->dx(), m_grid->dy(),
                                                             Nx, Ny));
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();
}

OrographicPrecipitation::~OrographicPrecipitation() {
  // empty
}

const array::Scalar &OrographicPrecipitation::precipitation_impl() const {
  return *m_precipitation;
}

void OrographicPrecipitation::init_impl(const Geometry &geometry) {
  (void)geometry;

  m_input_model->init(geometry);

  m_log->message(2, "* Initializing the atmosphere model computing precipitation using the\n"
                    "  Linear Theory of Orographic Precipitation model with scalar wind speeds...\n");

  m_reference = "R. B. Smith and I. Barstad, 2004.\n"
                "A Linear Theory of Orographic Precipitation. J. Atmos. Sci. 61, 1377-1391.";

  m_precipitation->metadata()["source"] = m_reference;
}


void OrographicPrecipitation::update_impl(const Geometry &geometry, double t, double dt) {
  m_input_model->update(geometry, t, dt);

  geometry.ice_surface_elevation.put_on_proc0(*m_work0);

  ParallelSection rank0(m_grid->com);
  try {
    if (m_grid->rank() == 0) { // processor zero updates the precipitation
      m_serial_model->update(*m_work0);

      PetscErrorCode ierr = VecCopy(m_serial_model->precipitation(), *m_work0);
      PISM_CHK(ierr, "VecCopy");
    }
  } catch (...) {
    rank0.failed();
  }
  rank0.check();

  m_precipitation->get_from_proc0(*m_work0);

  // convert from mm/s to kg / (m^2 s):
  double water_density = m_config->get_number("constants.fresh_water.density");
  m_precipitation->scale(1e-3 * water_density);
}

void OrographicPrecipitation::precip_time_series_impl(int i, int j,
                                                      std::vector<double> &result) const {

  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = (*m_precipitation)(i, j);
  }
}

void OrographicPrecipitation::begin_pointwise_access_impl() const {
  m_input_model->begin_pointwise_access();
  m_precipitation->begin_access();
}

void OrographicPrecipitation::end_pointwise_access_impl() const {
  m_precipitation->end_access();
  m_input_model->end_pointwise_access();
}

} // end of namespace atmosphere
} // end of namespace pism
