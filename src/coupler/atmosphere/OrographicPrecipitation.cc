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

#include "OrographicPrecipitation.hh"

#include "pism/coupler/util/options.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Time.hh"
#include "pism/geometry/Geometry.hh"
#include "OrographicPrecipitationSerial.hh"

namespace pism {
namespace atmosphere {

OrographicPrecipitation::OrographicPrecipitation(IceGrid::ConstPtr g)
  : AtmosphereModel(g, nullptr), m_surface(g, "surface_altitude", WITHOUT_GHOSTS) {
  ForcingOptions opt(*m_grid->ctx(), "atmosphere.orographic_precipitation");
  
  {

    unsigned int buffer_size = m_config->get_double("climate_forcing.buffer_size");
    unsigned int evaluations_per_year = m_config->get_double("climate_forcing.evaluations_per_year");
    bool periodic = opt.period > 0;

    PIO file(m_grid->com, "netcdf3", opt.filename, PISM_READONLY);

    m_air_temp = IceModelVec2T::ForcingField(m_grid,
                                             file,
                                             "air_temp",
                                             "", // no standard name
                                             buffer_size,
                                             evaluations_per_year,
                                             periodic);

  }

  {
    m_air_temp->set_attrs("diagnostic", "mean annual near-surface air temperature",
                          "Kelvin", "", 0);
    m_air_temp->metadata(0).set_double("valid_min", 0.0);
    m_air_temp->metadata(0).set_double("valid_max", 323.15); // 50 C
  }

  m_work0 = m_surface.allocate_proc0_copy();

  const int
    Mx = m_grid->Mx(),
    My = m_grid->My(),
    Z  = m_config->get_double("bed_deformation.lc.grid_size_factor"),
    Nx = Z*(Mx - 1) + 1,
    Ny = Z*(My - 1) + 1;

  const double
    Lx = Z * (m_grid->x0() - m_grid->x(0)),
    Ly = Z * (m_grid->y0() - m_grid->y(0));

  m_extended_grid = IceGrid::Shallow(m_grid->ctx(),
                                     Lx, Ly,
                                     m_grid->x0(), m_grid->y0(),
                                     Nx, Ny, CELL_CORNER, NOT_PERIODIC);
}

OrographicPrecipitation::~OrographicPrecipitation() {
  // empty
}

void OrographicPrecipitation::init_impl(const Geometry &geometry) {
  m_log->message(2,
             "* Initializing the atmosphere model computing precipitation using the\n"
             "  Linear Theory of Orographic Precipitation model\n"
             "  and reading near-surface air temperature from a file...\n");

  m_reference =
    "R. B. Smith and I. Barstad, 2004. "
    "A Linear Theory of Orographic Precipitation. J. Atmos. Sci. 61, 1377-1391.";

  ForcingOptions opt(*m_grid->ctx(), "atmosphere.orographic_precipitation");

  m_air_temp->init(opt.filename, opt.period, opt.reference_time);

  // read time-independent data right away:
  if (m_air_temp->n_records() == 1) {
    update(geometry, m_grid->ctx()->time()->current(), 0); // dt is irrelevant
  }
}

void OrographicPrecipitation::update_impl(const Geometry &geometry, double t, double dt) {

  const IceModelVec2S
    &h        = geometry.ice_surface_elevation;
  
  m_air_temp->update(t, dt);

  m_air_temp->average(t, dt);
}

const IceModelVec2S& OrographicPrecipitation::mean_annual_temp_impl() const {
  return *m_air_temp;
}

void OrographicPrecipitation::begin_pointwise_access_impl() const {

  m_air_temp->begin_access();
}

void OrographicPrecipitation::end_pointwise_access_impl() const {

  m_air_temp->end_access();
}

void OrographicPrecipitation::temp_time_series_impl(int i, int j, std::vector<double> &result) const {

  m_air_temp->interp(i, j, result);
}

void OrographicPrecipitation::init_timeseries_impl(const std::vector<double> &ts) const {

  m_air_temp->init_interpolation(ts);

  m_ts_times = ts;
}


} // end of namespace atmosphere
} // end of namespace pism
