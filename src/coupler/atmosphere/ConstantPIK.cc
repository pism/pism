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

#include "ConstantPIK.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Vars.hh"
#include "pism/util/IceGrid.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/MaxTimestep.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {
namespace atmosphere {

PIK::PIK(IceGrid::ConstPtr g)
  : AtmosphereModel(g) {

  m_precipitation.create(m_grid, "precipitation", WITHOUT_GHOSTS);
  m_precipitation.set_attrs("model_state", "precipitation rate",
                            "kg m-2 second-1",
                            "precipitation_flux", 0);
  m_precipitation.metadata(0).set_string("glaciological_units", "kg m-2 year-1");
  m_precipitation.set_time_independent(true);

  m_air_temp.create(m_grid, "air_temp", WITHOUT_GHOSTS);
  m_air_temp.set_attrs("model_state", "mean annual near-surface air temperature",
                           "Kelvin", "", 0);
  m_air_temp.set_time_independent(true);
}

const IceModelVec2S& PIK::mean_precipitation_impl() const {
  return m_precipitation;
}

const IceModelVec2S& PIK::mean_annual_temp_impl() const {
  return m_air_temp;
}

void PIK::begin_pointwise_access_impl() const {
  m_precipitation.begin_access();
  m_air_temp.begin_access();
}

void PIK::end_pointwise_access_impl() const {
  m_precipitation.end_access();
  m_air_temp.end_access();
}

void PIK::temp_time_series_impl(int i, int j, std::vector<double> &result) const {
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = m_air_temp(i,j);
  }
}

void PIK::precip_time_series_impl(int i, int j, std::vector<double> &result) const {
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = m_precipitation(i,j);
  }
}

void PIK::define_model_state_impl(const PIO &output) const {
  m_precipitation.define(output);
}

void PIK::write_model_state_impl(const PIO &output) const {
  m_precipitation.write(output);
}

void PIK::init_impl(const Geometry &geometry) {
  (void) geometry;

  m_log->message(2,
                 "* Initializing the constant-in-time atmosphere model PIK.\n"
                 "  It reads a precipitation field directly from the file and holds it constant.\n"
                 "  Near-surface air temperature is parameterized as in Martin et al. 2011, Eqn. 2.0.2.\n");

  InputOptions opts = process_input_options(m_grid->com, m_config);

  // read snow precipitation rate and air_temps from file
  m_log->message(2,
                 "    reading mean annual ice-equivalent precipitation rate 'precipitation'\n"
                 "    from %s ... \n",
                 opts.filename.c_str());
  if (opts.type == INIT_BOOTSTRAP) {
    m_precipitation.regrid(opts.filename, CRITICAL);
  } else {
    m_precipitation.read(opts.filename, opts.record); // fails if not found!
  }
}

MaxTimestep PIK::max_timestep_impl(double t) const {
  (void) t;
  return MaxTimestep("atmosphere PIK");
}

void PIK::update_impl(const Geometry &geometry, double, double) {
  // Compute near-surface air temperature using a latitude- and
  // elevation-dependent parameterization:

  const IceModelVec2S
    &elevation = geometry.ice_surface_elevation,
    &latitude  = geometry.latitude;

  IceModelVec::AccessList list{&m_air_temp, &elevation, &latitude};
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    m_air_temp(i, j) = 273.15 + 30 - 0.0075 * elevation(i, j) - 0.68775 * latitude(i, j)*(-1.0) ;
  }
}

void PIK::init_timeseries_impl(const std::vector<double> &ts) const {
  m_ts_times = ts;
}


} // end of namespace atmosphere
} // end of namespace pism
