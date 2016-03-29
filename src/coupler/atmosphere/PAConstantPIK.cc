// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016 PISM Authors
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

#include <gsl/gsl_math.h>

#include "PAConstantPIK.hh"
#include "base/util/pism_utilities.hh"
#include "base/util/PISMVars.hh"
#include "base/util/IceGrid.hh"
#include "base/util/PISMConfigInterface.hh"
#include "base/util/pism_const.hh"
#include "base/util/io/io_helpers.hh"
#include "base/util/MaxTimestep.hh"

namespace pism {
namespace atmosphere {

PIK::PIK(IceGrid::ConstPtr g)
  : AtmosphereModel(g),
    m_air_temp_snapshot(m_sys, "air_temp_snapshot") {

  // allocate IceModelVecs for storing temperature and precipitation fields:

  // create mean annual ice equivalent precipitation rate (before separating
  // rain, and before melt, etc. in SurfaceModel)
  m_precipitation.create(m_grid, "precipitation", WITHOUT_GHOSTS);
  m_precipitation.set_attrs("climate_state",
                          "mean annual ice-equivalent precipitation rate",
                          "m s-1",
                          ""); // no CF standard_name ??
  m_precipitation.metadata().set_string("glaciological_units", "m year-1");
  m_precipitation.write_in_glaciological_units = true;
  m_precipitation.set_time_independent(true);

  m_air_temp.create(m_grid, "air_temp", WITHOUT_GHOSTS);
  m_air_temp.set_attrs("climate_state",
                     "mean annual near-surface (2 m) air temperature",
                     "K",
                     "");
  m_air_temp.set_time_independent(true);

  // initialize metadata for "air_temp_snapshot"
  m_air_temp_snapshot.set_string("pism_intent", "diagnostic");
  m_air_temp_snapshot.set_string("long_name",
                               "snapshot of the near-surface air temperature");
  m_air_temp_snapshot.set_string("units", "K");
}

void PIK::mean_precipitation(IceModelVec2S &result) {
  result.copy_from(m_precipitation);
}

void PIK::mean_annual_temp(IceModelVec2S &result) {
  result.copy_from(m_air_temp);
}

void PIK::begin_pointwise_access() {
  m_precipitation.begin_access();
  m_air_temp.begin_access();
}

void PIK::end_pointwise_access() {
  m_precipitation.end_access();
  m_air_temp.end_access();
}

void PIK::temp_time_series(int i, int j, std::vector<double> &result) {
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = m_air_temp(i,j);
  }
}

void PIK::precip_time_series(int i, int j, std::vector<double> &result) {
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = m_precipitation(i,j);
  }
}

void PIK::temp_snapshot(IceModelVec2S &result) {
  mean_annual_temp(result);
}

void PIK::add_vars_to_output_impl(const std::string &keyword, std::set<std::string> &result) {
  result.insert("precipitation");
  result.insert("air_temp");

  if (keyword == "big" || keyword == "2dbig") {
    result.insert("air_temp_snapshot");
  }
}

void PIK::define_variables_impl(const std::set<std::string> &vars, const PIO &nc,
                                            IO_Type nctype) {
  if (set_contains(vars, "air_temp_snapshot")) {
    std::string order = m_config->get_string("output_variable_order");
    io::define_spatial_variable(m_air_temp_snapshot, *m_grid, nc, nctype, order, false);
  }

  if (set_contains(vars, "precipitation")) {
    m_precipitation.define(nc, nctype);
  }

  if (set_contains(vars, "air_temp")) {
    m_air_temp.define(nc, nctype);
  }
}

void PIK::write_variables_impl(const std::set<std::string> &vars, const PIO &nc) {
  if (set_contains(vars, "air_temp_snapshot")) {
    IceModelVec2S tmp;
    tmp.create(m_grid, "air_temp_snapshot", WITHOUT_GHOSTS);
    tmp.metadata() = m_air_temp_snapshot;

    temp_snapshot(tmp);

    tmp.write(nc);
  }

  if (set_contains(vars, "precipitation")) {
    m_precipitation.write(nc);
  }

  if (set_contains(vars, "air_temp")) {
    m_air_temp.write(nc);
  }
}

void PIK::init() {
  bool do_regrid = false;
  int start = -1;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  m_log->message(2,
             "* Initializing the constant-in-time atmosphere model PIK.\n"
             "  It reads a precipitation field directly from the file and holds it constant.\n"
             "  Near-surface air temperature is parameterized as in Martin et al. 2011, Eqn. 2.0.2.\n");

  // find PISM input file to read data from:
  find_pism_input(m_input_file, do_regrid, start);

  // read snow precipitation rate and air_temps from file
  m_log->message(2,
             "    reading mean annual ice-equivalent precipitation rate 'precipitation'\n"
             "    from %s ... \n",
             m_input_file.c_str());
  if (do_regrid) {
    m_precipitation.regrid(m_input_file, CRITICAL);
  } else {
    m_precipitation.read(m_input_file, start); // fails if not found!
  }
}

MaxTimestep PIK::max_timestep_impl(double t) {
  (void) t;
  return MaxTimestep();
}

void PIK::update_impl(double, double) {
  // Compute near-surface air temperature using a latitude- and
  // elevation-dependent parameterization:

  const IceModelVec2S
    *usurf = m_grid->variables().get_2d_scalar("surface_altitude"),
    *lat   = m_grid->variables().get_2d_scalar("latitude");

  IceModelVec::AccessList list;
  list.add(m_air_temp);
  list.add(*usurf);
  list.add(*lat);
  for (Points p(*m_grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    m_air_temp(i, j) = 273.15 + 30 - 0.0075 * (*usurf)(i,j) - 0.68775 * (*lat)(i,j)*(-1.0) ;
  }
}

void PIK::init_timeseries(const std::vector<double> &ts) {
  m_ts_times = ts;
}


} // end of namespace atmosphere
} // end of namespace pism
