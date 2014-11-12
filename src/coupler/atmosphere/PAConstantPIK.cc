// Copyright (C) 2011, 2012, 2013, 2014 PISM Authors
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

#include "PAConstantPIK.hh"
#include "PISMVars.hh"
#include "IceGrid.hh"

#include <stdexcept>

namespace pism {

PAConstantPIK::PAConstantPIK(IceGrid &g, const Config &conf)
  : AtmosphereModel(g, conf),
    air_temp_snapshot(g.get_unit_system(), "air_temp_snapshot", g) {
  PetscErrorCode ierr = allocate_PAConstantPIK(); CHKERRCONTINUE(ierr);
  if (ierr != 0) {
    throw std::runtime_error("PAConstantPIK allocation failed");
  }
}

PetscErrorCode PAConstantPIK::allocate_PAConstantPIK() {
  PetscErrorCode ierr;
  // allocate IceModelVecs for storing temperature and precipitation fields:

  // create mean annual ice equivalent precipitation rate (before separating
  // rain, and before melt, etc. in SurfaceModel)
  ierr = precipitation.create(grid, "precipitation", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = precipitation.set_attrs("climate_state",
                                 "mean annual ice-equivalent precipitation rate",
                                 "m s-1",
                                 ""); CHKERRQ(ierr); // no CF standard_name ??
  ierr = precipitation.set_glaciological_units("m year-1"); CHKERRQ(ierr);
  precipitation.write_in_glaciological_units = true;
  precipitation.set_time_independent(true);

  ierr = air_temp.create(grid, "air_temp", WITHOUT_GHOSTS); CHKERRQ(ierr);
  ierr = air_temp.set_attrs("climate_state",
                                   "mean annual near-surface (2 m) air temperature",
                                   "K",
                                   ""); CHKERRQ(ierr);
  air_temp.set_time_independent(true);

  // initialize metadata for "air_temp_snapshot"
  air_temp_snapshot.set_string("pism_intent", "diagnostic");
  air_temp_snapshot.set_string("long_name",
                               "snapshot of the near-surface air temperature");
  ierr = air_temp_snapshot.set_units("K"); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PAConstantPIK::mean_precipitation(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = precipitation.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAConstantPIK::mean_annual_temp(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = air_temp.copy_to(result); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAConstantPIK::begin_pointwise_access() {
  PetscErrorCode ierr;
  ierr = precipitation.begin_access(); CHKERRQ(ierr);
  ierr = air_temp.begin_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAConstantPIK::end_pointwise_access() {
  PetscErrorCode ierr;
  ierr = precipitation.end_access(); CHKERRQ(ierr);
  ierr = air_temp.end_access(); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode PAConstantPIK::temp_time_series(int i, int j, std::vector<double> &result) {
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = air_temp(i,j);
  }
  return 0;
}

PetscErrorCode PAConstantPIK::precip_time_series(int i, int j, std::vector<double> &result) {
  for (unsigned int k = 0; k < m_ts_times.size(); k++) {
    result[k] = precipitation(i,j);
  }
  return 0;
}

PetscErrorCode PAConstantPIK::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr;

  ierr = mean_annual_temp(result); CHKERRQ(ierr);

  return 0;
}

void PAConstantPIK::add_vars_to_output(const std::string &keyword, std::set<std::string> &result) {
  result.insert("precipitation");
  result.insert("air_temp");

  if (keyword == "big") {
    result.insert("air_temp_snapshot");
  }
}

PetscErrorCode PAConstantPIK::define_variables(const std::set<std::string> &vars, const PIO &nc,
                                            IO_Type nctype) {
  PetscErrorCode ierr;

  if (set_contains(vars, "air_temp_snapshot")) {
    ierr = air_temp_snapshot.define(nc, nctype, false); CHKERRQ(ierr);
  }

  if (set_contains(vars, "precipitation")) {
    ierr = precipitation.define(nc, nctype); CHKERRQ(ierr);
  }

  if (set_contains(vars, "air_temp")) {
    ierr = air_temp.define(nc, nctype); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAConstantPIK::write_variables(const std::set<std::string> &vars, const PIO &nc) {
  PetscErrorCode ierr;

  if (set_contains(vars, "air_temp_snapshot")) {
    IceModelVec2S tmp;
    ierr = tmp.create(grid, "air_temp_snapshot", WITHOUT_GHOSTS); CHKERRQ(ierr);
    tmp.metadata() = air_temp_snapshot;

    ierr = temp_snapshot(tmp); CHKERRQ(ierr);

    ierr = tmp.write(nc); CHKERRQ(ierr);
  }

  if (set_contains(vars, "precipitation")) {
    ierr = precipitation.write(nc); CHKERRQ(ierr);
  }

  if (set_contains(vars, "air_temp")) {
    ierr = air_temp.write(nc); CHKERRQ(ierr);
  }

  return 0;
}

PetscErrorCode PAConstantPIK::init(Vars &vars) {
  PetscErrorCode ierr;
  bool do_regrid = false;
  int start = -1;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  ierr = verbPrintf(2, grid.com,
     "* Initializing the constant-in-time atmosphere model PAConstantPIK.\n"
     "  It reads a precipitation field directly from the file and holds it constant.\n"
     "  Near-surface air temperature is parameterized as in Martin et al. 2011, Eqn. 2.0.2.\n"); CHKERRQ(ierr);

  // find PISM input file to read data from:
  ierr = find_pism_input(input_file, do_regrid, start); CHKERRQ(ierr);

  // read snow precipitation rate and air_temps from file
  ierr = verbPrintf(2, grid.com,
                    "    reading mean annual ice-equivalent precipitation rate 'precipitation'\n"
                    "    from %s ... \n",
                    input_file.c_str()); CHKERRQ(ierr); 
  if (do_regrid) {
    ierr = precipitation.regrid(input_file, CRITICAL); CHKERRQ(ierr); // fails if not found!
  } else {
    ierr = precipitation.read(input_file, start); CHKERRQ(ierr); // fails if not found!
  }

  usurf = vars.get_2d_scalar("surface_altitude");
  lat   = vars.get_2d_scalar("latitude");

  return 0;
}

PetscErrorCode PAConstantPIK::update(double, double) {
  // Compute near-surface air temperature using a latitude- and
  // elevation-dependent parameterization:

  IceModelVec::AccessList list;
  list.add(air_temp);
  list.add(*usurf);
  list.add(*lat);
  for (Points p(grid); p; p.next()) {
    const int i = p.i(), j = p.j();

    air_temp(i, j) = 273.15 + 30 - 0.0075 * (*usurf)(i,j) - 0.68775 * (*lat)(i,j)*(-1.0) ;
  }

  return 0;
}

PetscErrorCode PAConstantPIK::init_timeseries(const std::vector<double> &ts) {
  m_ts_times = ts;
  return 0;
}


} // end of namespace pism
