// Copyright (C) 2012, 2013, 2014 PISM Authors
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

#include "PACosineYearlyCycle.hh"
#include "Timeseries.hh"
#include "PISMTime.hh"
#include "pism_options.hh"
#include "PISMConfig.hh"


PACosineYearlyCycle::PACosineYearlyCycle(IceGrid &g, const PISMConfig &conf)
  : PAYearlyCycle(g, conf), A(NULL) {
}

PACosineYearlyCycle::~PACosineYearlyCycle() {
  if (A != NULL)
    delete A;
}

PetscErrorCode PACosineYearlyCycle::init(PISMVars &vars) {
  PetscErrorCode ierr;
  bool input_file_flag, scaling_flag;
  std::string input_file, scaling_file;

  m_t = m_dt = GSL_NAN;  // every re-init restarts the clock

  variables = &vars;

  ierr = verbPrintf(2, grid.com,
		    "* Initializing the 'cosine yearly cycle' atmosphere model (-atmosphere yearly_cycle)...\n");
  CHKERRQ(ierr);

  ierr = PetscOptionsBegin(grid.com, "", "Options controlling '-atmosphere yearly_cycle'",
                           ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsString("-atmosphere_yearly_cycle_file",
                             "PACosineYearlyCycle input file name",
                             input_file, input_file_flag); CHKERRQ(ierr);
    ierr = PISMOptionsString("-atmosphere_yearly_cycle_scaling_file",
                             "PACosineYearlyCycle amplitude scaling input file name",
                             scaling_file, scaling_flag); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (input_file_flag == false) {
    PetscPrintf(grid.com,
                "PISM ERROR: Please specify an '-atmosphere yearly_cycle' input file\n"
                "            using the -atmosphere_yearly_cycle_file option.\n");
    PISMEnd();
  }

  ierr = verbPrintf(2, grid.com,
                    "  Reading mean annual air temperature, mean July air temperature, and\n"
                    "  precipitation fields from '%s'...\n", input_file.c_str()); CHKERRQ(ierr);

  ierr = air_temp_mean_annual.regrid(input_file, CRITICAL); CHKERRQ(ierr);
  ierr = air_temp_mean_july.regrid(input_file, CRITICAL); CHKERRQ(ierr);
  ierr = precipitation.regrid(input_file, CRITICAL); CHKERRQ(ierr);

  if (scaling_flag) {

    if (A == NULL) {
      A = new Timeseries(&grid, "amplitude_scaling",
                         config.get_string("time_dimension_name"));
      A->set_units("1", "1");
      A->set_dimension_units(grid.time->units_string(), "");
      A->set_attr("long_name", "cosine yearly cycle amplitude scaling");
    }

    ierr = verbPrintf(2, grid.com,
                      "  Reading cosine yearly cycle amplitude scaling from '%s'...\n",
                      scaling_file.c_str()); CHKERRQ(ierr);

    PIO nc(grid, "netcdf3");    // OK to use netcdf3
    ierr = nc.open(scaling_file, PISM_NOWRITE); CHKERRQ(ierr);
    {
      ierr = A->read(nc, grid.time); CHKERRQ(ierr);
    }
    ierr = nc.close(); CHKERRQ(ierr);

  } else {
    if (A != NULL)
      delete A;
    A = NULL;
  }

  return 0;
}


PetscErrorCode PACosineYearlyCycle::update(double my_t, double my_dt) {
  m_t = my_t;
  m_dt = my_dt;
  return 0;
}

PetscErrorCode PACosineYearlyCycle::temp_snapshot(IceModelVec2S &result) {
  PetscErrorCode ierr;

  const double
    julyday_fraction = grid.time->day_of_the_year_to_day_fraction(snow_temp_july_day),
    T                = grid.time->year_fraction(m_t + 0.5 * m_dt) - julyday_fraction,
    cos_T            = cos(2.0 * M_PI * T);

  double scaling = 1.0;
  if (A != NULL) {
    scaling = (*A)(m_t + 0.5 * m_dt);
  }

  ierr = result.begin_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.begin_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_july.begin_access(); CHKERRQ(ierr);

  for (int   i = grid.xs; i < grid.xs+grid.xm; ++i) {
    for (int j = grid.ys; j < grid.ys+grid.ym; ++j) {
      result(i,j) = air_temp_mean_annual(i,j) + (air_temp_mean_july(i,j) - air_temp_mean_annual(i,j)) * scaling * cos_T;
    }
  }

  ierr = air_temp_mean_july.end_access(); CHKERRQ(ierr);
  ierr = air_temp_mean_annual.end_access(); CHKERRQ(ierr);
  ierr = result.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode PACosineYearlyCycle::init_timeseries(double *ts, unsigned int N) {
  PetscErrorCode ierr;

  ierr = PAYearlyCycle::init_timeseries(ts, N); CHKERRQ(ierr);

  if (A != NULL) {
    for (unsigned int k = 0; k < N; ++k)
      m_cosine_cycle[k] *= (*A)(ts[k]);
  }

  return 0;
}
