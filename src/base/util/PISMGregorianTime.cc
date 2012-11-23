// Copyright (C) 2012 PISM Authors
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

#include "PISMGregorianTime.hh"
#include "pism_options.hh"
#include "PIO.hh"

PISMGregorianTime::PISMGregorianTime(MPI_Comm c, const NCConfigVariable &conf)
  : PISMTime(c, conf) {

  calendar_string = "gregorian";  // only "gregorian" is supported by this class
}

PetscErrorCode PISMGregorianTime::init() {
  PetscErrorCode ierr;
  string time_file;
  bool flag;

  ierr = PISMTime::init(); CHKERRQ(ierr);

  ierr = PISMOptionsString("-time_file", "Reads time information from a file",
                           time_file, flag); CHKERRQ(ierr);

  if (flag == true) {
    ierr = verbPrintf(2, com,
                      "* Setting time from '%s'...\n",
                      time_file.c_str()); CHKERRQ(ierr);

    ierr = ignore_option(com, "-y"); CHKERRQ(ierr);
    ierr = ignore_option(com, "-ys"); CHKERRQ(ierr);
    ierr = ignore_option(com, "-ye"); CHKERRQ(ierr);
    ierr = ignore_option(com, "-reference_date"); CHKERRQ(ierr);

    ierr = init_from_file(time_file); CHKERRQ(ierr);

  }

  // initialize the units object:
  ierr = utScan(this->units().c_str(), &ut_units); CHKERRQ(ierr);

  return 0;
}

//! \brief Sets the time from a NetCDF with forcing data.
/*!
 * This allows running PISM for the duration of the available forcing.
 */
PetscErrorCode PISMGregorianTime::init_from_file(string filename) {
  PetscErrorCode ierr;
  NCTimeseries time_axis;
  NCTimeBounds bounds;
  PetscMPIInt rank;
  vector<double> time, time_bounds;
  string time_units, time_bounds_name,
    time_name = config.get_string("time_dimension_name");
  bool exists;

  ierr = MPI_Comm_rank(com, &rank); CHKERRQ(ierr);
  PIO nc(com, rank, "netcdf3"); // OK to use netcdf3

  ierr = nc.open(filename, PISM_NOWRITE); CHKERRQ(ierr);
  ierr = nc.inq_var(time_name, exists); CHKERRQ(ierr);
  if (exists == false) {
    ierr = nc.close(); CHKERRQ(ierr);

    PetscPrintf(com, "PISM ERROR: '%s' variable is not present in '%s'.\n",
                time_name.c_str(), filename.c_str());
    PISMEnd();
  }
  ierr = nc.get_att_text(time_name, "units", time_units); CHKERRQ(ierr);
  ierr = nc.get_att_text(time_name, "bounds", time_bounds_name); CHKERRQ(ierr);

  if (time_bounds_name.empty() == false) {
    ierr = nc.inq_var(time_bounds_name, exists); CHKERRQ(ierr);

    if (exists == false) {
      ierr = nc.close(); CHKERRQ(ierr);

      PetscPrintf(com, "PISM ERROR: '%s' variable is not present in '%s'.\n",
                  time_bounds_name.c_str(), filename.c_str());
      PISMEnd();
    }
  }

  // set the reference date:
  {
    size_t position = time_units.find("since");
    if (position == string::npos) {
      PetscPrintf(com, "PISM ERROR: time units string '%s' does not contain a reference date.\n",
                  time_units.c_str());
      PISMEnd();
    }

    reference_date = time_units.substr(position + 6); // 6 is the length of "since "
  }

  // set the time
  if (time_bounds_name.empty() == false) {
    // use the time bounds
    bounds.init(time_bounds_name, time_name, com, rank);
    bounds.set_units("seconds");

    // do *not* use the reference date
    ierr = bounds.read(nc, false, time); CHKERRQ(ierr);
  } else {
    // use the time axis
    time_axis.init(time_name, time_name, com, rank);
    time_axis.set_units("seconds");

    // do *not* use the reference date
    ierr = time_axis.read(nc, false, time); CHKERRQ(ierr);
  }

  run_start = time.front();
  run_end = time.back();
  time_in_seconds = run_start;

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

double PISMGregorianTime::mod(double time, double) {
  // This class does not support the "mod" operation.
  return time;
}

double PISMGregorianTime::year_fraction(double T) {
  int year, month, day, hour, minute;
  float second;
  double year_start, next_year_start;

  utCalendar(T, &ut_units,
             &year, &month, &day, &hour, &minute, &second);

  utInvCalendar(year,
                1, 1,            // month, day
                0, 0, 0,         // hour, minute, second
                &ut_units,
                &year_start);

  utInvCalendar(year + 1,
                1, 1,           // month, day
                0, 0, 0,        // hour, minute, second
                &ut_units,
                &next_year_start);

  return (T - year_start) / (next_year_start - year_start);
}

string PISMGregorianTime::date(double T) {
  char tmp[256];
  int year, month, day, hour, minute;
  float second;

  utCalendar(T, &ut_units,
             &year, &month, &day, &hour, &minute, &second);

  snprintf(tmp, 256, "%04d-%02d-%02d", year, month, day);

  return string(tmp);
}

string PISMGregorianTime::date() {
  return this->date(time_in_seconds);
}

string PISMGregorianTime::start_date() {
  return this->date(run_start);
}

string PISMGregorianTime::end_date() {
  return this->date(run_end);
}

