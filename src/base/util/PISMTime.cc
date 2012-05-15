// Copyright (C) 2011, 2012 Constantine Khroulev
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

#include "PISMTime.hh"
#include "pism_options.hh"

PISMTime::PISMTime(MPI_Comm c, const NCConfigVariable &conf)
  : com(c), config(conf) {

  secpera = convert(1.0, "year", "seconds");

  reference_date = config.get_string("reference_date");
  calendar_string = "365_day";  // only 365_day is supported by this class

  run_start = config.get("start_year", "years", "seconds");
  run_end   = run_start + config.get("run_length_years", "years", "seconds");

  time_in_seconds = run_start;
}

PetscErrorCode PISMTime::init() {
  PetscErrorCode ierr;
  bool y_set, ys_set, ye_set;
  PetscReal y = config.get("run_length_years"),
    ys = config.get("start_year"),
    ye = ys + y;

  ierr = PetscOptionsBegin(com, "", "PISM model time options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsReal("-ys", "Start year", ys, ys_set); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-ye", "End year", ye, ye_set); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-y", "Run length, in years", y, y_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (ys_set && ye_set && y_set) {
    ierr = PetscPrintf(com, "PISM ERROR: all of -y, -ys, -ye are set. Exiting...\n");
    CHKERRQ(ierr);
    PISMEnd();
  }

  if (y_set && ye_set) {
    ierr = PetscPrintf(com, "PISM ERROR: using -y and -ye together is not allowed. Exiting...\n"); CHKERRQ(ierr);
    PISMEnd();
  }

  // Set the start year if -ys is set, use the default otherwise.
  if (ys_set == true) {
    run_start = ys * secpera;
  }

  time_in_seconds = run_start;

  if (ye_set == true) {
    if (ye < seconds_to_years(time_in_seconds)) {
      ierr = PetscPrintf(com,
			"PISM ERROR: -ye (%3.3f) is less than -ys (%3.3f) (or input file year or default).\n"
			"PISM cannot run backward in time.\n",
			 ye, seconds_to_years(run_start)); CHKERRQ(ierr);
      PISMEnd();
    }
    run_end = ye * secpera;
  } else if (y_set == true) {
    run_end = run_start + y * secpera;
  } else {
    run_end = run_start + config.get("run_length_years", "years", "seconds");
  }

  return 0;
}

string PISMTime::date(double T) {
  char tmp[256];
  snprintf(tmp, 256, "%3.3f", seconds_to_years(T));
  return string(tmp);
}

string PISMTime::date() {
  return date(time_in_seconds);
}

string PISMTime::start_date() {
  return date(run_start);
}

string PISMTime::end_date() {
  return date(run_end);
}

string PISMTime::run_length() {
  char tmp[256];
  snprintf(tmp, 256, "%3.3f", seconds_to_years(run_end - run_start));
  return string(tmp);
}

double PISMTime::mod(double time, double period) {
  if (period <= 0)
    return time;

  double tmp = time - floor(time / period) * period;

  if (fabs(tmp - period) < 1)
    tmp = 0;

  return tmp;
}

double PISMTime::year_fraction(double T) {
  return seconds_to_years(T) - floor(seconds_to_years(T));
}
