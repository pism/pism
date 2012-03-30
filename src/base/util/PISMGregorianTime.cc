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

PISMGregorianTime::PISMGregorianTime(MPI_Comm c, const NCConfigVariable &conf)
  : PISMTime(c, conf) {

  calendar_string = "gregorian";  // only "gregorian" is supported by this class
}

PetscErrorCode PISMGregorianTime::init() {
  PetscErrorCode ierr;

  ierr = PISMTime::init(); CHKERRQ(ierr);

  // initialize the units object:
  ierr = utScan(this->units().c_str(), &ut_units); CHKERRQ(ierr);

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

