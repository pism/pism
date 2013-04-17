// Copyright (C) 2012, 2013 PISM Authors
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

#include "pism_options.hh"
#include "NCVariable.hh"
#include <algorithm>
#include <sstream>
#include "PISMUnits.hh"
#include "utCalendar2_cal.h"

//! \brief Compute
/*!
 * Here a and b are integer years.
 */
vector<double> compute_times(MPI_Comm com, const NCConfigVariable &config,
                             PISMUnitSystem unit_system,
                             int a, int b, string keyword) {
  PISMUnit time_units;
  string unit_str = "seconds since " + config.get_string("reference_date"),
    calendar = config.get_string("calendar");
  vector<double> result;
  double a_offset, b_offset;

  // scan the units:
  if (time_units.parse(unit_system, unit_str) != 0) {
    PetscPrintf(com, "PISM ERROR: invalid units specification: %s\n",
                unit_str.c_str());
    PISMEnd();
  }

  // get the 'year' out of the reference date:
  int reference_year, tmp1;
  double tmp2;
  utCalendar2_cal(0, time_units.get(), &reference_year,
                  &tmp1, &tmp1, &tmp1, &tmp1, &tmp2,
                  calendar.c_str());

  // compute the number of seconds-since-the-reference date 'a' corresponds to:
  utInvCalendar2_cal(reference_year + a, // year
                     1, 1,               // month, day
                     0, 0, 0,            // hour, minute, second
                     time_units.get(),
                     &a_offset,
                     calendar.c_str());

  // compute the number of seconds-since-the-reference date 'b' corresponds to:
  utInvCalendar2_cal(reference_year + b, // year
                     1, 1,               // month, day
                     0, 0, 0,            // hour, minute, second
                     time_units.get(),
                     &b_offset,
                     calendar.c_str());

  if (keyword == "hourly" || keyword == "daily") {

    double t = a_offset, delta = 60*60; // seconds per hour
    int year;

    if (keyword == "daily")
      delta *= 24;              // seconds per day

    do {
      result.push_back(t);
      t += delta;
      utCalendar2_cal(t, time_units.get(), &year, &tmp1, &tmp1, &tmp1, &tmp1, &tmp2,
                      calendar.c_str());
    } while (year <= reference_year + b);

    // add the last record:
    result.push_back(t);

  } else if (keyword == "monthly") {

    double t;
    int y, m;
    for (y = a; y <= b; y++) {
      for (m = 1; m <= 12; m++) {
        utInvCalendar2_cal(reference_year + y,   // year
                           m, 1,                 // month, day
                           0, 0, 0,              // hour, minute, second
                           time_units.get(),
                           &t,
                           calendar.c_str());
        result.push_back(t);
      }
    }

    // add the last record:
    utInvCalendar2_cal(reference_year + b + 1,   // year
                       1, 1,                     // month, day
                       0, 0, 0,                  // hour, minute, second
                       time_units.get(),
                       &t,
                       calendar.c_str());
    result.push_back(t);

  } else if (keyword == "yearly") {

    double t;
    for (int y = a; y <= b+1; y++) {    // note the "b + 1"
      utInvCalendar2_cal(reference_year + y,   // year
                         1, 1,                 // month, day
                         0, 0, 0,              // hour, minute, second
                         time_units.get(),
                         &t,
                         calendar.c_str());
      result.push_back(t);
    }

  } else {
    PetscPrintf(com,
                "PISM ERROR: unknown time-step keyword: %s\n"
                "            (only 'daily', 'monthly' and 'yearly' are implemented).\n",
                keyword.c_str());
    PISMEnd();
  }

  return result;
}

//! Parses a time specification.
/*!
  If it is a MATLAB-style range, then calls parse_range and computes all the points.

  If it is a comma-separated list, converts to double (with error-checking).
 */
PetscErrorCode parse_times(MPI_Comm com, const NCConfigVariable &config,
                           PISMUnitSystem unit_system, string str,
                           double run_start, double run_end, vector<double> &result) {
  PetscErrorCode ierr;
  int N;
  string calendar = config.get_string("calendar");

  if (str.find(':') != string::npos) { // it's a range specification

    double a, delta, b;
    string keyword;
    ierr = parse_range(com, str, &a, &delta, &b, keyword);
    if (ierr != 0) return 1;

    if (a >= b) {
      ierr = PetscPrintf(com, "PISM ERROR: a >= b in the range specification %s.\n",
			 str.c_str()); CHKERRQ(ierr);
      return 1;
    }

    if (keyword != "simple") {

      result = compute_times(com, config, unit_system,
                             (int)a, (int)b, keyword);

    } else {
      if (delta <= 0) {
        ierr = PetscPrintf(com, "PISM ERROR: delta <= 0 in the range specification %s.\n",
                           str.c_str()); CHKERRQ(ierr);
        return 1;
      }

      N = (int)floor((b - a)/delta) + 1; // number of points in the range
      result.resize(N);

      for (int j = 0; j < N; ++j)
        result[j] = convert(a + delta*j, unit_system, "years", "seconds");
    }

  } else if (str == "hourly"  ||
             str == "daily"   ||
             str == "monthly" ||
             str == "yearly") { // it is a keyword without the range

    PISMUnit time_units;
    string unit_str = "seconds since " + config.get_string("reference_date");

    // scan the units:
    if (time_units.parse(unit_system, unit_str) != 0) {
      PetscPrintf(com, "PISM ERROR: invalid units specification: %s\n",
                  unit_str.c_str());
      PISMEnd();
    }

    // get the 'year' out of the reference date:
    int reference_year, start_year, end_year, tmp1;
    double tmp2;
    utCalendar2_cal(0, time_units.get(), &reference_year,
                    &tmp1, &tmp1, &tmp1, &tmp1, &tmp2,
                    calendar.c_str());

    // get the year at the start of the run
    utCalendar2_cal(run_start, time_units.get(), &start_year,
                    &tmp1, &tmp1, &tmp1, &tmp1, &tmp2,
                    calendar.c_str());

    // get the year at the end of the run
    utCalendar2_cal(run_end, time_units.get(), &end_year,
                    &tmp1, &tmp1, &tmp1, &tmp1, &tmp2,
                    calendar.c_str());

    result = compute_times(com, config, unit_system,
                           start_year - reference_year,
                           end_year - reference_year,
                           str);
    
  } else if (str.find(',') != string::npos) {			// it's a list of times
    // a list will always contain a comma because at least two numbers are
    // needed to specify reporting intervals
    istringstream arg(str);
    string tmp;

    result.clear();
    while(getline(arg, tmp, ',')) {
      double d;
      char *endptr;

      d = strtod(tmp.c_str(), &endptr);
      if (*endptr != '\0') {
	ierr = PetscPrintf(com, "PISM ERROR: Can't parse %s (%s is not a number).\n",
			   str.c_str(), tmp.c_str()); CHKERRQ(ierr);
	return 1;
      }
      else
	result.push_back(convert(d, unit_system, "years", "seconds"));
    }
    sort(result.begin(), result.end());
  } else {
    PetscPrintf(com, "PISM ERROR: Can't parse %s", str.c_str());
    PISMEnd();
  }

  return 0;
}


//! Parses a MATLAB-style range (a:delta:b).
PetscErrorCode parse_range(MPI_Comm com, string str, double *a, double *delta, double *b, string &keyword) {
  PetscErrorCode ierr;
  istringstream arg(str);
  vector<string> numbers;
  double doubles[3];

  // Split the string:
  string tmp;
  while (getline(arg, tmp, ':'))
    numbers.push_back(tmp);

  // Check if we have 3 numbers:
  if (numbers.size() != 3) {
      ierr = PetscPrintf(com,
			 "PISM ERROR: A range has to consist of exactly three numbers, separated by colons.\n");
      CHKERRQ(ierr);
      return 1;
  }

  keyword = "simple";

  // Convert each number from a string to double:
  for (int j = 0; j < 3; ++j) {
    double d;
    char *endptr;

    // take care of daily, monthly and yearly keywords:
    if (j == 1 && (numbers[j] == "hourly"  ||
                   numbers[j] == "daily"   ||
                   numbers[j] == "monthly" ||
                   numbers[j] == "yearly")) {
      keyword = numbers[j];
      doubles[j] = 0;
      continue;
    }

    d = strtod(numbers[j].c_str(), &endptr);
    if (*endptr != '\0') {
      ierr = PetscPrintf(com, "PISM ERROR: Can't parse %s (%s is not a number).\n",
			 str.c_str(), numbers[j].c_str()); CHKERRQ(ierr);
      return 1;
    }
    else
      doubles[j] = d;
  }

  if (a) *a = doubles[0];
  if (delta) *delta = doubles[1];
  if (b) *b = doubles[2];

  return 0;
}

/*!
 * Converts time interval specifications like "10days", "5hours", "day",
 * "month" to seconds.
 */
PetscErrorCode interval_to_seconds(MPI_Comm com,
                                   PISMUnitSystem unit_system,
                                   string interval, double &result) {
  PISMUnit interval_units, seconds;
  cv_converter *c;

  if (seconds.parse(unit_system, "seconds") != 0)
    SETERRQ(PETSC_COMM_SELF, 1, "parsing 'seconds' failed");

  if (interval_units.parse(unit_system, interval)) {
    PetscPrintf(com, "PISM ERROR: can't parse '%s'\n", interval.c_str());
    PISMEnd();
  }

  if (ut_are_convertible(interval_units.get(), seconds.get()) != 0) {
    c = ut_get_converter(interval_units.get(), seconds.get());
    if (c == NULL) {
      PetscPrintf(com, "PISM ERROR: can't convert '%s' to seconds\n", interval.c_str());
      PISMEnd();
    }

    result = cv_convert_double(c, 1.0);
    cv_free(c);
  } else {
    PISMUnit one;
    if (one.parse(unit_system, "1") != 0)
      SETERRQ(PETSC_COMM_SELF, 1, "ut_parse(..., \"1\", ...) failed");

    c = ut_get_converter(interval_units.get(), one.get());
    if (c == NULL) {
      PetscPrintf(com, "PISM ERROR: can't parse '%s'\n", interval.c_str());
      PISMEnd();
    }

    result = cv_convert_double(c, 1.0);
  }
  cv_free(c);
  
  return 0;
}
