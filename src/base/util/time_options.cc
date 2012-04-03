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

#include "pism_options.hh"
#include "NCVariable.hh"
#include <algorithm>
#include <sstream>

//! \brief Compute
/*!
 * Here a and b are integer years.
 */
vector<double> compute_times(MPI_Comm com, const NCConfigVariable &config,
                             int a, int b, string keyword) {
  utUnit unit;
  string unit_str = "seconds since " + config.get_string("reference_date");
  vector<double> result;
  double a_offset, b_offset;

  // scan the units:
  int err = utScan(unit_str.c_str(), &unit);
  if (err != 0) {
    PetscPrintf(com, "PISM ERROR: invalid units specification: %s\n",
                unit_str.c_str());
    PISMEnd();
  }

  // get the 'year' out of the reference date:
  int reference_year, tmp1;
  float tmp2;
  utCalendar(0, &unit, &reference_year,
             &tmp1, &tmp1, &tmp1, &tmp1, &tmp2);

  // compute the number of seconds-since-the-reference date 'a' corresponds to:
  utInvCalendar(reference_year + a, // year
                1, 1,               // month, day
                0, 0, 0,            // hour, minute, second
                &unit,
                &a_offset);

  // compute the number of seconds-since-the-reference date 'b' corresponds to:
  utInvCalendar(reference_year + b, // year
                1, 1,               // month, day
                0, 0, 0,            // hour, minute, second
                &unit,
                &b_offset);

  if (keyword == "daily") {

    double t = a_offset, delta = 60*60*24; // seconds per day
    int year;

    do {
      result.push_back(t);
      t += delta;
      utCalendar(t, &unit, &year, &tmp1, &tmp1, &tmp1, &tmp1, &tmp2);
    } while (year <= reference_year + b);

    // add the last record:
    result.push_back(t);

  } else if (keyword == "monthly") {

    double t;
    int y, m;
    for (y = a; y <= b; y++) {
      for (m = 1; m <= 12; m++) {
        utInvCalendar(reference_year + y,   // year
                      m, 1,                 // month, day
                      0, 0, 0,              // hour, minute, second
                      &unit,
                      &t);
        result.push_back(t);
      }
    }

    // add the last record:
    utInvCalendar(reference_year + b + 1,   // year
                  1, 1,                     // month, day
                  0, 0, 0,                  // hour, minute, second
                  &unit,
                  &t);
    result.push_back(t);

  } else if (keyword == "yearly") {

    double t;
    for (int y = a; y <= b+1; y++) {    // note the "b + 1"
      utInvCalendar(reference_year + y,   // year
                    1, 1,                 // month, day
                    0, 0, 0,              // hour, minute, second
                    &unit,
                    &t);
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
PetscErrorCode parse_times(MPI_Comm com, const NCConfigVariable &config, string str, vector<double> &result) {
  PetscErrorCode ierr;
  int N;

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

      result = compute_times(com, config, (int)a, (int)b, keyword);

    } else {
      if (delta <= 0) {
        ierr = PetscPrintf(com, "PISM ERROR: delta <= 0 in the range specification %s.\n",
                           str.c_str()); CHKERRQ(ierr);
        return 1;
      }

      N = (int)floor((b - a)/delta) + 1; // number of points in the range
      result.resize(N);

      for (int j = 0; j < N; ++j)
        result[j] = convert(a + delta*j, "years", "seconds");
    }

  } else {			// it's a list of times
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
	result.push_back(convert(d, "years", "seconds"));
    }
    sort(result.begin(), result.end());
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
    if (j == 1 && (numbers[j] == "daily" ||
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
PetscErrorCode interval_to_seconds(MPI_Comm com, string interval, double &result) {
  utUnit interval_units, seconds;
  double slope, intercept;
  int errcode;

  errcode = utScan("seconds", &seconds);
  if (errcode != 0)
    SETERRQ(PETSC_COMM_SELF, 1, "utScan(\"seconds\", ...) failed");

  errcode = utScan(interval.c_str(), &interval_units);
  if (errcode != 0) {
    PetscPrintf(com, "PISM ERROR: can't parse '%s'\n", interval.c_str());
    PISMEnd();
  }

  if (utIsTime(&interval_units) == 1) {
    errcode = utConvert(&interval_units, &seconds, &slope, &intercept);
    if (errcode != 0) {
      PetscPrintf(com, "PISM ERROR: can't convert '%s' to seconds\n", interval.c_str());
      PISMEnd();
    }

    result = slope;
  } else {
    errcode = utScan("1", &seconds);
    if(errcode != 0)
      SETERRQ(PETSC_COMM_SELF, 1, "utScan(\"1\", ...) failed");

    errcode = utConvert(&interval_units, &seconds, &slope, &intercept);
    if (errcode != 0) {
      PetscPrintf(com, "PISM ERROR: can't parse '%s'\n", interval.c_str());
      PISMEnd();
    }

    result = slope;
  }

  return 0;
}
