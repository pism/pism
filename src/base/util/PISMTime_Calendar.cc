// Copyright (C) 2012, 2013, 2014 PISM Authors
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

#include <assert.h>
#include <sstream>
#include <stdlib.h>

#include "PISMTime_Calendar.hh"
#include "pism_options.hh"
#include "PIO.hh"
#include "utCalendar2_cal.h"
#include "calcalcs.h"
#include "PISMConfig.hh"

namespace pism {

static inline std::string string_strip(std::string input) {
  if (input.empty() == true)
    return input;

  // strip leading spaces
  input.erase(0, input.find_first_not_of(" \t"));

  // strip trailing spaces
  input.substr(input.find_last_not_of(" \t"));

  return input;
}

/*!

  See http://meteora.ucsd.edu/~pierce/calcalcs/index.html and
  http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.6/cf-conventions.html#calendar

  for more details about supported calendars.
 */
Time_Calendar::Time_Calendar(MPI_Comm c, const Config &conf,
                                     std::string calendar_string,
                                     UnitSystem units_system)
  : Time(c, conf, calendar_string, units_system) {

  // init_calendar() was called by the constructor of Time.
  if (pism_is_valid_calendar_name(m_calendar_string) == false) {
    PetscPrintf(m_com, "PISM ERROR: unsupported calendar: %s\n", m_calendar_string.c_str());
    PISMEnd();
  }

  std::string ref_date = m_config.get_string("reference_date");

  int errcode = parse_date(ref_date, NULL);
  if (errcode != 0) {
    PetscPrintf(m_com, "PISM ERROR: reference date %s is invalid.\n",
                ref_date.c_str());
    PISMEnd();
  }

  std::string tmp = "seconds since " + ref_date;
  errcode = m_time_units.parse(tmp);
  if (errcode != 0) {
    PetscPrintf(m_com, "PISM ERROR: time units '%s' are invalid.\n",
                tmp.c_str());
    PISMEnd();
  }

  m_run_start = increment_date(0, (int)m_config.get("start_year"));
  m_run_end   = increment_date(m_run_start, (int)m_config.get("run_length_years"));

  m_time_in_seconds = m_run_start;
}

Time_Calendar::~Time_Calendar() {
}

PetscErrorCode Time_Calendar::process_ys(double &result, bool &flag) {
  PetscErrorCode ierr;
  std::string tmp;
  result = m_config.get("start_year", "years", "seconds");

  ierr = OptionsString("-ys", "Start date", tmp, flag); CHKERRQ(ierr);

  if (flag) {
    ierr = parse_date(tmp, &result);
    if (ierr != 0) {
      PetscPrintf(m_com, "PISM ERROR: processing -ys option failed.\n");
      PISMEnd();
    }
  }

  return 0;
}

PetscErrorCode Time_Calendar::process_y(double &result, bool &flag) {
  PetscErrorCode ierr;
  int tmp;
  result = m_config.get("run_length_years", "years", "seconds");

  ierr = OptionsInt("-y", "Run length, in years (integer)", tmp, flag); CHKERRQ(ierr);

  if (flag) {
    result = years_to_seconds(tmp);
  }

  return 0;
}


PetscErrorCode Time_Calendar::process_ye(double &result, bool &flag) {
  PetscErrorCode ierr;
  std::string tmp;
  result = (m_config.get("start_year", "years", "seconds") +
            m_config.get("run_length_years", "years", "seconds"));

  ierr = OptionsString("-ye", "Start date", tmp, flag); CHKERRQ(ierr);

  if (flag) {
    ierr = parse_date(tmp, &result);
    if (ierr != 0) {
      PetscPrintf(m_com, "PISM ERROR: processing -ye option failed.\n");
      PISMEnd();
    }
  }

  return 0;
}


PetscErrorCode Time_Calendar::init() {
  PetscErrorCode ierr;
  std::string time_file;
  bool flag;

  ierr = Time::init(); CHKERRQ(ierr);

  ierr = OptionsString("-time_file", "Reads time information from a file",
                           time_file, flag); CHKERRQ(ierr);

  if (flag == true) {
    ierr = verbPrintf(2, m_com,
                      "* Setting time from '%s'...\n",
                      time_file.c_str()); CHKERRQ(ierr);

    ierr = ignore_option(m_com, "-y"); CHKERRQ(ierr);
    ierr = ignore_option(m_com, "-ys"); CHKERRQ(ierr);
    ierr = ignore_option(m_com, "-ye"); CHKERRQ(ierr);

    ierr = init_from_file(time_file); CHKERRQ(ierr);
  }

  return 0;
}

//! \brief Sets the time from a NetCDF with forcing data.
/*!
 * This allows running PISM for the duration of the available forcing.
 */
PetscErrorCode Time_Calendar::init_from_file(std::string filename) {
  PetscErrorCode ierr;
  int rank = 0;
  std::vector<double> time, time_bounds;
  std::string time_units, time_bounds_name, new_calendar,
    time_name = m_config.get_string("time_dimension_name");
  bool exists;

  NCTimeseries time_axis(time_name, time_name, m_unit_system);
  time_axis.set_units(m_time_units.format());

  ierr = MPI_Comm_rank(m_com, &rank); CHKERRQ(ierr);
  PIO nc(m_com, "netcdf3", m_unit_system); // OK to use netcdf3

  ierr = nc.open(filename, PISM_READONLY); CHKERRQ(ierr);
  ierr = nc.inq_var(time_name, exists); CHKERRQ(ierr);
  if (exists == false) {
    ierr = nc.close(); CHKERRQ(ierr);

    PetscPrintf(m_com, "PISM ERROR: '%s' variable is not present in '%s'.\n",
                time_name.c_str(), filename.c_str());
    PISMEnd();
  }
  ierr = nc.get_att_text(time_name, "units",  time_units);       CHKERRQ(ierr);
  ierr = nc.get_att_text(time_name, "bounds", time_bounds_name); CHKERRQ(ierr);

  ierr = nc.get_att_text(time_name, "calendar", new_calendar); CHKERRQ(ierr);
  if (new_calendar.empty() == false) {
    if (pism_is_valid_calendar_name(new_calendar) == false) {
      PetscPrintf(m_com,
                  "PISM ERROR: unsupported calendar name '%s' found in a -time_file '%s'.\n",
                  new_calendar.c_str(), filename.c_str());
      PISMEnd();
    }
    init_calendar(new_calendar);
  }

  if (time_bounds_name.empty() == false) {
    ierr = nc.inq_var(time_bounds_name, exists); CHKERRQ(ierr);

    if (exists == false) {
      ierr = nc.close(); CHKERRQ(ierr);

      PetscPrintf(m_com, "PISM ERROR: variable '%s' is not present in '%s'.\n",
                  time_bounds_name.c_str(), filename.c_str());
      PISMEnd();
    }
  }

  {
    // Check if the time_file has a reference date set:
    size_t position = time_units.find("since");
    if (position == std::string::npos) {
      PetscPrintf(m_com, "PISM ERROR: time units string '%s' does not contain a reference date.\n",
                  time_units.c_str());
      PISMEnd();
    }

    std::string tmp = "seconds " + time_units.substr(position);
    ierr = m_time_units.parse(tmp);
    if (ierr != 0) {
      PetscPrintf(m_com,
                  "PISM ERROR: units specification '%s' is invalid (processing -time_file).\n",
                  tmp.c_str());
      PISMEnd();

    }
  }

  // set the time
  if (time_bounds_name.empty() == false) {
    // use the time bounds
    NCTimeBounds bounds(time_bounds_name, time_name, m_unit_system);
    bounds.set_units(m_time_units.format());

    ierr = nc.read_time_bounds(bounds, this, time); CHKERRQ(ierr);
  } else {
    // use the time axis

    ierr = nc.read_timeseries(time_axis, this, time); CHKERRQ(ierr);
  }

  m_run_start       = time.front();
  m_run_end         = time.back();
  m_time_in_seconds = m_run_start;

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;
}

double Time_Calendar::mod(double time, unsigned int) {
  // This class does not support the "mod" operation.
  return time;
}

double Time_Calendar::year_fraction(double T) {
  int year, month, day, hour, minute;
  double second, year_start, next_year_start;

  utCalendar2_cal(T, m_time_units.get(),
                  &year, &month, &day, &hour, &minute, &second,
                  m_calendar_string.c_str());

  utInvCalendar2_cal(year,
                     1, 1,            // month, day
                     0, 0, 0,         // hour, minute, second
                     m_time_units.get(),
                     &year_start,
                     m_calendar_string.c_str());

  utInvCalendar2_cal(year + 1,
                     1, 1,           // month, day
                     0, 0, 0,        // hour, minute, second
                     m_time_units.get(),
                     &next_year_start,
                     m_calendar_string.c_str());

  return (T - year_start) / (next_year_start - year_start);
}

std::string Time_Calendar::date(double T) {
  char tmp[256];
  int year, month, day, hour, minute;
  double second;

  utCalendar2_cal(T, m_time_units.get(),
                  &year, &month, &day, &hour, &minute, &second,
                  m_calendar_string.c_str());

  snprintf(tmp, 256, "%04d-%02d-%02d", year, month, day);

  return std::string(tmp);
}

std::string Time_Calendar::date() {
  return this->date(m_time_in_seconds);
}

std::string Time_Calendar::start_date() {
  return this->date(m_run_start);
}

std::string Time_Calendar::end_date() {
  return this->date(m_run_end);
}

double Time_Calendar::calendar_year_start(double T) {
  int year, month, day, hour, minute;
  double second, result;

  // Get the date corresponding to time T:
  utCalendar2_cal(T, m_time_units.get(),
                  &year, &month, &day, &hour, &minute, &second,
                  m_calendar_string.c_str());

  // Get the time in seconds corresponding to the beginning of the
  // year.
  utInvCalendar2_cal(year,
                     1, 1, 0, 0, 0.0, // month, day, hour, minute, second
                     m_time_units.get(), &result,
                     m_calendar_string.c_str());

  return result;
}


// FIXME: this feeds invalid dates to utInvCalendar2_cal! (step 1 year from Feb 29...)
double Time_Calendar::increment_date(double T, int years) {
  int year, month, day, hour, minute;
  double second, result;

  // Get the date corresponding ti time T:
  utCalendar2_cal(T, m_time_units.get(),
                  &year, &month, &day, &hour, &minute, &second,
                  m_calendar_string.c_str());

  calcalcs_cal *cal = NULL;
  int errcode, leap = 0;
  cal = ccs_init_calendar(m_calendar_string.c_str());
  assert(cal != NULL);
  errcode = ccs_isleap(cal, year + years, &leap);
  assert(errcode == 0);
  ccs_free_calendar(cal);

  if (leap == 0 && month == 2 && day == 29) {
    PetscPrintf(m_com,
                "PISM WARNING: date %d year(s) since %d-%d-%d does not exist. Using %d-%d-%d instead of %d-%d-%d.\n",
                years,
                year, month, day,
                year + years, month, day-1,
                year + years, month, day);
    day -= 1;
  }

  // Get the time in seconds corresponding to the new date.
  utInvCalendar2_cal(year + years, month, day,
                     hour, minute, second, m_time_units.get(), &result,
                     m_calendar_string.c_str());

  return result;
}

/**
 * Parses the date string `spec`.
 *
 * @param[in] spec date string of the form "Y-M-D", where year, month,
 *                 and day can use as many digits as necessary
 * @param[out] result
 *
 * If `result` is NULL, then it validates the date string, but does
 * not try to convert to seconds since the reference date. This makes
 * it possible to validate the reference date itself.
 *
 * @return 0 on success, 1 otherwise
 */
PetscErrorCode Time_Calendar::parse_date(std::string spec, double *result) {
  int errcode, dummy;
  calcalcs_cal *cal = NULL;
  std::vector<int> numbers;
  bool year_is_negative = false;
  std::string tmp;

  spec = string_strip(spec);

  std::istringstream arg(spec);

  if (spec.empty() == true)
    goto failure;

  // ignore negative years
  if (spec[0] == '-') {
    year_is_negative = true;
    spec.substr(1);
  }

  while(getline(arg, tmp, '-')) {

    // an empty part in the date specification (corresponds to "--")
    if (tmp.empty() == true)
      goto failure;

    // check if strtol can parse it:
    char *endptr = NULL;
    long int n = strtol(tmp.c_str(), &endptr, 10);
    if (*endptr != '\0')
      goto failure;

    // FIXME: this may overflow!
    numbers.push_back((int)n);
  }

  // wrong number of parts in the YYYY-MM-DD date:
  if (numbers.size() != 3)
    goto failure;

  if (year_is_negative)
    numbers[0] *= -1;

  cal = ccs_init_calendar(m_calendar_string.c_str());
  assert(cal != NULL);
  errcode = ccs_date2jday(cal, numbers[0], numbers[1], numbers[2], &dummy);
  ccs_free_calendar(cal);

  if (result != NULL) {

    if (errcode != 0) {
      PetscPrintf(m_com, "PISM ERROR: date %s is invalid in the %s calendar\n",
                  spec.c_str(), m_calendar_string.c_str());
      goto failure;
    }

    double time;
    errcode = utInvCalendar2_cal(numbers[0], numbers[1], numbers[2], 0, 0, 0.0,
                                 m_time_units.get(), &time, m_calendar_string.c_str());
    if (errcode != 0)
      goto failure;

    *result = time;
  }


  return 0;

 failure:
  PetscPrintf(m_com, "PISM ERROR: invalid date: '%s'\n",
              spec.c_str());

  return 1;
}

int Time_Calendar::parse_interval_length(std::string spec, std::string &keyword, double *result) {

  int ierr;

  ierr = Time::parse_interval_length(spec, keyword, result);
  if (ierr != 0)
    return 1;

  // This is called *only* if the 'spec' is *not* one of "monthly",
  // "yearly", "daily", "hourly", so we don't need to worry about
  // spec.find("...") finding "year" in "yearly".

  // do not allow intervals specified in terms of "fuzzy" units
  if (spec.find("year") != std::string::npos || spec.find("month") != std::string::npos) {
    PetscPrintf(m_com, "PISM ERROR: interval length '%s' with the calendar '%s' is not supported.\n",
                spec.c_str(), m_calendar_string.c_str());
    return 1;
  }

  return 0;
}


PetscErrorCode Time_Calendar::compute_times_monthly(std::vector<double> &result) {
  int errcode;

  int year, month, day, hour, minute;
  double second;

  double time = m_run_start;
  // get the date corresponding to the current time
  errcode = utCalendar2_cal(time, m_time_units.get(),
                            &year, &month, &day,
                            &hour, &minute, &second,
                            m_calendar_string.c_str());
  assert(errcode == 0);

  result.clear();
  while (true) {
    // find the time corresponding to the beginning of the current
    // month
    errcode = utInvCalendar2_cal(year, month, 1, // year, month, day
                                 0, 0, 0.0,      // hour, minute, second
                                 m_time_units.get(), &time,
                                 m_calendar_string.c_str());

    if (time > m_run_end)
      break;

    if (time >= m_run_start && time <= m_run_end)
      result.push_back(time);

    if (month == 12) {
      year  += 1;
      month  = 1;
    } else {
      month += 1;
    }

  }

  return 0;
}

PetscErrorCode Time_Calendar::compute_times_yearly(std::vector<double> &result) {
  int errcode;

  int year, month, day, hour, minute;
  double second;

  double time = m_run_start;
  // get the date corresponding to the current time
  errcode = utCalendar2_cal(time, m_time_units.get(),
                            &year, &month, &day,
                            &hour, &minute, &second,
                            m_calendar_string.c_str());
  assert(errcode == 0);

  result.clear();
  while (true) {
    // find the time corresponding to the beginning of the current
    // year
    errcode = utInvCalendar2_cal(year, 1, 1, // year, month, day
                                 0, 0, 0.0,  // hour, minute, second
                                 m_time_units.get(), &time,
                                 m_calendar_string.c_str());

    if (time > m_run_end)
      break;

    if (time >= m_run_start && time <= m_run_end)
      result.push_back(time);

    year  += 1;
  }

  return 0;
}

PetscErrorCode Time_Calendar::compute_times(double time_start, double delta, double time_end,
                                                std::string keyword,
                                                std::vector<double> &result) {
  if (keyword == "simple") {
    return compute_times_simple(time_start, delta, time_end, result);
  }

  if (keyword == "monthly") {
    return compute_times_monthly(result);
  }

  if (keyword == "yearly") {
    return compute_times_yearly(result);
  }

  PetscErrorCode ierr = PetscPrintf(m_com,
                                    "PISM ERROR: '%s' reporting is not implemented.\n");
  CHKERRQ(ierr);

  return 1;
}

} // end of namespace pism
