// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017 PISM Authors
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

#include <cassert>
#include <sstream>
#include <cstdlib>
#include <petscsys.h>

#include "error_handling.hh"

#include "Time_Calendar.hh"
#include "pism_options.hh"
#include "pism/util/io/PIO.hh"
#include "pism/external/calcalcs/utCalendar2_cal.h"
#include "pism/external/calcalcs/calcalcs.h"
#include "ConfigInterface.hh"
#include "VariableMetadata.hh"
#include "io/io_helpers.hh"
#include "pism/util/Logger.hh"

namespace pism {

static inline std::string string_strip(std::string input) {
  if (input.empty() == true) {
    return input;
  }

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
Time_Calendar::Time_Calendar(MPI_Comm c, Config::ConstPtr conf,
                             const std::string &calendar_string,
                             units::System::Ptr units_system)
  : Time(conf, calendar_string, units_system),
    m_com(c) {

  std::string ref_date = m_config->get_string("time.reference_date");

  try {
    parse_date(ref_date, NULL);
  } catch (RuntimeError &e) {
    e.add_context("validating the reference date");
    throw;
  }

  try {
    m_time_units = units::Unit(m_time_units.get_system(), "seconds since " + ref_date);
  } catch (RuntimeError &e) {
    e.add_context("setting time units");
    throw;
  }

  m_run_start = increment_date(0, (int)m_config->get_double("time.start_year"));
  m_run_end   = increment_date(m_run_start, (int)m_config->get_double("time.run_length"));

  m_time_in_seconds = m_run_start;
}

Time_Calendar::~Time_Calendar() {
}

bool Time_Calendar::process_ys(double &result) {

  options::String ys("-ys", "Start date");

  if (ys.is_set()) {
    try {
      parse_date(ys, &result);
    } catch (RuntimeError &e) {
      e.add_context("processing the -ys option");
      throw;
    }
  } else {
    result = m_config->get_double("time.start_year", "seconds");
  }
  return ys.is_set();
}

bool Time_Calendar::process_y(double &result) {

  options::Integer y("-y", "Run length, in years (integer)", 0);

  if (y.is_set()) {
    if (y < 0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "-y %d is not allowed (run length can't be negative)",
                                    y.value());
    }
    result = years_to_seconds(y);
  } else {
    result = m_config->get_double("time.run_length", "seconds");
  }
  return y.is_set();
}


bool Time_Calendar::process_ye(double &result) {

  options::String ye("-ye", "Start date");

  if (ye.is_set()) {
    try {
      parse_date(ye, &result);
    } catch (RuntimeError &e) {
      e.add_context("processing the -ye option");
      throw;
    }
  } else {
    result = (m_config->get_double("time.start_year", "seconds") +
              m_config->get_double("time.run_length", "seconds"));
  }
  return ye.is_set();
}

void Time_Calendar::init_from_input_file(const PIO &nc,
                                         const std::string &time_name,
                                         const Logger &log) {
  try {
    // Set the calendar name from file, unless we are re-starting from a PISM run using the "none"
    // calendar.
    std::string new_calendar = nc.get_att_text(time_name, "calendar");
    if (not new_calendar.empty() and
        not (new_calendar == "none")) {
      init_calendar(new_calendar);
    }

    // Set the reference date of internal units.
    {
      std::string date_string = reference_date_from_file(nc, time_name);
      m_time_units = units::Unit(m_unit_system, "seconds " + date_string);
    }

    // Read time information from the file. (PISM output files don't have time bounds, so we don't
    // bother checking for them.)
    std::vector<double> time;
    {
      TimeseriesMetadata time_axis(time_name, time_name, m_unit_system);
      time_axis.set_string("units", m_time_units.format());

      io::read_timeseries(nc, time_axis, *this, log, time);
    }

    // Set time.
    this->set_start(time.front());
    this->set(time.front());
    log.message(2,
                "* Time t = %s (calendar: %s) found in '%s'; setting current time\n",
                this->date().c_str(),
                this->calendar().c_str(),
                nc.inq_filename().c_str());
  } catch (RuntimeError &e) {
    e.add_context("initializing model time from \"%s\"", nc.inq_filename().c_str());
    throw;
  }
}

void Time_Calendar::init(const Logger &log) {

  // process command-line options -y, -ys, -ye
  Time::init(log);

  options::String time_file("-time_file", "Reads time information from a file");

  if (time_file.is_set()) {
    log.message(2, 
                "* Setting time from '%s'...\n",
                time_file->c_str());

    options::ignored(log, "-y");
    options::ignored(log, "-ys");
    options::ignored(log, "-ye");

    bool continue_run = options::Bool("-time_file_continue_run",
                                      "continue a run using start time in the -i file");
    init_from_file(time_file, log, not continue_run);
  }
}

//! \brief Sets the time from a NetCDF file with a time dimension (`-time_file`).
/*!
 * This allows running PISM for the duration of the available forcing.
 */
void Time_Calendar::init_from_file(const std::string &filename, const Logger &log,
                                   bool set_start_time) {
  try {
    std::string time_name = m_config->get_string("time.dimension_name");

    PIO nc(m_com, "netcdf3", filename, PISM_READONLY); // OK to use netcdf3

    // Set the calendar name from file.
    std::string new_calendar = nc.get_att_text(time_name, "calendar");
    if (not new_calendar.empty()) {
      init_calendar(new_calendar);
    }

    // Set the reference date of internal units.
    {
      std::string date_string = reference_date_from_file(nc, time_name);
      m_time_units = units::Unit(m_unit_system, "seconds " + date_string);
    }

    // Read time information from the file.
    std::vector<double> time;
    std::string time_bounds_name = nc.get_att_text(time_name, "bounds");
    if (not time_bounds_name.empty()) {
      // use the time bounds
      TimeBoundsMetadata bounds(time_bounds_name, time_name, m_unit_system);
      bounds.set_string("units", m_time_units.format());

      io::read_time_bounds(nc, bounds, *this, log, time);
    } else {
      // use the time axis
      TimeseriesMetadata time_axis(time_name, time_name, m_unit_system);
      time_axis.set_string("units", m_time_units.format());

      io::read_timeseries(nc, time_axis, *this, log, time);
    }

    // Set time.
    if (set_start_time) {
      this->set_start(time.front());
      this->set(time.front());
    } else {
      log.message(2, "* Using start time from an -i file to continue an interrupted run.\n");
    }
    this->set_end(time.back());
  } catch (RuntimeError &e) {
    e.add_context("initializing model time from \"%s\"", filename.c_str());
    throw;
  }
}

double Time_Calendar::mod(double time, unsigned int) const {
  // This class does not support the "mod" operation.
  return time;
}

double Time_Calendar::year_fraction(double T) const {
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

std::string Time_Calendar::date(double T) const {
  char tmp[256];
  int year, month, day, hour, minute;
  double second;

  utCalendar2_cal(T, m_time_units.get(),
                  &year, &month, &day, &hour, &minute, &second,
                  m_calendar_string.c_str());

  snprintf(tmp, 256, "%04d-%02d-%02d", year, month, day);

  return std::string(tmp);
}

std::string Time_Calendar::date() const {
  return this->date(m_time_in_seconds);
}

std::string Time_Calendar::start_date() const {
  return this->date(m_run_start);
}

std::string Time_Calendar::end_date() const {
  return this->date(m_run_end);
}

double Time_Calendar::calendar_year_start(double T) const {
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
double Time_Calendar::increment_date(double T, int years) const {
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
    PetscErrorCode ierr = PetscPrintf(m_com,
                                      "PISM WARNING: date %d year(s) since %d-%d-%d does not exist."
                                      " Using %d-%d-%d instead of %d-%d-%d.\n",
                                      years,
                                      year, month, day,
                                      year + years, month, day-1,
                                      year + years, month, day);
    PISM_CHK(ierr, "PetscPrintf");
    day -= 1;
  }

  // Get the time in seconds corresponding to the new date.
  utInvCalendar2_cal(year + years, month, day,
                     hour, minute, second, m_time_units.get(), &result,
                     m_calendar_string.c_str());

  return result;
}

/**
 * Parses the date string `input`.
 *
 * @param[in] input date string of the form "Y-M-D", where year, month,
 *                 and day can use as many digits as necessary
 * @param[out] result
 *
 * If `result` is NULL, then it validates the date string, but does
 * not try to convert to seconds since the reference date. This makes
 * it possible to validate the reference date itself.
 *
 */
void Time_Calendar::parse_date(const std::string &input, double *result) const {
  std::vector<int> numbers;
  bool year_is_negative = false;
  std::string tmp, spec = input;

  spec = string_strip(spec);

  std::istringstream arg(spec);

  if (spec.empty() == true) {
    throw RuntimeError(PISM_ERROR_LOCATION, "got an empty date specification");
  }

  // If the string starts with "-" then the year is negative. This
  // would confuse the code below, which treats "-" as a separator, so
  // we remember that the year is negative and remove "-".
  if (spec[0] == '-') {
    year_is_negative = true;
    spec.substr(1);
  }

  while(getline(arg, tmp, '-')) {

    // an empty part in the date specification (corresponds to "--")
    if (tmp.empty() == true) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "date specification '%s' is invalid (can't have two '-' in a row)",
                                    spec.c_str());
    }

    // check if strtol can parse it:
    char *endptr = NULL;
    long int n = strtol(tmp.c_str(), &endptr, 10);
    if (*endptr != '\0') {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "date specification '%s' is invalid ('%s' is not an integer)",
                                    spec.c_str(), tmp.c_str());
    }

    // FIXME: this may overflow!
    numbers.push_back((int)n);
  }

  // wrong number of parts in the YYYY-MM-DD date:
  if (numbers.size() != 3) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "date specification '%s' is invalid (should have 3 parts: YYYY-MM-DD, got %d)",
                                  spec.c_str(), (int)numbers.size());
  }

  if (year_is_negative) {
    numbers[0] *= -1;
  }

  calcalcs_cal *cal = ccs_init_calendar(m_calendar_string.c_str());
  if (cal == NULL) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "calendar string '%s' is invalid", m_calendar_string.c_str());
  }

  int dummy = 0;
  int errcode = ccs_date2jday(cal, numbers[0], numbers[1], numbers[2], &dummy);
  ccs_free_calendar(cal);
  if (errcode != 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "date %s is invalid in the %s calendar",
                                  spec.c_str(), m_calendar_string.c_str());
  }

  // result is the *output* argument. If it is not NULL, then the user
  // asked us to convert a date to seconds since the reference time.
  if (result != NULL) {
    double time = 0.0;
    errcode = utInvCalendar2_cal(numbers[0], numbers[1], numbers[2], 0, 0, 0.0,
                                 m_time_units.get(), &time, m_calendar_string.c_str());
    if (errcode != 0) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "cannot convert '%s' to '%s'",
                                    spec.c_str(), m_time_units.format().c_str());
    }

    *result = time;
  }
}

void Time_Calendar::parse_interval_length(const std::string &spec,
                                          std::string &keyword, double *result) const {

  Time::parse_interval_length(spec, keyword, result);

  // This is called *only* if the 'spec' is *not* one of "monthly",
  // "yearly", "daily", "hourly", so we don't need to worry about
  // spec.find("...") finding "year" in "yearly".

  // do not allow intervals specified in terms of "fuzzy" units
  if (spec.find("year") != std::string::npos || spec.find("month") != std::string::npos) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "interval length '%s' with the calendar '%s' is not supported",
                                  spec.c_str(), m_calendar_string.c_str());
  }
}


void Time_Calendar::compute_times_monthly(std::vector<double> &result) const {
  int errcode;

  int year, month, day, hour, minute;
  double second;

  double time = m_run_start;
  // get the date corresponding to the current time
  errcode = utCalendar2_cal(time, m_time_units.get(),
                            &year, &month, &day,
                            &hour, &minute, &second,
                            m_calendar_string.c_str());
  PISM_C_CHK(errcode, 0, "utCalendar2_cal");

  result.clear();
  while (true) {
    // find the time corresponding to the beginning of the current
    // month
    errcode = utInvCalendar2_cal(year, month, 1, // year, month, day
                                 0, 0, 0.0,      // hour, minute, second
                                 m_time_units.get(), &time,
                                 m_calendar_string.c_str());
    PISM_C_CHK(errcode, 0, "utCalendar2_cal");

    if (time > m_run_end) {
      break;
    }

    if (time >= m_run_start && time <= m_run_end) {
      result.push_back(time);
    }

    if (month == 12) {
      year  += 1;
      month  = 1;
    } else {
      month += 1;
    }

  }
}

void Time_Calendar::compute_times_yearly(std::vector<double> &result) const {
  int errcode;

  int year, month, day, hour, minute;
  double second;

  double time = m_run_start;
  // get the date corresponding to the current time
  errcode = utCalendar2_cal(time, m_time_units.get(),
                            &year, &month, &day,
                            &hour, &minute, &second,
                            m_calendar_string.c_str());
  PISM_C_CHK(errcode, 0, "utCalendar2_cal");

  result.clear();
  while (true) {
    // find the time corresponding to the beginning of the current
    // year
    errcode = utInvCalendar2_cal(year, 1, 1, // year, month, day
                                 0, 0, 0.0,  // hour, minute, second
                                 m_time_units.get(), &time,
                                 m_calendar_string.c_str());
    PISM_C_CHK(errcode, 0, "utInvCalendar2_cal");

    if (time > m_run_end) {
      break;
    }

    if (time >= m_run_start && time <= m_run_end) {
      result.push_back(time);
    }

    year  += 1;
  }
}

void Time_Calendar::compute_times(double time_start, double delta, double time_end,
                                  const std::string &keyword,
                                  std::vector<double> &result) const {
  if (keyword == "simple") {
    compute_times_simple(time_start, delta, time_end, result);
  } else if (keyword == "monthly") {
    compute_times_monthly(result);
  } else if (keyword == "yearly") {
    compute_times_yearly(result);
  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "'%s' reporting is not implemented", keyword.c_str());
  }
}

} // end of namespace pism
