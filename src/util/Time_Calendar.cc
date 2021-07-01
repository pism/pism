// Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017, 2019, 2020, 2021 PISM Authors
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

#include <cassert>              // assert()
#include <cstdlib>              // strtol()
#include <petscsys.h>

#include "error_handling.hh"

#include "Time_Calendar.hh"
#include "pism_options.hh"
#include "pism/util/io/File.hh"
#include "pism/external/calcalcs/calcalcs.h"
#include "ConfigInterface.hh"
#include "VariableMetadata.hh"
#include "io/io_helpers.hh"
#include "pism/util/Logger.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {

//! Get the reference date from a file.
static std::string reference_date_from_file(const File &file,
                                            const std::string &time_name,
                                            const std::string &default_value,
                                            bool stop_on_error) {

  if (file.find_variable(time_name)) {
    std::string time_units = file.read_text_attribute(time_name, "units");

    if (not time_units.empty()) {
      // Check if the time_units includes a reference date.
      size_t position = time_units.find("since");

      if (position != std::string::npos) {

        std::string since = "since";
        return string_strip(time_units.substr(position + since.size()));

      } else if (stop_on_error) {

        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "%s:units = \"%s\" in '%s' does not contain a reference date",
                                      time_name.c_str(),
                                      time_units.c_str(),
                                      file.filename().c_str());
      }
    } else if (stop_on_error) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "the '%s' variable in '%s' has no units",
                                    time_name.c_str(), file.filename().c_str());
    }
  } else if (stop_on_error) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "'%s' variable is not present in '%s'.",
                                    time_name.c_str(), file.filename().c_str());
  }

  return default_value;
}

//! Get the reference date from a file.
static std::string calendar_from_file(const File &file,
                                      const std::string &time_name,
                                      const std::string &default_value,
                                      bool stop_on_error) {

  if (file.find_variable(time_name)) {
    std::string calendar_name = file.read_text_attribute(time_name, "calendar");

    if (not calendar_name.empty()) {
      return calendar_name;
    } else if (stop_on_error) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "the '%s' variable in '%s' has no calendar attribute",
                                    time_name.c_str(), file.filename().c_str());
    }

  } else if (stop_on_error) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "'%s' variable is not present in '%s'.",
                                  time_name.c_str(), file.filename().c_str());
  }

  return default_value;
}

static std::string reference_date(const File *input_file,
                                  const Config &config,
                                  const Logger &log) {

  auto default_reference_date = config.get_string("time.reference_date");

  if (input_file != nullptr) {
    // input file is not empty

    auto time = config.get_string("time.dimension_name");

    if (not config.get_flag("input.bootstrap")) {
      // restarting from a file: use the reference date in this file
      bool stop_on_error = true;
      return reference_date_from_file(*input_file, time, default_reference_date, stop_on_error);
    }

    // Bootstrapping: use the configuration parameter and warn about mismatches
    bool stop_on_error = false;
    auto ref_date = reference_date_from_file(*input_file, time, default_reference_date, stop_on_error);

    if (ref_date != default_reference_date) {
      log.message(2,
                  "WARNING: Using reference date %s\n"
                  "         instead of the one present in the input file '%s' (%s)",
                  default_reference_date.c_str(), input_file->filename().c_str(), ref_date.c_str());
    }

    return ref_date;
  }

  return default_reference_date;
}

static std::string calendar(const File *input_file,
                            const Config &config,
                            const Logger &log) {
  auto default_calendar = config.get_string("time.calendar");

  if (input_file != nullptr) {
    // input file is not empty

    auto time = config.get_string("time.dimension_name");

    if (not config.get_flag("input.bootstrap")) {
      // restarting from a file: use the calendar in this file
      bool stop_on_error = true;
      return calendar_from_file(*input_file, time, default_calendar, stop_on_error);
    }

    // Bootstrapping: use the configuration parameter and warn about mismatches
    bool stop_on_error = false;
    auto calendar = calendar_from_file(*input_file, time, default_calendar, stop_on_error);

    if (calendar != default_calendar) {
      log.message(2,
                  "WARNING: Using calendar %s\n"
                  "         instead of the one present in the input file '%s' (%s)",
                  default_calendar.c_str(), input_file->filename().c_str(), calendar.c_str());
    }

    return calendar;
  }

  return default_calendar;
}

static double increment_date(const units::Unit &time_units,
                             const std::string &calendar,
                             double T, double years) {
  assert(years >= 0.0);

  int whole_years = static_cast<int>(std::floor(years));
  double year_fraction = years - whole_years;
  const double day_length = 86400.0;

  // Get the date corresponding to time T:
  auto date = time_units.date(T, calendar);

  // shift the date by the number of whole years requested
  date.year += whole_years;

  // check if the resulting year is a leap year:
  int leap = 0;
  {
    calcalcs_cal *cal = ccs_init_calendar(calendar.c_str());
    assert(cal != NULL);
    int errcode = ccs_isleap(cal, date.year, &leap);
    assert(errcode == 0);
    ccs_free_calendar(cal);
  }

  double result = 0.0;
  if (leap == 0 and date.month == 2 and date.day == 29) {
    // avoid passing an impossible date to UDUNITS:
    date.day -= 1;
    result = time_units.time(date, calendar);
    // add back the day we substracted above
    result += day_length;
  } else {
    result = time_units.time(date, calendar);
  }

  int year_length = (leap == 1) ? 366 : 365;

  result += year_fraction * (year_length * day_length);

  return result;
}

/*!
 * Parse the date.
 *
 * `input` can be
 *
 * - a YYYY-MM-DD date (YYYY can be negative)
 *
 * - a number (interpreted as the number of years since the reference date in
 *   `time_units`)
 *
 * - a number with units attached ("1 day", etc) interpreted as time since the reference
 *   date in `time_units`
 */
static double parse_date(const std::string &input,
                         const units::Unit &time_units,
                         const std::string &calendar) {
  using units::Converter;
  using units::DateTime;
  using units::Unit;

  std::string spec = string_strip(input);

  if (spec.empty()) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "got an empty date specification");
  }

  // If the string starts with "-" then the year is negative. This
  // would confuse the code below, which treats "-" as a separator, so
  // we remember that the year is negative and remove "-".
  bool year_is_negative = false;
  if (spec[0] == '-') {
    year_is_negative = true;
    spec.substr(1);
  }

  auto parts = split(spec, '-');

  if (parts.size() == 3) {

    std::vector<int> numbers;
    for (const auto &p : parts) {
      // check if strtol can parse it:
      char *endptr = NULL;
      long int n = strtol(p.c_str(), &endptr, 10);
      if (*endptr != '\0') {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "date specification '%s' is invalid ('%s' is not an integer)",
                                      spec.c_str(), p.c_str());
      }

      // FIXME: this may overflow!
      numbers.push_back((int)n);
    }

    if (year_is_negative) {
      numbers[0] *= -1;
    }

    // Validate the calendar string and the date in this calendar:
    {
      calcalcs_cal *cal = ccs_init_calendar(calendar.c_str());
      if (cal == NULL) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "calendar string '%s' is invalid",
                                      calendar.c_str());
      }

      int dummy = 0;
      int errcode = ccs_date2jday(cal, numbers[0], numbers[1], numbers[2], &dummy);
      if (errcode != 0) {
        throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                      "date %s is invalid in the %s calendar",
                                      spec.c_str(), calendar.c_str());
      }
      ccs_free_calendar(cal);
    }

    DateTime d{numbers[0], numbers[1], numbers[2], 0, 0, 0.0};

    return time_units.time(d, calendar);
  } else {
    // "spec" must be a number or a number with units attached to it
    double t = 0.0;
    try {
      // check if strtod() can parse it:
      char *endptr = NULL;
      t = strtod(spec.c_str(), &endptr);
      if (*endptr != '\0') {
        // strtod() failed -- assume that this is a number followed by units compatible
        // with seconds
        auto system = time_units.system();

        t = units::convert(system, 1.0, spec, "seconds");
      }

      return increment_date(time_units, calendar, 0, t);
    } catch (RuntimeError &e) {
      e.add_context("parsing the date " + spec);
      throw;
    }
  }
}

static double start_time(const Config &config,
                         const Logger &log,
                         const File *file,
                         const std::string &reference_date,
                         const std::string &calendar,
                         const units::Unit &time_units) {

  auto time_start = config.get_string("time.start");

  if (not time_start.empty()) {
    return parse_date(time_start, time_units, calendar);
  }

  if (file == nullptr) {
    // 0.0 corresponds to the reference date
    return 0.0;
  }

  // get the calendar in this file
  auto time_name     = config.get_string("time.dimension_name");
  bool stop_on_error = false;
  auto file_calendar = calendar_from_file(*file, time_name, calendar, stop_on_error);
  auto ref_date      = reference_date_from_file(*file, time_name, reference_date, stop_on_error);

  if (file_calendar != calendar) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "calendar in '%s' (%s) does not match the selected calendar (%s)",
                                  file->filename().c_str(), file_calendar.c_str(), calendar.c_str());
  }

  if (ref_date != reference_date) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "reference date in '%s' (%s) does not match the selected date (%s)",
                                  file->filename().c_str(), ref_date.c_str(), reference_date.c_str());
  }

  // FIXME: it would make sense to get the length of the time dimension and read the last
  // number instead.
  if (file->dimension_length(time_name) > 0) {
    VariableMetadata time_axis(time_name, time_units.system());
    time_axis.set_string("units", time_units.format());

    std::vector<double> time{};
    io::read_timeseries(*file, time_axis, log, time);

    return time.back();
  }

  return 0.0;
}

static double end_time(const Config &config,
                       double time_start,
                       const std::string &calendar,
                       const units::Unit &time_units) {
  auto time_end = config.get_string("time.end");

  if (not time_end.empty()) {
    // parse use time_end and use it
    return parse_date(time_end, time_units, calendar);
  } else {
    auto run_length = config.get_number("time.run_length", "seconds");
    // use time_start and run_length
    return time_start + run_length;
  }
}

/*!

  See http://meteora.ucsd.edu/~pierce/calcalcs/index.html and
  http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html#calendar

  for more details about supported calendars.
 */
Time_Calendar::Time_Calendar(MPI_Comm c, Config::ConstPtr config,
                             const Logger &log,
                             units::System::Ptr unit_system)
  : Time(config, unit_system), // FIXME calendar
    m_com(c) {

  // FIXME: implement -time_file

  auto input_file = config->get_string("input.file");

  std::unique_ptr<File> file{};
  if (not input_file.empty()) {
    file.reset(new File(m_com, input_file, PISM_NETCDF3, PISM_READONLY));
  }

  // set the reference date
  auto ref_date = reference_date(file.get(), *config, log);

  try {
    // this will validate the reference date
    m_time_units = units::Unit(m_unit_system, "seconds since " + ref_date);
  } catch (RuntimeError &e) {
    e.add_context("setting time units");
    throw;
  }

  m_calendar_string = ::pism::calendar(file.get(), *config, log);
  init_calendar(m_calendar_string);
  m_simple_calendar = member(m_calendar_string, {"360_day", "365_day", "no_leap"});

  m_run_start = start_time(*config,
                           log,
                           file.get(),
                           ref_date,
                           m_calendar_string,
                           m_time_units);

  m_run_end = end_time(*config,
                       m_run_start,
                       m_calendar_string,
                       m_time_units);

  m_time_in_seconds = m_run_start;
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
    result = m_config->get_number("time.run_length", "seconds");
  }
  return y.is_set();
}


bool Time_Calendar::process_ye(double &result) {

  options::String ye("-ye", "Start date");

  if (ye.is_set()) {
    try {
      result = parse_date(ye);
    } catch (RuntimeError &e) {
      e.add_context("processing the -ye option");
      throw;
    }
  } else {
    result = (m_config->get_number("time.start_year", "seconds") +
              m_config->get_number("time.run_length", "seconds"));
  }
  return ye.is_set();
}

/*!
 * Set the calendar, reference date, the start time, and the current time from a file.
 */
void Time_Calendar::init_from_input_file(const File &,
                                         const std::string &,
                                         const Logger &) {
  // FIXME: remove
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
 * Sets
 * - calendar
 * - reference date
 * - start time
 * - current time
 * - end time
 *
 * This allows running PISM for the duration of the available forcing.
 */
void Time_Calendar::init_from_file(const std::string &filename, const Logger &log,
                                   bool set_start_time) {
  try {
    std::string time_name = m_config->get_string("time.dimension_name");

    File file(m_com, filename, PISM_NETCDF3, PISM_READONLY); // OK to use netcdf3

    // Set the calendar name from file.
    std::string new_calendar = file.read_text_attribute(time_name, "calendar");
    if (not new_calendar.empty()) {
      init_calendar(new_calendar);
    }

    // Set the reference date of internal units.
    {
      bool stop_on_error = true;
      std::string date_string = reference_date_from_file(file,
                                                         time_name,
                                                         "FIXME",
                                                         stop_on_error);
      m_time_units = units::Unit(m_unit_system, "seconds " + date_string);
    }

    // Read time information from the file.
    std::vector<double> time;
    std::string time_bounds_name = file.read_text_attribute(time_name, "bounds");
    if (not time_bounds_name.empty()) {
      // use the time bounds
      VariableMetadata bounds(time_bounds_name, m_unit_system);
      bounds.set_string("units", m_time_units.format());

      io::read_time_bounds(file, bounds, *this, log, time);
    } else {
      // use the time axis
      VariableMetadata time_axis(time_name, m_unit_system);
      time_axis.set_string("units", m_time_units.format());

      io::read_timeseries(file, time_axis, log, time);
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

double Time_Calendar::mod(double time, unsigned int period_years) const {

  if (m_simple_calendar) {
    if (period_years == 0) {
      return time;
    }

    double period_seconds = years_to_seconds(period_years);

    double tmp = time - floor(time / period_seconds) * period_seconds;

    if (fabs(tmp - period_seconds) < 1) {
      tmp = 0;
    }

    return tmp;
  }

  // Other calendars don't have a consistent year length and do not support this
  // operation.
  return time;
}

double Time_Calendar::year_fraction(double T) const {

  auto D = m_time_units.date(T, m_calendar_string);

  units::DateTime D2{D.year, 1, 1, 0, 0, 0.0};

  auto year_start = m_time_units.time(D2, m_calendar_string);

  D2.year += 1;
  auto next_year_start = m_time_units.time(D2, m_calendar_string);

  return (T - year_start) / (next_year_start - year_start);
}

std::string Time_Calendar::date(double T) const {
  auto date = m_time_units.date(T, m_calendar_string);

  double hour = date.hour + date.minute / 60.0 + date.second / 3600.0;

  return pism::printf("%04d-%02d-%02d %04.1fh",
                      date.year, date.month, date.day, hour);
}

double Time_Calendar::calendar_year_start(double T) const {
  auto D = m_time_units.date(T, m_calendar_string);

  units::DateTime D2{D.year, 1, 1, 0, 0, 0.0};

  return m_time_units.time(D2, m_calendar_string);
}


double Time_Calendar::increment_date(double T, double years) const {
  return ::pism::increment_date(m_time_units, m_calendar_string, T, years);
}

auto Time_Calendar::parse_interval_length(const std::string &spec) const -> Interval {
  // do not allow intervals specified in terms of "fuzzy" units

  // This is called *only* if the 'spec' is *not* one of "monthly",
  // "yearly", "daily", "hourly", so we don't need to worry about
  // spec.find("...") finding "year" in "yearly".
  if (spec.find("year") != std::string::npos or spec.find("month") != std::string::npos) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "interval length '%s' with the calendar '%s' is not supported",
                                  spec.c_str(), m_calendar_string.c_str());
  }

  return Time::parse_interval_length(spec);
}

void Time_Calendar::compute_times_monthly(std::vector<double> &result) const {

  double time = m_run_start;

  // get the date corresponding to the current time
  auto date = m_time_units.date(time, m_calendar_string);

  // beginning of the current month
  units::DateTime d{date.year, date.month, 1, 0, 0, 0.0};

  result.clear();
  while (true) {
    // find the time corresponding to the beginning of the current month
    time = m_time_units.time(d, m_calendar_string);

    if (time > m_run_end) {
      break;
    }

    if (time >= m_run_start and time <= m_run_end) {
      result.push_back(time);
    }

    if (d.month == 12) {
      d.year  += 1;
      d.month  = 1;
    } else {
      d.month += 1;
    }
  }
}

void Time_Calendar::compute_times_yearly(std::vector<double> &result) const {

  double time = m_run_start;
  // get the date corresponding to the current time
  auto date = m_time_units.date(time, m_calendar_string);

  units::DateTime d{date.year, 1, 1, 0, 0, 0.0};

  result.clear();
  while (true) {
    // find the time corresponding to the beginning of the current year
    time = m_time_units.time(d, m_calendar_string);

    if (time > m_run_end) {
      break;
    }

    if (time >= m_run_start and time <= m_run_end) {
      result.push_back(time);
    }

    d.year += 1;
  }
}

void Time_Calendar::compute_times(double time_start, double time_end,
                                  const Interval &interval,
                                  std::vector<double> &result) const {
  switch (interval.type) {
  case SIMPLE:
    compute_times_simple(time_start, interval.dt, time_end, result);
    break;
  case MONTHLY:
    compute_times_monthly(result);
    break;
  case YEARLY:
    compute_times_yearly(result);
    break;
  }
}

} // end of namespace pism
