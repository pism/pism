// Copyright (C) 2011-2021 Constantine Khroulev
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

#include <cmath>
#include <cassert>              // assert
#include <cstdlib>              // strtol(), strtod()
#include <cstring>              // strlen()

#include "Time.hh"

#include "pism/external/calcalcs/calcalcs.h"

#include "ConfigInterface.hh"
#include "VariableMetadata.hh"
#include "pism_utilities.hh"
#include "error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/io_helpers.hh"
#include "pism/util/Logger.hh"

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

        return string_strip(time_units.substr(position + strlen("since")));

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
                  "         instead of the one present in the input file '%s' (%s)\n",
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
                  "         instead of the one present in the input file '%s' (%s)\n",
                  default_calendar.c_str(), input_file->filename().c_str(), calendar.c_str());
    }

    return calendar;
  }

  return default_calendar;
}

/*!
 * Increment the date corresponding to `T` by `years` years.
 */
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
    // avoid passing an impossible date to UDUNITS (no February 29 in non-leap years):
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

    units::DateTime d{numbers[0], numbers[1], numbers[2], 0, 0, 0.0};

    return time_units.time(d, calendar);
  } else {
    // "spec" must be a number or a number with units attached to it
    double t = 0.0;
    try {
      // check if strtod() can parse it:
      char *endptr = NULL;
      t = strtod(spec.c_str(), &endptr);
      if (*endptr == '\0') {
        // strtod() parsed it successfully: assume that it is in years. This will return
        // time in seconds.
        return increment_date(time_units, calendar, 0, t);
      } else {
        // strtod() failed -- assume that this is a number followed by units compatible
        // with seconds
        auto system = time_units.system();

        // Convert to seconds:
        return units::convert(system, 1.0, spec, "seconds");
      }
    } catch (RuntimeError &e) {
      e.add_context("parsing the date " + spec);
      throw;
    }
  }
}

/*!
 * Return the start time.
 */
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

//! Convert model years into seconds using the year length
//! corresponding to the current calendar.
/*! Do not use this to convert quantities other than time intervals!
 */
double Time::years_to_seconds(double input) const {
  return input * m_year_length;
}

//! Convert seconds into model years using the year length
//! corresponding to the current calendar.
/*! Do not use this to convert quantities other than time intervals!
 */
double Time::seconds_to_years(double input) const {
  return input / m_year_length;
}

void Time::init_calendar(const std::string &calendar_string) {

  if (not pism_is_valid_calendar_name(calendar_string)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "unsupported calendar: %s", calendar_string.c_str());
  }

  m_calendar_string = calendar_string;

  double seconds_per_day = convert(m_unit_system, 1.0, "day", "seconds");
  if (calendar_string == "360_day") {
    m_year_length = 360 * seconds_per_day;
  } else if (member(calendar_string, {"365_day", "noleap"})) {
    m_year_length = 365 * seconds_per_day;
  } else {
    // use the ~365.2524-day year
    m_year_length = convert(m_unit_system, 1.0, "year", "seconds");
  }
}

//! \brief Sets the current time (in seconds since the reference time).
void Time::set(double new_time) {
  m_time_in_seconds = new_time;
}

void Time::set_start(double new_start) {
  m_run_start = new_start;
}

void Time::set_end(double new_end) {
  m_run_end = new_end;
}

//! \brief Current time, in seconds.
double Time::current() const {
  return m_time_in_seconds;
}

double Time::start() const {
  return m_run_start;
}

double Time::end() const {
  return m_run_end;
}

std::string Time::units_string() const {
  return m_time_units.format();
}

//! \brief Returns the calendar string.
std::string Time::calendar() const {
  return m_calendar_string;
}

void Time::step(double delta_t) {
  m_time_in_seconds += delta_t;

  // If we are less than 0.001 second from the end of the run, reset
  // m_time_in_seconds to avoid taking a very small (and useless) time step.
  if (m_run_end > m_time_in_seconds &&
      m_run_end - m_time_in_seconds < 1e-3) {
    m_time_in_seconds = m_run_end;
  }
}

std::string Time::date() const {
  return date(current());
}

std::string Time::start_date() const {
  return date(m_run_start);
}

std::string Time::end_date() const {
  return date(m_run_end);
}

std::string Time::run_length() const {
  return pism::printf("%3.3f", seconds_to_years(m_run_end - m_run_start));
}

double Time::day_of_the_year_to_day_fraction(unsigned int day) const {
  const double sperd = 86400.0;
  return (sperd / m_year_length) * (double) day;
}

std::vector<double> Time::parse_times(const std::string &spec) const {
  if (spec.find(',') != std::string::npos) {
    // a list will always contain a comma because at least two numbers are
    // needed to specify reporting intervals
    return parse_list(spec);
  } else {
    // it must be a range specification
    return parse_range(spec);
  }
}

std::vector<double> Time::parse_list(const std::string &spec) const {
  std::vector<double> result;

  try {
    for (const auto &s : split(spec, ',')) {
      result.emplace_back(parse_date(s, m_time_units, m_calendar_string));
    }
  } catch (RuntimeError &e) {
      e.add_context("parsing a list of dates %s", spec.c_str());
      throw;
  }

  return result;
}

/**
 * Parses an interval specification string.
 *
 * @param[in] spec specification string
 *
 */
auto Time::parse_interval_length(const std::string &spec) const -> Interval {

  // check if it is a keyword
  if (spec == "hourly") {
    return {3600.0, SIMPLE};
  }

  if (spec == "daily") {
    return {86400.0, SIMPLE};
  }

  if (spec == "monthly") {
    return {1.0, MONTHLY};
  }

  if (spec == "yearly") {
    return {1.0, YEARLY};
  }

  if (not m_simple_calendar) {
    if (spec.find("year") != std::string::npos or spec.find("month") != std::string::npos) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "interval length '%s' with the calendar '%s' is not supported",
                                    spec.c_str(), m_calendar_string.c_str());
    }
  }

  try {
    units::Unit seconds(m_unit_system, "seconds"),
      one(m_unit_system, "1"),
      tmp(m_unit_system, spec);

    // Check if these units are compatible with "seconds" or "1". The
    // latter allows intervals of the form "0.5", which stands for "half
    // of a model year". This also discards interval specs such as "days
    // since 1-1-1", even though "days" is compatible with "seconds".
    if (units::are_convertible(tmp, seconds)) {
      units::Converter c(tmp, seconds);

      return {c(1.0), SIMPLE};
    } else if (units::are_convertible(tmp, one)) {
      units::Converter c(tmp, one);

      // convert from years to seconds without using UDUNITS-2 (this
      // way we handle 360-day and 365-day years correctly)
      return {years_to_seconds(c(1.0)), SIMPLE};
    }
  } catch (RuntimeError &e) {
    e.add_context("processing interval length " + spec);
    throw;
  }

  throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                "interval length '%s' is invalid", spec.c_str());
}


std::vector<double> Time::parse_range(const std::string &spec) const {
  double
    time_start   = m_run_start,
    time_end     = m_run_end;

  Interval I{0.0, SIMPLE};

  if (spec == "hourly") {
    I = {3600.0, SIMPLE};
  } else if (spec == "daily") {
    I = {86400.0, SIMPLE};
  } else if (spec == "monthly") {
    I = {1.0, MONTHLY};
  } else if (spec == "yearly") {
    I = {1.0, YEARLY};
  } else {

    auto parts = pism::split(spec, ':');

    if (parts.size() == 1) {
      I = parse_interval_length(parts[0]);

    } else if (parts.size() == 3) {
      time_start = parse_date(parts[0], m_time_units, m_calendar_string);
      I          = parse_interval_length(parts[1]);
      time_end   = parse_date(parts[2], m_time_units, m_calendar_string);
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "a time range must consist of exactly 3 parts separated by colons (got '%s').",
                                    spec.c_str());
    }
  }

  std::vector<double> result;
  compute_times(time_start, time_end, I, result);
  return result;
}

/**
 * Compute times corresponding to a "simple" time range.
 *
 * This is a time range with a fixed step ("hourly" is an example).
 *
 * @param[in] time_start beginning of the time interval, in seconds
 * @param[in] delta step, in seconds
 * @param[in] time_end end of the interval, in seconds
 * @param[out] result list of model times
 *
 */
void Time::compute_times_simple(double time_start, double delta, double time_end,
                                std::vector<double> &result) const {
  if (time_start >= time_end) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "a >= b in time range a:dt:b (got a = %s, b = %s)",
                                  this->date(time_start).c_str(), this->date(time_end).c_str());
  }

  if (delta <= 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "dt <= 0 in time range a:dt:b (got dt = %f seconds)", delta);
  }

  int k = 0;
  double t = time_start;

  result.clear();
  do {
    if (t >= this->start() && t <= this->end()) {
      result.push_back(t);
    }
    k += 1;
    t = time_start + k * delta;
  } while (t <= time_end);
}

double Time::convert_time_interval(double T, const std::string &units) const {
  if (member(units, {"year", "years", "yr", "a"})) {
    return this->seconds_to_years(T); // uses year length here
  }
  return convert(m_unit_system, T, "seconds", units);
}

/*!

  See http://meteora.ucsd.edu/~pierce/calcalcs/index.html and
  http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html#calendar

  for more details about supported calendars.
 */
Time::Time(MPI_Comm com, Config::ConstPtr config,
           const Logger &log,
           units::System::Ptr unit_system)
  : m_config(config),
    m_unit_system(unit_system),
    m_time_units(unit_system, "seconds since 1-1-1") {

  auto input_file = config->get_string("input.file");

  std::unique_ptr<File> file{};
  if (not input_file.empty()) {
    file.reset(new File(com, input_file, PISM_NETCDF3, PISM_READONLY));
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

  auto time_file = config->get_string("time.file");
  bool continue_run = config->get_flag("time.file.continue");

  if (not time_file.empty()) {
    log.message(2,
                "* Setting time from '%s'...\n",
                time_file.c_str());
    init_from_file(com, time_file, log, not continue_run);
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
void Time::init_from_file(MPI_Comm com,
                          const std::string &filename,
                          const Logger &log,
                          bool set_start_time) {
  try {
    std::string time_name = m_config->get_string("time.dimension_name");

    File file(com, filename, PISM_NETCDF3, PISM_READONLY); // OK to use netcdf3

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

      io::read_time_bounds(file, bounds, log, time);
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

double Time::year_fraction(double T) const {

  auto D = m_time_units.date(T, m_calendar_string);

  units::DateTime D2{D.year, 1, 1, 0, 0, 0.0};

  auto year_start = m_time_units.time(D2, m_calendar_string);

  D2.year += 1;
  auto next_year_start = m_time_units.time(D2, m_calendar_string);

  return (T - year_start) / (next_year_start - year_start);
}

std::string Time::date(double T) const {
  auto date = m_time_units.date(T, m_calendar_string);

  double hour = date.hour + date.minute / 60.0 + date.second / 3600.0;

  return pism::printf("%04d-%02d-%02d %04.1fh",
                      date.year, date.month, date.day, hour);
}

double Time::calendar_year_start(double T) const {
  auto D = m_time_units.date(T, m_calendar_string);

  units::DateTime D2{D.year, 1, 1, 0, 0, 0.0};

  return m_time_units.time(D2, m_calendar_string);
}


double Time::increment_date(double T, double years) const {
  return ::pism::increment_date(m_time_units, m_calendar_string, T, years);
}

void Time::compute_times_monthly(std::vector<double> &result) const {

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

void Time::compute_times_yearly(std::vector<double> &result) const {

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

void Time::compute_times(double time_start, double time_end,
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
