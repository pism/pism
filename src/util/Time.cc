// Copyright (C) 2011, 2012, 2013, 2014, 2015, 2016, 2017 Constantine Khroulev
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
#include <sstream>

#include "Time.hh"

#include "ConfigInterface.hh"
#include "Time_Calendar.hh"
#include "pism_options.hh"
#include "pism_utilities.hh"
#include "error_handling.hh"
#include "pism/util/io/PIO.hh"
#include "pism/util/Logger.hh"

namespace pism {

/**
 * Select a calendar using the "time.calendar" configuration parameter, the
 * "-calendar" command-line option, or the "calendar" attribute of the
 * "time" variable in the file specified using "-time_file".
 *
 */
std::string calendar_from_options(MPI_Comm com, const Config &config) {
  // Set the default calendar using the config. parameter or the
  // "-calendar" option:
  std::string result = config.get_string("time.calendar");

  // Check if -time_file was set and override the setting above if the
  // "calendar" attribute is found.
  options::String time_file("-time_file", "name of the file specifying the run duration");
  if (time_file.is_set()) {
    PIO nc(com, "netcdf3", time_file, PISM_READONLY);    // OK to use netcdf3

    std::string time_name = config.get_string("time.dimension_name");
    if (nc.inq_var(time_name)) {
      std::string tmp = nc.get_att_text(time_name, "calendar");
      if (not tmp.empty()) {
        result = tmp;
      }
    }
  }
  return result;
}

Time::Ptr time_from_options(MPI_Comm com, Config::ConstPtr config, units::System::Ptr system) {
  try {
    std::string calendar = calendar_from_options(com, *config);

    if (calendar == "360_day" or
        calendar == "365_day" or
        calendar == "noleap" or
        calendar == "none") {
      return Time::Ptr(new Time(config, calendar, system));
    } else {
      return Time::Ptr(new Time_Calendar(com, config, calendar, system));
    }

  } catch (RuntimeError &e) {
    e.add_context("initializing Time from options");
    throw;
  }
}

//! Initialize model time using command-line options and (possibly) files.
void initialize_time(MPI_Comm com, const std::string &dimension_name,
                     const Logger &log, Time &time) {

  // Check if we are initializing from a PISM output file:
  options::String input_file("-i", "Specifies a PISM input file");

  if (input_file.is_set()) {
    PIO nc(com, "netcdf3", input_file, PISM_READONLY);     // OK to use netcdf3
    time.init_from_input_file(nc, dimension_name, log);
  }

  time.init(log);
}

//! Get the reference date from a file.
std::string reference_date_from_file(const PIO &nc,
                                     const std::string &time_name) {

  if (not nc.inq_var(time_name)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "'%s' variable is not present in '%s'.",
                                  time_name.c_str(), nc.inq_filename().c_str());
  }
  std::string time_units = nc.get_att_text(time_name, "units");

  // Check if the time_units includes a reference date.
  size_t position = time_units.find("since");
  if (position == std::string::npos) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "time units string '%s' does not contain a reference date",
                                  time_units.c_str());
  }

  return time_units.substr(position);
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


Time::Time(Config::ConstPtr conf,
           const std::string &calendar_string,
           units::System::Ptr unit_system)
  : m_config(conf),
    m_unit_system(unit_system),
    m_time_units(m_unit_system, "seconds") {

  init_calendar(calendar_string);

  m_run_start = years_to_seconds(m_config->get_double("time.start_year"));
  m_run_end   = increment_date(m_run_start, (int)m_config->get_double("time.run_length"));

  m_time_in_seconds = m_run_start;
}

Time::~Time() {
}

void Time::init_calendar(const std::string &calendar_string) {

  if (not pism_is_valid_calendar_name(calendar_string)) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "unsupported calendar: %s", calendar_string.c_str());
  }

  m_calendar_string = calendar_string;

  double seconds_per_day = convert(m_unit_system, 1.0, "day", "seconds");
  if (calendar_string == "360_day") {
    m_year_length = 360 * seconds_per_day;
  } else if (calendar_string == "365_day" || calendar_string == "noleap") {
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

std::string Time::CF_units_string() const {
  return "seconds since " + m_config->get_string("time.reference_date");
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

std::string Time::units_string() const {
  return "seconds";
}


std::string Time::CF_units_to_PISM_units(const std::string &input) const {
  std::string units = input;
  size_t n = units.find("since");

  /*!
    \note This code finds the string "since" in the units_string and
    terminates it on the first 's' of "since", if this sub-string was found.
    This is done to ignore the reference date in the time units string (the
    reference date specification always starts with this word).
  */
  if (n != std::string::npos) {
    units.resize(n);
  }

  // strip trailing spaces
  while (ends_with(units, " ") && not units.empty()) {
    units.resize(units.size() - 1);  // this would fail on empty strings
  }

  return units;
}

bool Time::process_ys(double &result) {
  options::Real ys("-ys", "Start year", m_config->get_double("time.start_year"));
  result = years_to_seconds(ys);
  return ys.is_set();
}

bool Time::process_y(double &result) {
  options::Real y("-y", "Run length, in years", m_config->get_double("time.run_length"));
  result = years_to_seconds(y);
  return y.is_set();
}

bool Time::process_ye(double &result) {
  options::Real ye("-ye", "End year",
                      m_config->get_double("time.start_year") + m_config->get_double("time.run_length"));
  result = years_to_seconds(ye);
  return ye.is_set();
}


//! Set start time from a PISM input file.
/**
 * FIXME: This crude implementation does not use reference dates and does not convert units.
 */
void Time::init_from_input_file(const PIO &nc,
                                const std::string &time_name,
                                const Logger &log) {
  unsigned int time_length = nc.inq_dimlen(time_name);

  bool ys = options::Bool("-ys", "starting time");
  if (not ys and time_length > 0) {
    double T = 0.0;
    // Set the default starting time to be equal to the last time saved in the input file
    nc.inq_dim_limits(time_name, NULL, &T);
    this->set_start(T);
    this->set(T);
    log.message(2,
                "* Time t = %s found in '%s'; setting current time\n",
                this->date().c_str(), nc.inq_filename().c_str());
  }
}


void Time::init(const Logger &log) {

  (void) log;

  double y_seconds, ys_seconds, ye_seconds;

  // At this point the calendar and the year length are set (in the
  // constructor). The Time_Calendar class will (potentially)
  // override all this by using settings from -time_file, so that is
  // fine, too.

  bool y_set = process_y(y_seconds);
  bool ys_set = process_ys(ys_seconds);
  bool ye_set = process_ye(ye_seconds);

  if (ys_set && ye_set && y_set) {
    throw RuntimeError(PISM_ERROR_LOCATION, "all of -y, -ys, -ye are set.");
  }

  if (y_set && ye_set) {
    throw RuntimeError(PISM_ERROR_LOCATION, "using -y and -ye together is not allowed.");
  }

  // Set the start year if -ys is set, use the default otherwise.
  if (ys_set == true) {
    m_run_start = ys_seconds;
  }

  m_time_in_seconds = m_run_start;

  if (ye_set == true) {
    if (ye_seconds < m_time_in_seconds) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "-ye (%s) is less than -ys (%s) (or input file year or default).\n"
                                    "PISM cannot run backward in time.",
                                    date(ye_seconds).c_str(), date(m_run_start).c_str());
    }
    m_run_end = ye_seconds;
  } else if (y_set == true) {
    m_run_end = m_run_start + y_seconds;
  } else {
    m_run_end = increment_date(m_run_start, (int)m_config->get_double("time.run_length"));
  }
}

std::string Time::date(double T) const {
  char tmp[256];
  snprintf(tmp, 256, "%.3f", seconds_to_years(T));
  return std::string(tmp);
}

std::string Time::date() const {
  return date(m_time_in_seconds);
}

std::string Time::start_date() const {
  return date(m_run_start);
}

std::string Time::end_date() const {
  return date(m_run_end);
}

std::string Time::run_length() const {
  char tmp[256];
  snprintf(tmp, 256, "%3.3f", seconds_to_years(m_run_end - m_run_start));
  return std::string(tmp);
}

double Time::mod(double time, unsigned int period_years) const {
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

double Time::year_fraction(double T) const {
  double Y = seconds_to_years(T);
  return Y - floor(Y);
}

double Time::day_of_the_year_to_day_fraction(unsigned int day) const {
  const double sperd = 86400.0;
  return (sperd / m_year_length) * (double) day;
}

double Time::calendar_year_start(double T) const {
  return T - this->mod(T, 1);
}

double Time::increment_date(double T, int years) const {
  return T + years_to_seconds(years);
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
  std::istringstream arg(spec);
  std::string tmp;
  std::vector<double> result;

  while(getline(arg, tmp, ',')) {
    try {
      double d;
      parse_date(tmp, &d);
      result.push_back(d);
    } catch (RuntimeError &e) {
      e.add_context("parsing a list of dates %s", spec.c_str());
      throw;
    }
  }
  return result;
}

/**
 * Parses an interval specification string.
 *
 * @param[in] spec specification string
 * @param[out] keyword interval type keyword, one of "hourly",
 *                     "daily", "monthly", "yearly", "equal"
 * @param[out] result if `keyword` == "equal", `result` is set to
 *                    the interval length in seconds
 *
 * @return 0 on success, 1 otherwise
 */
void Time::parse_interval_length(const std::string &spec, std::string &keyword, double *result) const {

  // check if it is a keyword
  if (spec == "hourly") {
    keyword = "simple";
    if (result) {
      *result = 3600;
    }
    return;
  }

  if (spec == "daily") {
    keyword = "simple";
    if (result) {
      *result = 86400;
    }
    return;
  }

  if (spec == "monthly" || spec == "yearly") {
    keyword = spec;
    if (result) {
      *result = 0;
    }
    return;
  }

  units::Unit seconds(m_time_units.get_system(), "seconds"),
    one(m_time_units.get_system(), "1"),
    tmp = one;

  try {
    tmp = units::Unit(m_time_units.get_system(), spec);
  } catch (RuntimeError &e) {
    e.add_context("processing interval length " + spec);
    throw;
  }

  // Check if these units are compatible with "seconds" or "1". The
  // latter allows intervals of the form "0.5", which stands for "half
  // of a model year". This also discards interval specs such as "days
  // since 1-1-1", even though "days" is compatible with "seconds".
  if (units::are_convertible(tmp, seconds) == true) {
    units::Converter c(tmp, seconds);

    if (result) {
      *result = c(1.0);
    }

  } else if (units::are_convertible(tmp, one) == true) {
    units::Converter c(tmp, one);

    if (result) {
      // this is a rather convoluted way of turning a string into a
      // floating point number:
      *result = c(1.0);
      // convert from years to seconds without using UDUNITS-2 (this
      // way we handle 360-day and 365-day years correctly)
      *result = years_to_seconds(*result);
    }

  } else {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "interval length '%s' is invalid", spec.c_str());
  }
}


std::vector<double> Time::parse_range(const std::string &spec) const {
  double
    time_start   = m_run_start,
    time_end     = m_run_end,
    delta        = 0;
  std::string keyword = "simple";

  if (spec == "hourly") {
    delta = 3600;
  } else if (spec == "daily") {
    delta = 86400;
  } else if (spec == "monthly" || spec == "yearly") {
    keyword = spec;
    delta   = 0;
  } else {

    std::vector<std::string> parts = pism::split(spec, ':');

    if (parts.size() == 1) {
      parse_interval_length(parts[0], keyword, &delta);
    } else if (parts.size() == 3) {
      parse_date(parts[0], &time_start);
      parse_interval_length(parts[1], keyword, &delta);
      parse_date(parts[2], &time_end);
    } else {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "a time range must consist of exactly 3 parts separated by colons (got '%s').",
                                    spec.c_str());
    }
  }

  std::vector<double> result;
  compute_times(time_start, delta, time_end, keyword, result);
  return result;
}


void Time::parse_date(const std::string &spec, double *result) const {

  if (spec.empty() == true) {
    throw RuntimeError(PISM_ERROR_LOCATION, "got an empty date specification");
  }

  char *endptr = NULL;
  double d = strtod(spec.c_str(), &endptr);
  if (*endptr != '\0') {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "date specification '%s' is invalid ('%s' is not an number)",
                                  spec.c_str(), spec.c_str());
  }

  if (result) {
    *result = years_to_seconds(d);
  }
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

void Time::compute_times(double time_start, double delta, double time_end,
                         const std::string &keyword,
                         std::vector<double> &result) const {
  if (keyword == "yearly") {
    delta = years_to_seconds(1.0);
  } else if (keyword == "monthly") {
    delta = years_to_seconds(1.0/12.0);
  } else if (keyword != "simple") {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "unknown time range keyword: %s",
                                  keyword.c_str());
  }

  compute_times_simple(time_start, delta, time_end, result);
}

double Time::convert_time_interval(double T, const std::string &units) const {
  if (units == "year" || units == "years" || units == "yr" || units == "a") {
    return this->seconds_to_years(T); // uses year length here
  }
  return convert(m_unit_system, T, "seconds", units);
}

} // end of namespace pism
