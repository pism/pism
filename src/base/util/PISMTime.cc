// Copyright (C) 2011, 2012, 2013, 2014 Constantine Khroulev
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

#include "PISMConfig.hh"
#include "PISMTime.hh"
#include "pism_options.hh"
#include <sstream>
#include <assert.h>

namespace pism {

//! Convert model years into seconds using the year length
//! corresponding to the current calendar.
/*! Do not use this to convert quantities other than time intervals!
 */
double Time::years_to_seconds(double input) {
  return input * m_year_length;
}

//! Convert seconds into model years using the year length
//! corresponding to the current calendar.
/*! Do not use this to convert quantities other than time intervals!
 */
double Time::seconds_to_years(double input) {
  return input / m_year_length;
}


Time::Time(MPI_Comm c,
           const Config &conf,
           const std::string &calendar_string,
           const UnitSystem &unit_system)
  : m_com(c), m_config(conf), m_unit_system(unit_system),
    m_time_units(m_unit_system) {

  init_calendar(calendar_string);

  m_run_start = years_to_seconds(m_config.get("start_year"));
  m_run_end   = increment_date(m_run_start, (int)m_config.get("run_length_years"));

  m_time_units.parse("seconds");

  m_time_in_seconds = m_run_start;
}

Time::~Time() {
}

void Time::init_calendar(const std::string &calendar_string) {
  m_calendar_string = calendar_string;

  double seconds_per_day = m_unit_system.convert(1.0, "day", "seconds");
  if (calendar_string == "360_day") {
    m_year_length = 360 * seconds_per_day;
  } else if (calendar_string == "365_day" || calendar_string == "noleap") {
    m_year_length = 365 * seconds_per_day;
  } else {
    // use the ~365.2524-day year
    m_year_length = m_unit_system.convert(1.0, "year", "seconds");
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
double Time::current() {
  return m_time_in_seconds;
}

double Time::start() {
  return m_run_start;
}

double Time::end() {
  return m_run_end;
}

std::string Time::CF_units_string() {
  return "seconds since " + m_config.get_string("reference_date");
}

//! \brief Returns the calendar string.
std::string Time::calendar() {
  return m_calendar_string;
}

void Time::step(double delta_t) {
  m_time_in_seconds += delta_t;

  // If we are less than 0.001 second from the end of the run, reset
  // m_time_in_seconds to avoid taking a very small (and useless) time step.
  if (m_run_end > m_time_in_seconds &&
      m_run_end - m_time_in_seconds < 1e-3)
    m_time_in_seconds = m_run_end;
}

std::string Time::units_string() {
  return "seconds";
}


std::string Time::CF_units_to_PISM_units(const std::string &input) {
  std::string units = input;
  size_t n = units.find("since");

  /*!
    \note This code finds the string "since" in the units_string and
    terminates it on the first 's' of "since", if this sub-string was found.
    This is done to ignore the reference date in the time units string (the
    reference date specification always starts with this word).
  */
  if (n != std::string::npos)
    units.resize(n);

  // strip trailing spaces
  while (ends_with(units, " ") && units.empty() == false)
    units.resize(units.size() - 1); // this would fail on empty strings

  return units;
}

PetscErrorCode Time::process_ys(double &result, bool &flag) {
  PetscErrorCode ierr;
  result = m_config.get("start_year");

  ierr = OptionsReal("-ys", "Start year", result, flag); CHKERRQ(ierr);

  result = years_to_seconds(result);

  return 0;
}

PetscErrorCode Time::process_y(double &result, bool &flag) {
  PetscErrorCode ierr;
  result = m_config.get("run_length_years");

  ierr = OptionsReal("-y", "Run length, in years", result, flag); CHKERRQ(ierr);

  result = years_to_seconds(result);

  return 0;
}

PetscErrorCode Time::process_ye(double &result, bool &flag) {
  PetscErrorCode ierr;
  result = m_config.get("start_year") + m_config.get("run_length_years");

  ierr = OptionsReal("-ye", "End year", result, flag); CHKERRQ(ierr);

  result = years_to_seconds(result);

  return 0;
}

PetscErrorCode Time::init() {
  PetscErrorCode ierr;
  bool y_set, ys_set, ye_set;
  double y_seconds, ys_seconds, ye_seconds;

  // At this point the calendar and the year length are set (in the
  // constructor). The Time_Calendar class will (potentially)
  // override all this by using settings from -time_file, so that is
  // fine, too.

  ierr = PetscOptionsBegin(m_com, "", "PISM model time options", ""); CHKERRQ(ierr);
  {
    ierr = process_y(y_seconds, y_set); CHKERRQ(ierr);
    ierr = process_ys(ys_seconds, ys_set); CHKERRQ(ierr);
    ierr = process_ye(ye_seconds, ye_set); CHKERRQ(ierr);
  }
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  if (ys_set && ye_set && y_set) {
    ierr = PetscPrintf(m_com, "PISM ERROR: all of -y, -ys, -ye are set. Exiting...\n");
    CHKERRQ(ierr);
    PISMEnd();
  }

  if (y_set && ye_set) {
    ierr = PetscPrintf(m_com, "PISM ERROR: using -y and -ye together is not allowed. Exiting...\n"); CHKERRQ(ierr);
    PISMEnd();
  }

  // Set the start year if -ys is set, use the default otherwise.
  if (ys_set == true) {
    m_run_start = ys_seconds;
  }

  m_time_in_seconds = m_run_start;

  if (ye_set == true) {
    if (ye_seconds < m_time_in_seconds) {
      ierr = PetscPrintf(m_com,
                        "PISM ERROR: -ye (%s) is less than -ys (%s) (or input file year or default).\n"
                        "PISM cannot run backward in time.\n",
                         date(ye_seconds).c_str(),
                         date(m_run_start).c_str()); CHKERRQ(ierr);
      PISMEnd();
    }
    m_run_end = ye_seconds;
  } else if (y_set == true) {
    m_run_end = m_run_start + y_seconds;
  } else {
    m_run_end = increment_date(m_run_start, (int)m_config.get("run_length_years"));
  }

  return 0;
}

std::string Time::date(double T) {
  char tmp[256];
  snprintf(tmp, 256, "%012.3f", seconds_to_years(T));
  return std::string(tmp);
}

std::string Time::date() {
  return date(m_time_in_seconds);
}

std::string Time::start_date() {
  return date(m_run_start);
}

std::string Time::end_date() {
  return date(m_run_end);
}

std::string Time::run_length() {
  char tmp[256];
  snprintf(tmp, 256, "%3.3f", seconds_to_years(m_run_end - m_run_start));
  return std::string(tmp);
}

double Time::mod(double time, unsigned int period_years) {
  if (period_years == 0)
    return time;

  double period_seconds = years_to_seconds(period_years);

  double tmp = time - floor(time / period_seconds) * period_seconds;

  if (fabs(tmp - period_seconds) < 1)
    tmp = 0;

  return tmp;
}

double Time::year_fraction(double T) {
  double Y = seconds_to_years(T);
  return Y - floor(Y);
}

double Time::day_of_the_year_to_day_fraction(unsigned int day) {
  const double sperd = 86400.0;
  return (sperd / m_year_length) * (double) day;
}

double Time::calendar_year_start(double T) {
  return T - this->mod(T, 1);
}

double Time::increment_date(double T, int years) {
  return T + years_to_seconds(years);
}

PetscErrorCode Time::parse_times(const std::string &spec,
                                     std::vector<double> &result) {

  if (spec.find(',') != std::string::npos) {
    // a list will always contain a comma because at least two numbers are
    // needed to specify reporting intervals

    return parse_list(spec, result);

  } else {
    // it must be a range specification
    return parse_range(spec, result);
  }

  return 0;
}

PetscErrorCode Time::parse_list(const std::string &spec, std::vector<double> &result) {
  std::vector<std::string> parts;
  std::istringstream arg(spec);
  std::string tmp;

  while(getline(arg, tmp, ','))
    parts.push_back(tmp);

  result.clear();
  for (unsigned int k = 0; k < parts.size(); ++k) {
    double d;
    int errcode = parse_date(parts[k], &d);
    if (errcode != 0)
      return 1;

    result.push_back(d);
  }

  return 0;
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
int Time::parse_interval_length(const std::string &spec, std::string &keyword, double *result) {

  // check if it is a keyword
  if (spec == "hourly") {
    keyword = "simple";
    if (result)
      *result = 3600;
    return 0;
  }

  if (spec == "daily") {
    keyword = "simple";
    if (result)
      *result = 86400;
    return 0;
  }

  if (spec == "monthly" || spec == "yearly") {
    keyword = spec;
    if (result)
      *result = 0;
    return 0;
  }

  Unit tmp(m_time_units.get_system()),
    seconds(m_time_units.get_system()),
    one(m_time_units.get_system());

  int errcode;
  errcode = seconds.parse("seconds");
  assert(errcode == 0);
  errcode = one.parse("1");
  assert(errcode == 0);

  // check if the interval spec is a valid unit spec:
  if (tmp.parse(spec) != 0) {
    PetscPrintf(m_com, "PISM ERROR: invalid interval length: '%s'\n",
                spec.c_str());
    return 1;
  }

  // Check if these units are compatible with "seconds" or "1". The
  // latter allows intervals of the form "0.5", which stands for "half
  // of a model year". This also discards interval specs such as "days
  // since 1-1-1", even though "days" is compatible with "seconds".
  if (units_are_convertible(tmp, seconds) == true) {
    cv_converter *c = seconds.get_converter_from(tmp);
    assert(c != NULL);

    if (result)
      *result = cv_convert_double(c, 1.0);

    cv_free(c);

  } else if (units_are_convertible(tmp, one) == true) {
    cv_converter *c = one.get_converter_from(tmp);
    assert(c != NULL);

    if (result) {
      // this is a rather convoluted way of turning a string into a
      // floating point number:
      *result = cv_convert_double(c, 1.0);
      // convert from years to seconds without using UDUNITS-2 (this
      // way we handle 360-day and 365-day years correctly)
      *result = years_to_seconds(*result);
    }

    cv_free(c);
  } else {
    PetscPrintf(m_com, "PISM ERROR: invalid interval length: '%s'\n",
                spec.c_str());
    return 1;
  }

  return 0;
}


PetscErrorCode Time::parse_range(const std::string &spec, std::vector<double> &result) {
  double
    time_start   = m_run_start,
    time_end     = m_run_end,
    delta        = 0;
  std::string keyword = "simple";
  int    errcode;

  if (spec == "hourly") {
    delta = 3600;
  } else if (spec == "daily") {
    delta = 86400;
  } else if (spec == "monthly" || spec == "yearly") {
    keyword = spec;
    delta   = 0;
  } else {
    std::istringstream  arg(spec);
    std::vector<std::string> parts;
    std::string         tmp;

    while(getline(arg, tmp, ':'))
      parts.push_back(tmp);

    if (parts.size() == 1) {
      errcode = parse_interval_length(parts[0], keyword, &delta);
      if (errcode != 0)
        return 1;

    } else if (parts.size() == 3) {
      errcode = parse_date(parts[0], &time_start);
      if (errcode != 0)
        return 1;

      errcode = parse_interval_length(parts[1], keyword, &delta);
      if (errcode != 0)
        return 1;

      errcode = parse_date(parts[2], &time_end);
      if (errcode != 0)
        return 1;

    } else {
      PetscPrintf(m_com,
                  "PISM ERROR: A time range must consist of exactly 3 parts, separated by colons. (Got '%s'.)\n",
                  spec.c_str());
      return 1;
    }
  }

  return compute_times(time_start, delta, time_end, keyword, result);
}


PetscErrorCode Time::parse_date(const std::string &spec, double *result) {
  PetscErrorCode ierr;
  double d;
  char *endptr;

  if (spec.empty() == true) {
    ierr = PetscPrintf(m_com,
                       "PISM ERROR: got an empty string '%s'.\n",
                       spec.c_str()); CHKERRQ(ierr);
    return 1;
  }

  d = strtod(spec.c_str(), &endptr);
  if (*endptr != '\0') {
    ierr = PetscPrintf(m_com, "PISM ERROR: '%s' is not a number.\n",
                       spec.c_str()); CHKERRQ(ierr);
    return 1;
  }

  if (result)
    *result = years_to_seconds(d);

  return 0;
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
 * @return 0 on success, 1 otherwise
 */
PetscErrorCode Time::compute_times_simple(double time_start, double delta, double time_end,
                                              std::vector<double> &result) {
  if (time_start >= time_end) {
    PetscPrintf(m_com, "PISM ERROR: a >= b in time range a:dt:b.\n");
    return 1;
  }

  if (delta <= 0) {
    PetscPrintf(m_com, "PISM ERROR: dt <= 0 in time range a:dt:b.\n");
    return 1;
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
  } while (t < time_end);

  return 0;
}

PetscErrorCode Time::compute_times(double time_start, double delta, double time_end,
                                       const std::string &keyword,
                                       std::vector<double> &result) {
  if (keyword == "yearly") {
    delta = years_to_seconds(1.0);
  } else if (keyword == "monthly") {
    delta = years_to_seconds(1.0/12.0);
  } else if (keyword != "simple") {
    PetscPrintf(m_com, "PISM ERROR: Unknown time range keyword: %s.\n",
                keyword.c_str());
    return 1;
  }

  return compute_times_simple(time_start, delta, time_end, result);
}

double Time::convert_time_interval(double T, const std::string &units) {
  if (units == "year" || units == "years" || units == "yr" || units == "a") {
    return this->seconds_to_years(T); // uses year length here
  }
  return m_unit_system.convert(T, "seconds", units);
}

} // end of namespace pism
