// Copyright (C) 2011, 2012, 2013 Constantine Khroulev
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
#include <sstream>
#include <assert.h>

PISMTime::PISMTime(MPI_Comm c, const NCConfigVariable &conf, PISMUnitSystem unit_system)
  : m_com(c), m_config(conf), m_unit_system(unit_system) {

  m_secpera = convert(1.0, m_unit_system, "year", "seconds");

  m_reference_date = m_config.get_string("reference_date");
  m_calendar_string = "none";  // No calendar

  m_run_start = m_config.get("start_year", m_unit_system, "years", "seconds");
  m_run_end   = m_run_start + m_config.get("run_length_years", m_unit_system, "years", "seconds");

  m_time_in_seconds = m_run_start;
}

PISMTime::~PISMTime() {
}

//! \brief Sets the current time (in seconds since the reference time).
void PISMTime::set(double new_time) {
  m_time_in_seconds = new_time;
}

void PISMTime::set_start(double new_start) {
  m_run_start = new_start;
}

void PISMTime::set_end(double new_end) {
  m_run_end = new_end;
}

//! \brief Current time, in seconds.
double PISMTime::current() {
  return m_time_in_seconds;
}

double PISMTime::start() {
  return m_run_start;
}

double PISMTime::end() {
  return m_run_end;
}

double PISMTime::seconds_to_years(double T) {
  return T / secpera;
}

double PISMTime::years_to_seconds(double T) {
  return T * secpera;
}

string PISMTime::CF_units_string() {
  return string("seconds since ") + m_reference_date;
}

//! \brief Returns the calendar string.
string PISMTime::calendar() {
  return m_calendar_string;
}

//! \brief Sets the reference date string.
void PISMTime::set_reference_date(string str) {
  m_reference_date = str;
}

void PISMTime::step(double delta_t) {
  m_time_in_seconds += delta_t;

  // If we are less than 0.001 second from the end of the run, reset
  // m_time_in_seconds to avoid taking a very small (and useless) time step.
  if (m_run_end > m_time_in_seconds &&
      m_run_end - m_time_in_seconds < 1e-3)
    m_time_in_seconds = m_run_end;
}

string PISMTime::units_string() {
  return "seconds";
}

bool PISMTime::use_reference_date() {
  return false;
}


PetscErrorCode PISMTime::init() {
  PetscErrorCode ierr;
  bool y_set, ys_set, ye_set;
  PetscReal y = m_config.get("run_length_years"),
    ys = m_config.get("start_year"),
    ye = ys + y;

  ierr = PetscOptionsBegin(m_com, "", "PISM model time options", ""); CHKERRQ(ierr);
  {
    ierr = PISMOptionsReal("-ys", "Start year", ys, ys_set); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-ye", "End year", ye, ye_set); CHKERRQ(ierr);
    ierr = PISMOptionsReal("-y", "Run length, in years", y, y_set); CHKERRQ(ierr);
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
    m_run_start = ys * secpera;
  }

  m_time_in_seconds = m_run_start;

  if (ye_set == true) {
    if (ye < seconds_to_years(m_time_in_seconds)) {
      ierr = PetscPrintf(m_com,
			"PISM ERROR: -ye (%3.3f) is less than -ys (%3.3f) (or input file year or default).\n"
			"PISM cannot run backward in time.\n",
			 ye, seconds_to_years(m_run_start)); CHKERRQ(ierr);
      PISMEnd();
    }
    m_run_end = ye * secpera;
  } else if (y_set == true) {
    m_run_end = m_run_start + y * secpera;
  } else {
    m_run_end = m_run_start + m_config.get("run_length_years", m_unit_system, "years", "seconds");
  }

  return 0;
}

string PISMTime::date(double T) {
  char tmp[256];
  snprintf(tmp, 256, "%3.3f", seconds_to_years(T));
  return string(tmp);
}

string PISMTime::date() {
  return date(m_time_in_seconds);
}

string PISMTime::start_date() {
  return date(m_run_start);
}

string PISMTime::end_date() {
  return date(m_run_end);
}

string PISMTime::run_length() {
  char tmp[256];
  snprintf(tmp, 256, "%3.3f", seconds_to_years(m_run_end - m_run_start));
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

double PISMTime::calendar_year_start(double T) {
  return T - this->mod(T, secpera);
}

double PISMTime::increment_date(double T, int years) {
  return T + secpera * years;
}

PetscErrorCode PISMTime::parse_times(string spec,
                                     vector<double> &result) {
  
  if (spec.find(',') != string::npos) {
    // a list will always contain a comma because at least two numbers are
    // needed to specify reporting intervals

    return parse_list(spec, result);

  } else {
    // it must be a range specification
    return parse_range(spec, result);
  }

  return 0;
}

PetscErrorCode PISMTime::parse_list(string spec, vector<double> &result) {
  vector<string> parts;
  istringstream arg(spec);
  string tmp;

  while(getline(arg, tmp, ','))
    parts.push_back(tmp);

  result.clear();
  for (unsigned int k = 0; k < parts.size(); ++k) {
    double date;
    int errcode = parse_date(parts[k], &date);
    if (errcode != 0)
      return 1;

    result.push_back(date);
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
int PISMTime::parse_interval_length(string spec, string &keyword, double *result) {

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

  // do not allow intervals specified in terms of "fuzzy" units
  if (spec.find("year") != string::npos || spec.find("month") != string::npos) {
    PetscPrintf(m_com, "PISM ERROR: invalid interval length: '%s'\n",
                spec.c_str());
    return 1;
  }

  PISMUnit tmp, seconds, one;
  assert(seconds.parse(m_unit_system, "seconds") == 0);
  assert(one.parse(m_unit_system, "1") == 0);

  // check if the interval spec is a valid unit spec:
  if (tmp.parse(m_unit_system, spec) != 0) {
    PetscPrintf(m_com, "PISM ERROR: invalid interval length: '%s'\n",
                spec.c_str());
    return 1;
  }

  // Check if these units are compatible with "seconds" or "1". The
  // latter allows intervals of the form "0.5", which stands for "half
  // of a model year". This also discards interval specs such as "days
  // since 1-1-1", even though "days" is compatible with "seconds".
  if (ut_are_convertible(tmp.get(), seconds.get()) == 1) {
    cv_converter *c = ut_get_converter(tmp.get(), seconds.get());
    assert(c != NULL);

    if (result)
      *result = cv_convert_double(c, 1.0);

    cv_free(c);

  } else if (ut_are_convertible(tmp.get(), one.get()) == 1) {
    cv_converter *c = ut_get_converter(tmp.get(), one.get());
    assert(c != NULL);

    if (result)
      *result = cv_convert_double(c, 1.0);

    cv_free(c);

  } else {
    PetscPrintf(m_com, "PISM ERROR: invalid interval length: '%s'\n",
                spec.c_str());
    return 1;
  }

  return 0;
}


PetscErrorCode PISMTime::parse_range(string spec, vector<double> &result) {
  double
    time_start   = m_run_start,
    time_end     = m_run_end,
    delta        = 0;
  string keyword = "simple";
  int    errcode;

  if (spec == "hourly") {
    delta = 3600;
  } else if (spec == "daily") {
    delta = 86400;
  } else if (spec == "monthly" || spec == "yearly") {
    keyword = spec;
    delta   = 0;
  } else {
    istringstream  arg(spec);
    vector<string> parts;
    string         tmp;

    while(getline(arg, tmp, ':'))
      parts.push_back(tmp);

    if (parts.size() != 3) {
      PetscPrintf(m_com,
                  "PISM ERROR: A time range must consist of exactly 3 parts, separated by colons. (Got '%s'.)\n",
                  spec.c_str());
      return 1;
    }

    errcode = parse_date(parts[0], &time_start);
    if (errcode != 0)
      return 1;

    errcode = parse_interval_length(parts[1], keyword, &delta);
    if (errcode != 0)
      return 1;

    errcode = parse_date(parts[2], &time_end);
    if (errcode != 0)
      return 1;
    
  }

  return compute_times(time_start, delta, time_end, keyword, result);
}


PetscErrorCode PISMTime::parse_date(string spec, double *result) {
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
    *result = d;

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
PetscErrorCode PISMTime::compute_times_simple(double time_start, double delta, double time_end,
                                              vector<double> &result) {
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
    result.push_back(t);
    k += 1;
    t = time_start + k * delta;
  } while (t < time_end);

  return 0;
}

PetscErrorCode PISMTime::compute_times(double time_start, double delta, double time_end,
                                       string keyword,
                                       vector<double> &result) {
  PetscErrorCode ierr;
  if (keyword != "simple") {
    ierr = PetscPrintf(m_com, "PISM ERROR: only simple time ranges are supported.\n");
    return 1;
  }

  return compute_times_simple(time_start, delta, time_end, result);
}
