/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2020, 2023, 2025, 2026 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * PISM is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PISM; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include <limits>

#include "pism/util/Units.hh"

#include <memory>
#include <string>
#include <udunits2.h>

#include "pism/external/calcalcs/utCalendar2_cal.h"

#include "pism/external/calcalcs/calcalcs.h"

#include "pism/util/error_handling.hh"
#include "pism/util/pism_utilities.hh"

namespace pism {

namespace units {

struct System::Impl {
  Impl(const std::string &path) {
    ut_system *tmp;

    ut_set_error_message_handler(ut_ignore);

    if (not path.empty()) {
      tmp = ut_read_xml(path.c_str());
    } else {
      tmp = ut_read_xml(NULL);
    }

    if (tmp == NULL) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "ut_read_xml(%s) failed", path.c_str());
    }
    ut_set_error_message_handler(ut_write_to_stderr);

    system = tmp;
  }
  ~Impl() {
    ut_free_system(system);
  }
  ut_system *system;
};

/** Initialize the unit system by reading from an XML unit
 * definition file.
 */
System::System(const std::string &path) {
  m_impl.reset(new Impl(path));
}

//! \brief Convert a quantity from unit1 to unit2.
/*!
 * Example: convert(1, "m year-1", "m second-1").
 *
 * Please avoid using in computationally-intensive code.
 */
double convert(System::Ptr system, double input,
               const std::string &spec1, const std::string &spec2) {
  Converter c(Unit(system, spec1), Unit(system, spec2));

  return c(input);
}

struct Unit::Impl {
  Impl(System::Ptr sys, const std::string &spec)
    : system(sys), unit_string(spec) {
    unit = ut_parse(sys->m_impl->system, spec.c_str(), UT_ASCII);

    if (unit == NULL) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "unit specification '%s' is unknown or invalid",
                                    spec.c_str());
    }
  }
  Impl(const Unit::Impl &other) {

    unit = ut_clone(other.unit);
    if (unit == NULL) {
      throw RuntimeError(PISM_ERROR_LOCATION, "ut_clone failed");
    }

    system      = other.system;
    unit_string = other.unit_string;
  }
  ~Impl() {
    reset();
  }
  void reset() {
    ut_free(unit);
    unit_string = "";
  }

  System::Ptr system;
  std::string unit_string;
  ut_unit *unit;
};

Unit::Unit(System::Ptr system, const std::string &spec) {
  m_impl.reset(new Impl(system, spec));
}

Unit::Unit(const Unit &other) {
  m_impl.reset(new Impl(*other.m_impl));
}

Unit& Unit::operator=(const Unit& other) {
  if (this == &other) {
    return *this;
  }

  reset();

  m_impl->system      = other.m_impl->system;
  m_impl->unit_string = other.m_impl->unit_string;

  m_impl->unit = ut_clone(other.m_impl->unit);
  if (m_impl->unit == NULL) {
    throw RuntimeError(PISM_ERROR_LOCATION, "ut_clone failed");
  }

  return *this;
}

bool Unit::is_convertible(const Unit &other) const {
  return ut_are_convertible(m_impl->unit, other.m_impl->unit) != 0;
}

std::string Unit::string() const {
  return m_impl->unit_string;
}

Unit::operator std::string() const {
  return m_impl->unit_string;
}


void Unit::reset() {
  m_impl->reset();
}

System::Ptr Unit::system() const {
  return m_impl->system;
}

DateTime Unit::date(double T, const std::string &calendar) const {

  int year, month, day, hour, minute;
  double second;
  int errcode = utCalendar2_cal(T, m_impl->unit,
                                &year, &month, &day,
                                &hour, &minute, &second,
                                calendar.c_str());
  if (errcode != 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot convert time %f to a date in calendar %s",
                                  T, calendar.c_str());
  }
  return {year, month, day, hour, minute, second, calendar};
}

double Unit::time(const DateTime &d) const {
  double result;

  int errcode = utInvCalendar2_cal(d.year, d.month, d.day,
                                   d.hour, d.minute, d.second,
                                   m_impl->unit,
                                   &result,
                                   d.calendar.c_str());
  if (errcode != 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot convert date and time %d-%d-%d %d:%d:%f to time in calendar %s",
                                  d.year, d.month, d.day,
                                  d.hour, d.minute, d.second,
                                  d.calendar.c_str());
  }
  return result;
}

struct Converter::Impl {
  Impl() {
    converter = cv_get_trivial();
  }
  Impl(System::Ptr sys, const std::string &spec1, const std::string &spec2) {

    Unit u1(sys, spec1), u2(sys, spec2);

    if (not u1.is_convertible(u2)) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "cannot convert '%s' to '%s'",
                                    spec1.c_str(), spec2.c_str());
    }

    converter = ut_get_converter(u1.m_impl->unit, u2.m_impl->unit);
    if (converter == nullptr) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "cannot create a converter from %s to %s",
                                    spec1.c_str(), spec2.c_str());
    }

  }
  Impl(const Unit &u1, const Unit &u2) {
    if (not u1.is_convertible(u2)) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "cannot convert '%s' to '%s'",
                                    u1.string().c_str(), u2.string().c_str());
    }

    converter = ut_get_converter(u1.m_impl->unit, u2.m_impl->unit);
    if (converter == nullptr) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "failed to create a converter from '%s' to '%s'",
                                    u1.string().c_str(), u2.string().c_str());
    }

  }
  ~Impl() {
    cv_free(converter);
    converter = NULL;
  }
  cv_converter *converter;
};

Converter::Converter() {
  m_impl.reset(new Impl());
}

Converter::Converter(System::Ptr sys,
                     const std::string &spec1, const std::string &spec2) {
  m_impl.reset(new Impl(sys, spec1, spec2));
}

Converter::Converter(const Unit &u1, const Unit &u2) {
  m_impl.reset(new Impl(u1, u2));
}

bool are_convertible(const Unit &u1, const Unit &u2) {
  return u1.is_convertible(u2);
}

double Converter::operator()(double input) const {
  return cv_convert_double(m_impl->converter, input);
}

void Converter::convert_doubles(double *data, size_t length) const {
  cv_convert_doubles(m_impl->converter, data, length, data);
}

DateTime::DateTime(int year_, int month_, int day_, int hour_, int minute_, double second_,
                   const std::string &calendar_) {
  year     = year_;
  month    = month_;
  day      = day_;
  hour     = hour_;
  minute   = minute_;
  second   = second_;
  calendar = calendar_;
}

DateTime::DateTime(const std::string &input, const std::string &input_calendar) {

  auto spec = string_strip(input);

  // separate spec into date and time
  std::string date, time;
  {
    auto parts = split(spec, ' ');

    if (parts.empty()) {
      throw RuntimeError::formatted(
          PISM_ERROR_LOCATION,
          "expected a date-time spec consisting of date and time separated by spaces; got '%s'",
          spec.c_str());
    }

    date = parts[0];
    if (parts.size() == 1) {
      time = "0:0:0";
    } else if (parts.size() == 2) {
      time = parts[1];
    } else {
      throw RuntimeError::formatted(
          PISM_ERROR_LOCATION,
          "expected a date-time spec consisting of date and time separated by spaces; got '%s'",
          spec.c_str());
    }
  }

  std::vector<int> date_parts = { 0, 0, 0 };
  std::vector<int> time_parts = { 0, 0, 0 };
  try {
    // date
    {
      // We need to remember if the year was negative in the input string: split() will ignore
      // empty tokens separated by "-", so "-1000-1-1" will produce ["1000", "1", "1"].
      bool year_is_negative = (date[0] == '-');

      if (year_is_negative and not set_member(input_calendar, { "365_day", "360_day", "noleap" })) {
        throw RuntimeError::formatted(
            PISM_ERROR_LOCATION,
            "negative dates such as '%s' are not allowed with the '%s' calendar.\n"
            "Please submit a bug report if this is a problem.",
            date.c_str(), input_calendar.c_str());
      }

      auto parts = split(date, '-');
      if (parts.size() != 3) {
        throw RuntimeError::formatted(
            PISM_ERROR_LOCATION,
            "date string has to contain 3 parts separated by '-' (Y-M-D); got '%s'", spec.c_str());
      }

      for (int k = 0; k < 3; ++k) {
        auto n = parse_integer(parts[k]);

        if (n > std::numeric_limits<int>::max() or n < std::numeric_limits<int>::min()) {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION, "%ld in '%s' does not fit in an 'int'",
                                        n, spec.c_str());
        }

        date_parts[k] = static_cast<int>(n);
      }

      date_parts[0] *= year_is_negative ? -1 : 1;
    }

    // time
    {
      auto parts = split(time, ':');
      if (parts.size() < 2 or parts.size() > 3) {
        throw RuntimeError::formatted(
            PISM_ERROR_LOCATION,
            "time string has to contain 2 or 3 parts separated by ':' (H:M:S); got '%s'",
            time.c_str());
      }

      for (int k = 0; k < (int)parts.size(); ++k) {
        auto n = parse_integer(parts[k]);

        if (n > std::numeric_limits<int>::max() or n < std::numeric_limits<int>::min()) {
          throw RuntimeError::formatted(PISM_ERROR_LOCATION, "%ld in '%s' does not fit in an 'int'",
                                        n, spec.c_str());
        }

        time_parts[k] = static_cast<int>(n);
      }
    }
  } catch (RuntimeError &e) {
    e.add_context("parsing a date-time specification %s", spec.c_str());
    throw;
  }

  {
    year     = date_parts[0];
    month    = date_parts[1];
    day      = date_parts[2];
    hour     = time_parts[0];
    minute   = time_parts[1];
    second   = time_parts[2];
    calendar = input_calendar;
  }

  // Validate the calendar string and the date in this calendar:
  {
    calcalcs_cal *cal = ccs_init_calendar(calendar.c_str());
    if (cal == NULL) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "calendar string '%s' is invalid",
                                    input_calendar.c_str());
    }

    int jday_dummy = 0;
    int errcode    = ccs_date2jday(cal, year, month, day, &jday_dummy);
    if (errcode != 0) {
      ccs_free_calendar(cal);
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "date %s is invalid in the %s calendar",
                                    spec.c_str(), calendar.c_str());
    }
    ccs_free_calendar(cal);
  }

  // Validate hour
  if (hour < 0 or hour > 23) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid hour in the date-time spec '%s'",
                                  spec.c_str());
  }
  // Validate minute
  if (minute < 0 or minute > 59) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid minute in the date-time spec '%s'",
                                  spec.c_str());
  }
  // Validate second
  if (second < 0.0 or second >= 60.0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid second in the date-time spec '%s'",
                                  spec.c_str());
  }
}


} // end of namespace units

} // end of namespace pism
