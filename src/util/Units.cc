/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018, 2020, 2023 PISM Authors
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

#include "pism/util/Units.hh"

#include <udunits2.h>

#include "pism/external/calcalcs/utCalendar2_cal.h"

#include "pism/util/error_handling.hh"

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

std::string Unit::format() const {
  return m_impl->unit_string;
}

void Unit::reset() {
  m_impl->reset();
}

System::Ptr Unit::system() const {
  return m_impl->system;
}

DateTime Unit::date(double T, const std::string &calendar) const {
  DateTime result;

  int errcode = utCalendar2_cal(T, m_impl->unit,
                                &result.year, &result.month, &result.day,
                                &result.hour, &result.minute, &result.second,
                                calendar.c_str());
  if (errcode != 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot convert time %f to a date in calendar %s",
                                  T, calendar.c_str());
  }
  return result;
}

double Unit::time(const DateTime &d, const std::string &calendar) const {
  double result;

  int errcode = utInvCalendar2_cal(d.year, d.month, d.day,
                                   d.hour, d.minute, d.second,
                                   m_impl->unit,
                                   &result,
                                   calendar.c_str());
  if (errcode != 0) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "cannot convert date and time %d-%d-%d %d:%d:%f to time in calendar %s",
                                  d.year, d.month, d.day,
                                  d.hour, d.minute, d.second,
                                  calendar.c_str());
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
    if (not converter) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "cannot create a converter from %s to %s",
                                    spec1.c_str(), spec2.c_str());
    }

  }
  Impl(const Unit &u1, const Unit &u2) {
    if (not u1.is_convertible(u2)) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION, "cannot convert '%s' to '%s'",
                                    u1.format().c_str(), u2.format().c_str());
    }

    converter = ut_get_converter(u1.m_impl->unit, u2.m_impl->unit);
    if (not converter) {
      throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                    "failed to create a converter from '%s' to '%s'",
                                    u1.format().c_str(), u2.format().c_str());
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

} // end of namespace units

} // end of namespace pism
