/* Copyright (C) 2013, 2014 PISM Authors
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

#include "PISMUnits.hh"
#include <gsl/gsl_math.h>       // GSL_NAN
#include <petscsys.h>
#include <cstdio>

#include "pism_const.hh"

namespace pism {

class ut_system_deleter {
public:
  void operator()(ut_system* p) const {
    ut_free_system(p);
  }
};

/** Initialize the unit system by reading from an XML unit
 * definition file.
 */
UnitSystem::UnitSystem(const char *path) {
  ut_system *tmp;
  ut_set_error_message_handler(ut_ignore);
  tmp = ut_read_xml(path);
  if (tmp == NULL) {
    PetscPrintf(PETSC_COMM_SELF, "PISM ERROR: UDUNITS-2 initialization failed.\n");
    PetscEnd();
  }
  ut_set_error_message_handler(ut_write_to_stderr);

  m_system = UnitSystem::Ptr(tmp, ut_system_deleter());
}

UnitSystem::Ptr UnitSystem::get() const {
  return m_system;
}

//! \brief Convert a quantity from unit1 to unit2.
/*!
 * Example: convert(1, "m/year", "m/s").
 *
 * Please avoid using in computationally-intensive code.
 */
double UnitSystem::convert(double input, const std::string &spec1, const std::string &spec2) const {
  Unit unit1(*this), unit2(*this);

  if (unit1.parse(spec1) != 0) {
#if (PISM_DEBUG==1)
    fprintf(stderr, "UnitSystem::convert() failed trying to parse %s\n", spec1.c_str());
#endif
    return GSL_NAN;
  }

  if (unit2.parse(spec2) != 0) {
#if (PISM_DEBUG==1)
    fprintf(stderr, "UnitSystem::convert() failed trying to parse %s\n", spec2.c_str());
#endif
    return GSL_NAN;
  }

  cv_converter *c = unit2.get_converter_from(unit1);
  if (c == NULL) {
#if (PISM_DEBUG==1)
    fprintf(stderr, "UnitSystem::convert() failed trying to convert %s to %s\n",
            spec1.c_str(), spec2.c_str());
#endif
    return GSL_NAN;
  }

  double result = cv_convert_double(c, input);
  cv_free(c);

  return result;
}

Unit::Unit(const UnitSystem &system)
  : m_unit(NULL), m_system(system) {
  this->parse("1");
}

Unit::Unit(const Unit &other)
  : m_system(other.m_system) {
  if (other.m_unit == NULL)
    m_unit = NULL;
  else
    m_unit = ut_clone(other.m_unit);

  m_system      = other.m_system;
  m_unit_string = other.m_unit_string;
}

Unit& Unit::operator=(const Unit& other) {
  if (this == &other)
    return *this;

  reset();

  m_system      = other.m_system;
  m_unit_string = other.m_unit_string;

  if (other.m_unit == NULL)
    m_unit = NULL;
  else
    m_unit = ut_clone(other.m_unit);

  return *this;
}

Unit::~Unit() {
  reset();
}

int Unit::parse(const std::string &spec) {
  reset();
  m_unit = ut_parse(m_system.get().get(), spec.c_str(), UT_ASCII);
  m_unit_string = spec;
  if (m_unit == NULL)
    return 1;
  else
    return 0;
}

cv_converter* Unit::get_converter_from(const Unit &from) const {
  return ut_get_converter(from.get(), this->get());
}


std::string Unit::format() const {
  return m_unit_string;
}

void Unit::reset() {
  ut_free(m_unit);
  m_unit = NULL;
}

ut_unit* Unit::get() const {
  return m_unit;
}

UnitSystem Unit::get_system() const {
  return m_system;
}

bool units_are_convertible(Unit from, Unit to) {
  return ut_are_convertible(from.get(), to.get()) != 0;
}

//! Check if provided units are convertible and terminate PISM if they are not.
PetscErrorCode units_check(std::string name, Unit from, Unit to) {

  if (units_are_convertible(from, to) == false) {              // can't convert
    PetscPrintf(PETSC_COMM_SELF,
                "PISM ERROR: processing variable '%s': attempted to convert data from '%s' to '%s'.\n",
                name.c_str(), from.format().c_str(), to.format().c_str());
    PISMEnd();
  }
  return 0;
}

bool Unit::is_valid() const {
  return m_unit != NULL;
}

PetscErrorCode convert_vec(Vec v, Unit from, Unit to) {
  PetscErrorCode ierr;

  int data_size = 0;
  ierr = VecGetLocalSize(v, &data_size); CHKERRQ(ierr);

  double *data = NULL;
  ierr = VecGetArray(v, &data); CHKERRQ(ierr);
  ierr = convert_doubles(data, data_size, from, to); CHKERRQ(ierr);
  ierr = VecRestoreArray(v, &data); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode convert_doubles(double *data, size_t length, Unit from, Unit to) {
  std::string from_name, to_name;

  // Get string representations of units:
  from_name = from.format();
  to_name   = to.format();

  // Get the converter:
  cv_converter *c = to.get_converter_from(from);
  if (c == NULL) {              // can't convert
    SETERRQ2(MPI_COMM_SELF, 1, "cannot convert data from '%s' to '%s'.\n",
             from_name.c_str(), to_name.c_str());
  }

  cv_convert_doubles(c, data, length, data);
  cv_free(c);

  return 0;
}

} // end of namespace pism
