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

#ifndef _PISMUNITS_H_
#define _PISMUNITS_H_

#include <udunits2.h>
#include <string>
#include <petscvec.h>

#ifdef PISM_USE_TR1
#include <tr1/memory>
#else
#include <memory>
#endif

namespace pism {

/** @file PISMUnits.hh This file contains thin wrappers around
 * UDUNITS-2 objects. Nothing fancy. The only purpose is to simplify
 * memory management for objects that are stored as data members of
 * other C++ classes.
 *
 * The `cv_converter` object is *not* wrapped, because we deallocate
 * these right away.
 *
 * One thing is worth mentioning, though: in UDUNITS-2, every ut_unit
 * object contains a pointer to the unit system that was used to create it.
 *
 * We use C++ shared pointers to make sure that the system a
 * Unit instance needs is allocated during the whole life span of
 * this instance. (De-allocating the unit system too early results in
 * having a "dangling" pointer.)
 */

class UnitSystem {
  friend class Unit;
public:
  UnitSystem(const char *path);
#ifdef PISM_USE_TR1
  typedef std::tr1::shared_ptr<ut_system> Ptr;
#else
  typedef std::shared_ptr<ut_system> Ptr;
#endif

  UnitSystem::Ptr get() const;

  double convert(double input, const std::string &spec1, const std::string &spec2) const;
private:
  UnitSystem::Ptr m_system;
};

class Unit {
public:
  Unit(const UnitSystem &system);
  Unit(const Unit &other);
  ~Unit();

  void reset();

  Unit& operator=(const Unit& other);
  int parse(const std::string &spec);
  std::string format() const;

  ut_unit* get() const;
  UnitSystem get_system() const;
  cv_converter* get_converter_from(const Unit &from) const;

  bool is_valid() const;
private:
  ut_unit *m_unit;
  UnitSystem m_system;
  std::string m_unit_string;
};

/**
 * Check if two units are convertible.
 *
 * @param[in] from source units
 * @param[in] to destination units
 *
 * @return true if units are convertible, false otherwise
 */
bool units_are_convertible(Unit from, Unit to);

PetscErrorCode units_check(std::string name, Unit from, Unit to);

/**
 * Convert a PETSc Vec from the units in `from` into units in `to` (in place).
 *
 * @param[in,out] v data
 * @param[in] from source units
 * @param[in] to destination units
 *
 * @return 0 on success
 */
PetscErrorCode convert_vec(Vec v, Unit from, Unit to);

/**
 * Convert 
 *
 * @param[in,out] data data to convert
 * @param[in] length number of elements in `data`
 * @param[in] from source units
 * @param[in] to destination units
 *
 * @return 0 on success
 */
PetscErrorCode convert_doubles(double *data, size_t length, Unit from, Unit to);

} // end of namespace pism

#endif /* _PISMUNITS_H_ */
