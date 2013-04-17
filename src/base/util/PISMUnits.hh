/* Copyright (C) 2013 PISM Authors
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

#ifndef _PISMUNIT_H_
#define _PISMUNIT_H_

#include <udunits2.h>
#include <string>
#include <tr1/memory>

/** @file PISMUnits.hh This file contains thin wrappers around
 * UDUNITS-2 objects. Nothing fancy. The only purpose is to simplify
 * memory management for objects that are stored as data members other
 * C++ classes.
 *
 * The `cv_converter` object is *not* wrapped, because we deallocate
 * these right away.
 *
 * One thing is worth mentioning, though: in UDUNITS-2, every ut_unit
 * object contains a pointer to the unit system that was used to create it.
 */

class PISMUnitSystem {
  friend class PISMUnit;
public:
  PISMUnitSystem();
  PISMUnitSystem(const char *path);
  typedef std::tr1::shared_ptr<ut_system> Ptr;

  PISMUnitSystem::Ptr get() const;
private:
  PISMUnitSystem::Ptr m_system;
};

class PISMUnit {
public:
  PISMUnit();
  PISMUnit(const PISMUnit &other);
  ~PISMUnit();

  void reset();

  PISMUnit& operator=(const PISMUnit& other);
  int parse(PISMUnitSystem system, std::string spec);
  std::string format() const;

  ut_unit* get() const;
  PISMUnitSystem get_system() const;

  bool is_valid() const;
private:
  ut_unit *m_unit;
  PISMUnitSystem m_system;
};

#endif /* _PISMUNIT_H_ */
