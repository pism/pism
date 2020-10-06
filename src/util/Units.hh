/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2020 PISM Authors
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

#include <string>
#include <memory>

namespace pism {

namespace units {

/** @file Units.hh This file contains thin wrappers around
 * UDUNITS-2 objects. Nothing fancy. The only purpose is to simplify
 * memory management for objects that are stored as data members of
 * other C++ classes.
 *
 * One thing is worth mentioning, though: in UDUNITS-2, every ut_unit
 * object contains a pointer to the unit system that was used to create it.
 *
 * We use C++ shared pointers to make sure that the system a
 * Unit instance needs is allocated during the whole life span of
 * this instance. (De-allocating the unit system too early results in
 * having a "dangling" pointer.)
 */

class System {
public:
  System(const std::string &path = "");
  typedef std::shared_ptr<System> Ptr;
private:
  friend class Unit;

  struct Impl;
  std::shared_ptr<Impl> m_impl;

  System(const System &);
  System& operator=(System const &);
};

double convert(System::Ptr system, double input,
               const std::string &spec1, const std::string &spec2);

struct DateTime {
  int year, month, day, hour, minute;
  double second;
};

class Unit {
public:
  Unit(System::Ptr system, const std::string &spec);
  Unit(const Unit &other);

  bool is_convertible(const Unit &other) const;

  DateTime date(double T, const std::string &calendar) const;
  double time(const DateTime &d, const std::string &calendar) const;

  Unit& operator=(const Unit& other);
  std::string format() const;

  System::Ptr system() const;
private:
  friend class Converter;
  void reset();

  struct Impl;
  std::shared_ptr<Impl> m_impl;
};

/** Check if units are convertible without creating a converter.
 *
 * @param[in] u1 first Unit instance
 * @param[in] u2 second Unit instance
 *
 * @return true if units are convertible, false otherwise
 */
bool are_convertible(const Unit &u1, const Unit &u2);

/** Unit converter.
 * 
 * Throws pism::RuntimeError() if the conversion is not possible.
 *
 */
class Converter {
public:
  Converter();
  Converter(const Unit &u1, const Unit &u2);
  Converter(System::Ptr sys, const std::string &u1, const std::string &u2);
  /** Convert an array of doubles in place
   *
   * @param[in,out] data array to process
   * @param length length of the array
   */
  void convert_doubles(double *data, size_t length) const;
  double operator()(double input) const;
private:

  struct Impl;
  std::shared_ptr<Impl> m_impl;

  // hide copy constructor and the assignment operator
  Converter(const Converter &);
  Converter& operator=(Converter const &);
};

} // end of namespace units

} // end of namespace pism

#endif /* _PISMUNITS_H_ */
