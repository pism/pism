/* Copyright (C) 2015, 2016, 2021 PISM Authors
 *
 * This file is part of PISM.
 *
 * PISM is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 3 of the License, or (at your option) any later
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

#ifndef PISM_MAXTIMESTEP_HH
#define PISM_MAXTIMESTEP_HH

#include <string>

namespace pism {

//! @brief Combines the max. time step with the flag indicating if a
//! restriction is active. Makes is possible to use the less-than
//! operator (and `std::min`, etc) to choose the stricter of two
//! restrictions.
class MaxTimestep {
public:
  //! @brief Create an instance corresponding to an "inactive"
  //! time-step restriction (in other words, @f$ \Delta t = \infty @f$).
  MaxTimestep();
  //! Create an instance corresponding to a max time step of `value`.
  MaxTimestep(double value);

  MaxTimestep(const std::string &new_description);
  //! Create an instance and provide a description.
  MaxTimestep(double value, const std::string &new_description);
  //! Convert to `bool` to check if a time step restriction is "active".
  bool finite() const;
  bool infinite() const;

  //! Get the value of the maximum time step.
  double value() const;

  std::string description() const;
private:
  bool m_finite;
  double m_value;
  std::string m_description;
};

//! Greater than operator for MaxTimestep.
bool operator>(const MaxTimestep &a, const MaxTimestep &b);

//! Less than operator for MaxTimestep.
bool operator<(const MaxTimestep &a, const MaxTimestep &b);

//! Equality operator for MaxTimestep.
bool operator==(const MaxTimestep &a, const MaxTimestep &b);

} // end of namespace pism

#endif /* PISM_MAXTIMESTEP_HH */
