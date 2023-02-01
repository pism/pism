/* Copyright (C) 2022 PISM Authors
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

#ifndef PISM_GIVEN_CALVING_RATE_H
#define PISM_GIVEN_CALVING_RATE_H

#include "pism/util/Component.hh"
#include "pism/util/array/Forcing.hh"

namespace pism {

class Geometry;

namespace calving {

/*! @brief Calving mechanism removing ice at the shelf front using a provided
    time-dependent spatially-variable calving rate.
 */
class GivenRate : public Component {
public:
  GivenRate(IceGrid::ConstPtr grid);
  virtual ~GivenRate() = default;

  void init();
  void update(double t, double dt);

  const array::Scalar &calving_rate() const;

private:
  DiagnosticList diagnostics_impl() const;

  std::shared_ptr<array::Forcing> m_calving_rate;
};

} // end of namespace calving
} // end of namespace pism

#endif /* PISM_GIVEN_CALVING_RATE_H */
