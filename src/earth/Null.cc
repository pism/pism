/* Copyright (C) 2015, 2016, 2017, 2018, 2023, 2024 PISM Authors
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

#include "pism/earth/BedDef.hh"
#include "pism/util/Grid.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {
namespace bed {

Null::Null(std::shared_ptr<const Grid> g)
  : BedDef(g, "dummy (no-op)") {
  // empty
}

void Null::init_impl(const InputOptions & /*opts*/, const array::Scalar & /*ice_thickness*/,
                     const array::Scalar & /*sea_level_elevation*/) {
  m_uplift.set(0.0);
}

void Null::bootstrap_impl(const array::Scalar & /*bed_elevation*/,
                           const array::Scalar & /*bed_uplift*/,
                           const array::Scalar & /*ice_thickness*/,
                           const array::Scalar & /*sea_level_elevation*/) {
  // empty
}

void Null::update_impl(const array::Scalar &/*load*/,
                       double /*t*/, double /*dt*/) {
  // This model does not update bed topography
}

MaxTimestep Null::max_timestep_impl(double /*t*/) const {
  return {};
}

} // end of namespace bed
} // end of namespace pism
