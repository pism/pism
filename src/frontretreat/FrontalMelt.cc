/* Copyright (C) 2016, 2017, 2018, 2019 PISM Authors
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

#include "FrontalMelt.hh"

#include "pism/geometry/part_grid_threshold_thickness.hh"
#include "pism/geometry/Geometry.hh"

namespace pism {

FrontalMelt::FrontalMelt(IceGrid::ConstPtr grid)
  : Component(grid) {
  // empty
}

FrontalMelt::~FrontalMelt() {
  // empty
}

void FrontalMelt::init() {
  m_log->message(2,
                 "* Initializing the frontal melt parameterization...\n");
}

const IceModelVec2S &FrontalMelt::retreat_rate() const {
  return m_retreat_rate;
}

DiagnosticList FrontalMelt::diagnostics_impl() const {
  return {{"retreat_rate_due_to_frontal_melt", Diagnostic::wrap(m_retreat_rate)}};
}

/*!
 * Convert provided melt rate into the corresponding rate of retreat, considering which
 * part of the front is submerged.
 */
void FrontalMelt::update(const Geometry &geometry, const IceModelVec2S &frontal_melt_rate) {
}

} // end of namespace pism
