/* Copyright (C) 2015 PISM Authors
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

#include "PISMBedDef.hh"

namespace pism {

PBNull::PBNull(const IceGrid &g)
  : BedDef(g) {
  // empty
}

void PBNull::init() {
  BedDef::init();

  verbPrintf(2, m_grid.com,
             "* Initializing the dummy (no-op) bed deformation model...\n");
  m_uplift.set(0.0);
}

void PBNull::update(double t, double dt) {
  m_t  = t;
  m_dt = dt;
  // This model does not update bed topography or bed uplift.
}

} // end of namespace pism
