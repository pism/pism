/* Copyright (C) 2013, 2014 PISM Authors
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

#include "PISMOcean.hh"
#include "iceModelVec.hh"

namespace pism {

/** Set `result` to the melange back pressure fraction.
 *
 * This default implementation sets `result` to 0.0.
 *
 * @param[out] result back pressure fraction
 *
 * @return 0 on success
 */
PetscErrorCode PISMOceanModel::melange_back_pressure_fraction(IceModelVec2S &result) {
  PetscErrorCode ierr;
  ierr = result.set(0.0); CHKERRQ(ierr);
  return 0;
}

} // end of namespace pism
