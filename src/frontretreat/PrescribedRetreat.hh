/* Copyright (C) 2019, 2022 PISM Authors
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

#ifndef PRESCRIBED_RETREAT_H
#define PRESCRIBED_RETREAT_H

#include "pism/util/Component.hh"
#include "pism/util/array/Forcing.hh"
#include "pism/util/MaxTimestep.hh"

namespace pism {

/*! Implementation of the ISMIP6 parameterized retreat.
 *
 * Reads a space-and-time-dependent ice extent mask from a file. This mask contains values
 * from 0 to 1. Zero corresponds to "ice-free", one to "ice covered", values in between
 * correspond to cells that are partially covered.
 *
 * Each time update() is called, this module gets the mask corresponding to the provided
 * model time. Then, for each grid cell, if the mask is zero, ice is removed (updating ice
 * thickness and area specific volume). If the mask is between 0 and 1, remove the
 * corresponding fraction of ice volume in this cell. If the mask is 1 ice thickness is
 * not modified.
 */
class PrescribedRetreat : public Component {
public:
  PrescribedRetreat(std::shared_ptr<const Grid> grid);
  virtual ~PrescribedRetreat() = default;

  void init();

  void update(double t, double dt,
              array::Scalar &ice_thickness,
              array::Scalar &ice_area_specific_volume);

protected:
  MaxTimestep max_timestep_impl(double t) const;

  std::shared_ptr<array::Forcing> m_retreat_mask;
};

} // end of namespace pism

#endif /* PRESCRIBED_RETREAT_H */
