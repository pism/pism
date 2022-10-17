/* Copyright (C) 2016 PISM Authors
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

#ifndef BTU_MINIMAL_H
#define BTU_MINIMAL_H

#include "BedThermalUnit.hh"

namespace pism {
namespace energy {

class BTU_Minimal : public BedThermalUnit {
public:
  BTU_Minimal(IceGrid::ConstPtr g);

protected:
  void init_impl(const InputOptions &opts);

  double vertical_spacing_impl() const;
  double depth_impl() const;
  unsigned int Mz_impl() const;
  MaxTimestep max_timestep_impl(double t) const;

  using BedThermalUnit::update_impl;
  void update_impl(const array::Scalar &bedrock_top_temperature, double t, double dt);
};

} // end of namespace energy
} // end of namespace pism


#endif /* BTU_MINIMAL_H */
