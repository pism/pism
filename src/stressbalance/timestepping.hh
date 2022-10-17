/* Copyright (C) 2016, 2017, 2022 PISM Authors
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

#ifndef TIMESTEPPING_H
#define TIMESTEPPING_H

#include "pism/util/MaxTimestep.hh"

namespace pism {


namespace array {
class Array3D;
class CellType1;
class CellType2;
class CellType;
class Scalar;
class Vector;
} // end of namespace array

struct CFLData {
  CFLData();
  MaxTimestep dt_max;
  double u_max, v_max, w_max;
};

/*! @brief Compute the max. time step according to the CFL condition (within the volume of the
    ice). */
/*!
 * Returns the maximum time step along with maximum speeds along x, y, and z directions within the
 * ice. Note that PISM uses semi-implicit methods in energy balance and age models, so this code
 * does *not* use the w component of the velocity field in the computation of the max. time step.
 * The maximum of the speed along the z axis is computed for reporting.
 */
CFLData max_timestep_cfl_3d(const array::Scalar &ice_thickness,
                            const array::CellType &cell_type,
                            const array::Array3D &u3,
                            const array::Array3D &v3,
                            const array::Array3D &w3);

/*! @brief Compute the max. time step according to the CFL condition (within the ice, 2D). */
/*!
 * Returns the maximum time step along with maximum speeds along x and y directions within the
 * ice.
 */
CFLData max_timestep_cfl_2d(const array::Scalar &ice_thickness,
                            const array::CellType &cell_type,
                            const array::Vector &velocity);

} // end of namespace pism


#endif /* TIMESTEPPING_H */
