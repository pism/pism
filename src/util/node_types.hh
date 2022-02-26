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

#ifndef _NODE_TYPES_H_
#define _NODE_TYPES_H_

namespace pism {

//! Node types (in the map plane). These are used to implement boundary conditions in the Q1 FEM
//! code.
/**
   Note: The condition node_type > 0.5 is also used to detect "exterior" (ice-free) nodes where we
   prescribe homogeneous Dirichlet B.C. See DirichletData_*::fix_residual_homogeneous.

   This means that interior and boundary types should not use positive values.
 */
enum NodeType {
  NODE_INTERIOR = -1,
  NODE_BOUNDARY = 0,
  NODE_EXTERIOR = 1
};

namespace array { class Scalar; }

void compute_node_types(const array::Scalar &ice_thickness,
                        double thickness_threshold,
                        array::Scalar &result);

} // end of namespace pism


#endif /* _NODE_TYPES_H_ */
