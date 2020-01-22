/* Copyright (C) 2013, 2014, 2015, 2016, 2017, 2018 PISM Authors
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

#ifndef _PISMICEBERGREMOVER_H_
#define _PISMICEBERGREMOVER_H_

#include "pism/util/Component.hh"
#include "pism/util/iceModelVec.hh"

namespace pism {

class IceModelVec2CellType;

namespace calving {

/*! \brief PISM iceberg remover */
/*!
 * Identifies and removes free-floating icebergs, which cause
 * well-posedness problems for stress solvers.
 *
 * Icebergs are, in this context, floating regions that are _not_
 * attached, through a chain of positive thickness ice-filled cells,
 * to at least one grounded cell.
 *
 * They cause the SSA operator to have a nontrivial null space.
 *
 * They are observed to cause unrealistically large velocities that
 * may affect ice velocities elsewhere.
 *
 * This class uses a serial connected component labeling algorithm to
 * remove "icebergs".
 */
class IcebergRemover : public Component
{
public:
  IcebergRemover(IceGrid::ConstPtr g);
  virtual ~IcebergRemover();

  virtual void init();
  void update(const IceModelVec2Int &bc_mask,
              IceModelVec2CellType &pism_mask,
              IceModelVec2S &ice_thickness);
protected:
  IceModelVec2S m_iceberg_mask;
  petsc::Vec::Ptr m_mask_p0;
};

} // end of namespace calving
} // end of namespace pism

#endif /* _PISMICEBERGREMOVER_H_ */
