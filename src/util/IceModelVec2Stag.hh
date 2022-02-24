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

#ifndef PISM_ICEMODELVEC2STAG_H
#define PISM_ICEMODELVEC2STAG_H

#include "pism/util/IceModelVec3.hh"

namespace pism {

class IceModelVec2V;
class CellTypeArray1;

//! \brief A class for storing and accessing internal staggered-grid 2D fields.
//! Uses dof=2 storage. This class is identical to IceModelVec2V, except that
//! components are not called `u` and `v` (to avoid confusion).
class IceModelVec2Stag : public IceModelVec3 {
public:
  IceModelVec2Stag(IceGrid::ConstPtr grid, const std::string &name,
                   IceModelVecKind ghostedp, unsigned int stencil_width = 1);

  typedef std::shared_ptr<IceModelVec2Stag> Ptr;
  typedef std::shared_ptr<const IceModelVec2Stag> ConstPtr;

  //! Returns the values at interfaces of the cell i,j using the staggered grid.
  /*! The ij member of the return value is set to 0, since it has no meaning in
    this context.
  */
  inline stencils::Star<double> star(int i, int j) const;
};

inline stencils::Star<double> IceModelVec2Stag::star(int i, int j) const {
  const IceModelVec2Stag &self = *this;

  stencils::Star<double> result;

  result.ij = 0.0;             // has no meaning in this context
  result.e =  self(i, j, 0);
  result.w =  self(i-1, j, 0);
  result.n =  self(i, j, 1);
  result.s =  self(i, j-1, 1);

  return result;
}

std::array<double,2> absmax(const IceModelVec2Stag &input);

/*!
 * Average a scalar field from the staggered grid onto the regular grid by considering
 * only ice-covered grid.
 *
 * If `include_floating_ice` is true, include floating ice, otherwise consider grounded
 * icy cells only.
 */
void staggered_to_regular(const CellTypeArray1 &cell_type,
                          const IceModelVec2Stag &input,
                          bool include_floating_ice,
                          IceModelVec2S &result);

/*!
 * Average a vector field from the staggered grid onto the regular grid by considering
 * only ice-covered grid.
 *
 * If `include_floating_ice` is true, include floating ice, otherwise consider grounded
 * icy cells only.
 */
void staggered_to_regular(const CellTypeArray1 &cell_type,
                          const IceModelVec2Stag &input,
                          bool include_floating_ice,
                          IceModelVec2V &result);

} // end of namespace pism

#endif /* PISM_ICEMODELVEC2STAG_H */
