/* Copyright (C) 2022, 2023 PISM Authors
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

#ifndef PISM_STAGGERED_H
#define PISM_STAGGERED_H

#include "pism/util/array/Array3D.hh"

namespace pism {

namespace array {
class CellType1;
class Vector;

//! @brief A class for storing and accessing internal staggered-grid 2D fields.
//! Uses dof=2 storage. This class is identical to array::Vector, except that
//! components are not called `u` and `v` (to avoid confusion).
class Staggered : public Array {
public:
  Staggered(IceGrid::ConstPtr grid, const std::string &name);

  typedef std::shared_ptr<array::Staggered> Ptr;
  typedef std::shared_ptr<const array::Staggered> ConstPtr;

  inline double& operator() (int i, int j, int k);
  inline const double& operator() (int i, int j, int k) const;

  void copy_from(const array::Staggered &input);
protected:
  Staggered(IceGrid::ConstPtr grid, const std::string &name,
            unsigned int stencil_width);
};

inline double& array::Staggered::operator() (int i, int j, int k) {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(m_array)[j][i][k];
}

inline const double& array::Staggered::operator() (int i, int j, int k) const {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(m_array)[j][i][k];
}

class Staggered1 : public Staggered {
public:
  Staggered1(IceGrid::ConstPtr grid, const std::string &name);

  //! Returns the values at interfaces of the cell i,j using the staggered grid.
  /*! The ij member of the return value is set to 0, since it has no meaning in
    this context.
  */
  inline stencils::Star<double> star(int i, int j) const;
};

inline stencils::Star<double> Staggered1::star(int i, int j) const {
  const Staggered1 &self = *this;

  stencils::Star<double> result;

  result.c = 0.0;               // has no meaning in this context
  result.e = self(i, j, 0);
  result.w = self(i-1, j, 0);
  result.n = self(i, j, 1);
  result.s = self(i, j-1, 1);

  return result;
}

/*!
 * Computes maximums of absolute values of both components.
 */
std::array<double,2> absmax(const array::Staggered &input);

/*!
 * Average a scalar field from the staggered grid onto the regular grid by considering
 * only ice-covered grid.
 *
 * If `include_floating_ice` is true, include floating ice, otherwise consider grounded
 * icy cells only.
 */
void staggered_to_regular(const array::CellType1 &cell_type,
                          const array::Staggered1 &input,
                          bool include_floating_ice,
                          array::Scalar &result);

/*!
 * Average a vector field from the staggered grid onto the regular grid by considering
 * only ice-covered grid.
 *
 * If `include_floating_ice` is true, include floating ice, otherwise consider grounded
 * icy cells only.
 */
void staggered_to_regular(const array::CellType1 &cell_type,
                          const array::Staggered1 &input,
                          bool include_floating_ice,
                          array::Vector &result);

} // end of namespace array

} // end of namespace pism

#endif /* PISM_STAGGERED_H */
