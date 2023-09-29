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

#ifndef PISM_ARRAY3DCOLLECTION_H
#define PISM_ARRAY3DCOLLECTION_H

#include "pism/util/array/Array.hh"

namespace pism {
namespace array {

//! \brief A class for storing *collections* of fields in a 3D array similar to Array3D.
//
// The only difference is that the last dimension does not correspond to the "z"
// dimension.
class Array3DCollection : public Array {
public:

  // Three-dimensional array with a number of vertical levels
  Array3DCollection(std::shared_ptr<const Grid> grid,
          const std::string &name,
          Kind ghostedp,
          const std::vector<double> &levels,
          unsigned int stencil_width = 1);

  virtual ~Array3DCollection() = default;

  std::shared_ptr<Array3DCollection> duplicate() const;

  double* column(int i, int j);
  const double* column(int i, int j) const;

  void copy_from(const Array3DCollection &input);
private:
  void regrid_impl(const File &file, io::RegriddingFlag flag, double default_value = 0.0);
};

} // end of namespace array
} // end of namespace pism

#endif /* PISM_ARRAY3DCOLLECTION_H */
