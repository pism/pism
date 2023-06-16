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

#ifndef PISM_ARRAY3D_H
#define PISM_ARRAY3D_H

#include "Array.hh"

namespace pism {

namespace array {

class Scalar;

//! \brief A virtual class collecting methods common to ice and bedrock 3D
//! fields.
class Array3D : public Array {
public:

  // Three-dimensional array with a number of vertical levels
  Array3D(std::shared_ptr<const Grid> grid,
          const std::string &name,
          Kind ghostedp,
          const std::vector<double> &levels,
          unsigned int stencil_width = 1);

  virtual ~Array3D() = default;

  typedef std::shared_ptr<Array3D> Ptr;
  typedef std::shared_ptr<const Array3D> ConstPtr;

  std::shared_ptr<Array3D> duplicate() const;

  void set_column(int i, int j, double c);
  void set_column(int i, int j, const double *input);
  double* get_column(int i, int j);
  const double* get_column(int i, int j) const;

  double interpolate(int i, int j, double z) const;

  inline double& operator() (int i, int j, int k);
  inline const double& operator() (int i, int j, int k) const;

  void copy_from(const Array3D &input);
};

void extract_surface(const Array3D &data, double z, Scalar &output);
void extract_surface(const Array3D &data, const Scalar &z, Scalar &output);

void sum_columns(const Array3D &data, double A, double B, Scalar &output);

inline double& Array3D::operator() (int i, int j, int k) {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(m_array)[j][i][k];
}

inline const double& Array3D::operator() (int i, int j, int k) const {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(m_array)[j][i][k];
}

} // end of namespace array
} // end of namespace pism

#endif /* PISM_ARRAY3D_H */
