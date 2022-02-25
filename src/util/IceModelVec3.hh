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

#ifndef PISM_ICEMODELVEC3_H
#define PISM_ICEMODELVEC3_H

#include "pism/util/iceModelVec.hh"

namespace pism {

class IceModelVec2S;

//! \brief A virtual class collecting methods common to ice and bedrock 3D
//! fields.
class IceModelVec3 : public IceModelVec {
public:

  // Three-dimensional array with a number of vertical levels
  IceModelVec3(IceGrid::ConstPtr grid,
               const std::string &name,
               IceModelVecKind ghostedp,
               const std::vector<double> &levels,
               unsigned int stencil_width = 1);

  // A collection of two-dimensional arrays using three-dimensional indexing
  IceModelVec3(IceGrid::ConstPtr grid,
               const std::string &name,
               IceModelVecKind ghostedp,
               unsigned int dof,
               unsigned int stencil_width = 1);

  virtual ~IceModelVec3() = default;

  typedef std::shared_ptr<IceModelVec3> Ptr;
  typedef std::shared_ptr<const IceModelVec3> ConstPtr;

  std::shared_ptr<IceModelVec3> duplicate() const;

  void set_column(int i, int j, double c);
  void set_column(int i, int j, const double *input);
  double* get_column(int i, int j);
  const double* get_column(int i, int j) const;

  double interpolate(int i, int j, double z) const;

  inline double& operator() (int i, int j, int k);
  inline const double& operator() (int i, int j, int k) const;

  void copy_from(const IceModelVec3 &input);
};

void extract_surface(const IceModelVec3 &data, double z, IceModelVec2S &output);
void extract_surface(const IceModelVec3 &data, const IceModelVec2S &z, IceModelVec2S &output);

void sum_columns(const IceModelVec3 &data, double A, double B, IceModelVec2S &output);

inline double& IceModelVec3::operator() (int i, int j, int k) {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(m_array)[j][i][k];
}

inline const double& IceModelVec3::operator() (int i, int j, int k) const {
#if (Pism_DEBUG==1)
  check_array_indices(i, j, k);
#endif
  return static_cast<double***>(m_array)[j][i][k];
}

} // end of namespace pism

#endif /* PISM_ICEMODELVEC3_H */
