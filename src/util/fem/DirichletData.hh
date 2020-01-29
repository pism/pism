/* Copyright (C) 2020 PISM Authors
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
#ifndef PISM_DIRICHLETDATA_H
#define PISM_DIRICHLETDATA_H

#include "FEM.hh"
#include "pism/util/petscwrappers/Mat.hh" // Mat
#include "pism/util/Vector2.hh"

namespace pism {

class IceModelVec;
class IceModelVec2Int;
class IceModelVec2S;
class IceModelVec2V;

namespace fem {

class Element;

//* Parts shared by scalar and 2D vector Dirichlet data classes.
class DirichletData {
public:
  void constrain(Element &element);
  operator bool() {
    return m_indices != NULL;
  }
protected:
  DirichletData();
  ~DirichletData();

  void init(const IceModelVec2Int *indices, const IceModelVec *values, double weight = 1.0);
  void finish(const IceModelVec *values);

  const IceModelVec2Int *m_indices;
  double m_indices_e[q1::n_chi];
  double m_weight;
};

class DirichletData_Scalar : public DirichletData {
public:
  DirichletData_Scalar(const IceModelVec2Int *indices, const IceModelVec2S *values,
                       double weight = 1.0);
  ~DirichletData_Scalar();

  void enforce(const Element &element, double* x_e);
  void enforce_homogeneous(const Element &element, double* x_e);
  void fix_residual(double const *const *const x_global, double **r_global);
  void fix_residual_homogeneous(double **r_global);
  void fix_jacobian(Mat J);
protected:
  const IceModelVec2S *m_values;
};

class DirichletData_Vector : public DirichletData {
public:
  DirichletData_Vector(const IceModelVec2Int *indices, const IceModelVec2V *values,
                       double weight);
  ~DirichletData_Vector();

  void enforce(const Element &element, Vector2* x_e);
  void enforce_homogeneous(const Element &element, Vector2* x_e);
  void fix_residual(Vector2 const *const *const x_global, Vector2 **r_global);
  void fix_residual_homogeneous(Vector2 **r);
  void fix_jacobian(Mat J);
protected:
  const IceModelVec2V *m_values;
};

} // end of namespace fem
} // end of namespace pism

#endif /* PISM_DIRICHLETDATA_H */
