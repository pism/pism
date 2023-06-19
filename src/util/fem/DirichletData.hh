/* Copyright (C) 2020, 2022 PISM Authors
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

#include "pism/util/fem/FEM.hh"
#include "pism/util/petscwrappers/Mat.hh" // Mat
#include "pism/util/Vector2d.hh"

namespace pism {

namespace array {
class Array;
class Scalar;
class Vector;
}


namespace fem {

class Element2;

//* Parts shared by scalar and 2D vector Dirichlet data classes.
class DirichletData {
public:
  void constrain(Element2 &element);
  operator bool() {
    return m_indices != NULL;
  }
protected:
  DirichletData();
  ~DirichletData();

  void init(const array::Scalar *indices, const array::Array *values, double weight = 1.0);
  void finish(const array::Array *values);

  const array::Scalar *m_indices;
  double m_indices_e[q1::n_chi];
  double m_weight;
};

class DirichletData_Scalar : public DirichletData {
public:
  DirichletData_Scalar(const array::Scalar *indices, const array::Scalar *values,
                       double weight = 1.0);
  ~DirichletData_Scalar();

  void enforce(const Element2 &element, double* x_e);
  void enforce_homogeneous(const Element2 &element, double* x_e);
  void fix_residual(double const *const *const x_global, double **r_global);
  void fix_residual_homogeneous(double **r_global);
  void fix_jacobian(Mat J);
protected:
  const array::Scalar *m_values;
};

class DirichletData_Vector : public DirichletData {
public:
  DirichletData_Vector(const array::Scalar *indices, const array::Vector *values,
                       double weight);
  ~DirichletData_Vector();

  void enforce(const Element2 &element, Vector2d* x_e);
  void enforce_homogeneous(const Element2 &element, Vector2d* x_e);
  void fix_residual(Vector2d const *const *const x_global, Vector2d **r_global);
  void fix_residual_homogeneous(Vector2d **r);
  void fix_jacobian(Mat J);
protected:
  const array::Vector *m_values;
};

} // end of namespace fem
} // end of namespace pism

#endif /* PISM_DIRICHLETDATA_H */
