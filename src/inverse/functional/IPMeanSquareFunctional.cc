// Copyright (C) 2012, 2014  David Maxwell
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#include "IPMeanSquareFunctional.hh"

namespace pism {

//! Implicitly set the normalization constant for the functional.
/*! The normalization constant is selected so that if an input
IceModelVec2V has component vectors all of length \a scale, then the funtional value will be 1. I.e.
\f[
c_N^{-1} = \sum_{i} w_i {\tt scale}^2.
\f]*/
PetscErrorCode IPMeanSquareFunctional2V::normalize(double scale) {
  PetscErrorCode   ierr;

  // The local value of the weights
  double value = 0;

  if (m_weights) {
    IceModelVec::AccessList list;
    list.add(*m_weights);
    for (Points p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      value += (*m_weights)(i, j);
    }
  } else {
    for (Points p(m_grid); p; p.next()) {
      value += 1;
    }
  }

  ierr = GlobalSum(m_grid.com, &value,  &m_normalization); CHKERRQ(ierr);
  m_normalization *= (scale*scale);
  return 0;
}

PetscErrorCode IPMeanSquareFunctional2V::valueAt(IceModelVec2V &x, double *OUTPUT)  {
  PetscErrorCode   ierr;

  // The value of the objective
  double value = 0;

  IceModelVec::AccessList list;
  list.add(x);

  if (m_weights) {
    list.add(*m_weights);
    for (Points p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      Vector2 &x_ij = x(i, j);
      value += (x_ij.u*x_ij.u + x_ij.v*x_ij.v)*(*m_weights)(i, j);
    }
  } else {
    for (Points p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      Vector2 &x_ij = x(i, j);
      value += (x_ij.u*x_ij.u + x_ij.v*x_ij.v);
    }
  }
  value /= m_normalization;

  ierr = GlobalSum( m_grid.com, &value,  OUTPUT); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IPMeanSquareFunctional2V::dot(IceModelVec2V &a, IceModelVec2V &b, double *OUTPUT)  {
  PetscErrorCode   ierr;

  // The value of the objective
  double value = 0;

  IceModelVec::AccessList list;
  list.add(a);
  list.add(b);

  if (m_weights) {
    list.add(*m_weights);
    for (Points p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      Vector2 &a_ij = a(i, j);
      Vector2 &b_ij = b(i, j);
      value += (a_ij.u*b_ij.u + a_ij.v*b_ij.v)*(*m_weights)(i, j);
    }
  } else {
    for (Points p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      Vector2 &a_ij = a(i, j);
      Vector2 &b_ij = b(i, j);
      value += (a_ij.u*b_ij.u + a_ij.v*b_ij.v);
    }
  }
  value /= m_normalization;

  ierr = GlobalSum( m_grid.com, &value,  OUTPUT); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IPMeanSquareFunctional2V::gradientAt(IceModelVec2V &x, IceModelVec2V &gradient)  {
  gradient.set(0);

  IceModelVec::AccessList list;
  list.add(x);
  list.add(gradient);

  if (m_weights) {
    list.add(*m_weights);
    for (Points p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      gradient(i, j).u = 2*x(i, j).u*(*m_weights)(i, j) / m_normalization;
      gradient(i, j).v = 2*x(i, j).v*(*m_weights)(i, j) / m_normalization;
    }
  } else {
    for (Points p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      gradient(i, j).u = 2*x(i, j).u / m_normalization;
      gradient(i, j).v = 2*x(i, j).v / m_normalization;
    }
  }

  return 0;
}

//! Implicitly set the normalization constant for the functional.
/*! The normalization constant is selected so that if an input
IceModelVec2S has entries all equal to \a scale, then the funtional value will be 1. I.e.
\f[
c_N^{-1} = \sum_{i} w_i {\tt scale}^2.
\f]*/
PetscErrorCode IPMeanSquareFunctional2S::normalize(double scale) {
  PetscErrorCode   ierr;

  // The local value of the weights
  double value = 0;

  if (m_weights) {
    IceModelVec::AccessList list(*m_weights);
    for (Points p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      value += (*m_weights)(i, j);
    }
  } else {
    for (Points p(m_grid); p; p.next()) {
      value += 1;
    }
  }

  ierr = GlobalSum(m_grid.com, &value,  &m_normalization); CHKERRQ(ierr);
  m_normalization *= (scale*scale);
  return 0;
}

PetscErrorCode IPMeanSquareFunctional2S::valueAt(IceModelVec2S &x, double *OUTPUT)  {
  PetscErrorCode   ierr;

  // The value of the objective
  double value = 0;

  IceModelVec::AccessList list(x);

  if (m_weights) {
    list.add(*m_weights);
    for (Points p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double &x_ij = x(i, j);
      value += x_ij*x_ij*(*m_weights)(i, j);
    }
  } else {
    for (Points p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      double &x_ij = x(i, j);
      value += x_ij*x_ij;
    }
  }
  value /= m_normalization;

  ierr = GlobalSum(m_grid.com, &value,  OUTPUT); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IPMeanSquareFunctional2S::dot(IceModelVec2S &a, IceModelVec2S &b, double *OUTPUT)  {
  PetscErrorCode   ierr;

  // The value of the objective
  double value = 0;

  IceModelVec::AccessList list;
  list.add(a);
  list.add(b);

  if (m_weights) {
    list.add(*m_weights);
    for (Points p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      value += (a(i, j)*b(i, j))*(*m_weights)(i, j);
    }
  } else {
    for (Points p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      value += (a(i, j)*b(i, j));
    }
  }
  value /= m_normalization;

  ierr = GlobalSum(m_grid.com, &value,  OUTPUT); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IPMeanSquareFunctional2S::gradientAt(IceModelVec2S &x, IceModelVec2S &gradient)  {
  gradient.set(0);

  IceModelVec::AccessList list;
  list.add(x);
  list.add(gradient);

  if (m_weights) {
    list.add(*m_weights);
    for (Points p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      gradient(i, j) = 2*x(i, j)*(*m_weights)(i, j) / m_normalization;
    }
  } else {
    for (Points p(m_grid); p; p.next()) {
      const int i = p.i(), j = p.j();

      gradient(i, j) = 2*x(i, j) / m_normalization;
    }
  }


  return 0;
}

} // end of namespace pism
