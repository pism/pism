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
    ierr = m_weights->begin_access(); CHKERRQ(ierr);
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        value += (*m_weights)(i, j);
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        value += 1;
      }
    }
  }

  ierr = GlobalSum(&value, &m_normalization, m_grid.com); CHKERRQ(ierr);
  m_normalization *= (scale*scale);
  return 0;
}

PetscErrorCode IPMeanSquareFunctional2V::valueAt(IceModelVec2V &x, double *OUTPUT)  {
  PetscErrorCode   ierr;

  // The value of the objective
  double value = 0;

  ierr = x.begin_access(); CHKERRQ(ierr);

  if (m_weights) {
    ierr = m_weights->begin_access(); CHKERRQ(ierr);
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        Vector2 &x_ij = x(i, j);
        value += (x_ij.u*x_ij.u + x_ij.v*x_ij.v)*(*m_weights)(i, j);
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        Vector2 &x_ij = x(i, j);
        value += (x_ij.u*x_ij.u + x_ij.v*x_ij.v);
      }
    }
  }
  value /= m_normalization;

  ierr = GlobalSum(&value, OUTPUT, m_grid.com); CHKERRQ(ierr);

  ierr = x.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IPMeanSquareFunctional2V::dot(IceModelVec2V &a, IceModelVec2V &b, double *OUTPUT)  {
  PetscErrorCode   ierr;

  // The value of the objective
  double value = 0;

  ierr = a.begin_access(); CHKERRQ(ierr);

  ierr = b.begin_access(); CHKERRQ(ierr);

  if (m_weights) {
    ierr = m_weights->begin_access(); CHKERRQ(ierr);
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        Vector2 &a_ij = a(i, j);
        Vector2 &b_ij = b(i, j);
        value += (a_ij.u*b_ij.u + a_ij.v*b_ij.v)*(*m_weights)(i, j);
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        Vector2 &a_ij = a(i, j);
        Vector2 &b_ij = b(i, j);
        value += (a_ij.u*b_ij.u + a_ij.v*b_ij.v);
      }
    }
  }
  value /= m_normalization;

  ierr = GlobalSum(&value, OUTPUT, m_grid.com); CHKERRQ(ierr);

  ierr = a.end_access(); CHKERRQ(ierr);
  ierr = b.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IPMeanSquareFunctional2V::gradientAt(IceModelVec2V &x, IceModelVec2V &gradient)  {
  PetscErrorCode   ierr;

  gradient.set(0);

  ierr = x.begin_access(); CHKERRQ(ierr);

  ierr = gradient.begin_access(); CHKERRQ(ierr);

  if (m_weights) {
    ierr = m_weights->begin_access(); CHKERRQ(ierr);
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        gradient(i, j).u = 2*x(i, j).u*(*m_weights)(i, j) / m_normalization;
        gradient(i, j).v = 2*x(i, j).v*(*m_weights)(i, j) / m_normalization;
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        gradient(i, j).u = 2*x(i, j).u / m_normalization;
        gradient(i, j).v = 2*x(i, j).v / m_normalization;
      }
    }
  }

  ierr = x.end_access(); CHKERRQ(ierr);
  ierr = gradient.end_access(); CHKERRQ(ierr);

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
    ierr = m_weights->begin_access(); CHKERRQ(ierr);
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        value += (*m_weights)(i, j);
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        value += 1;
      }
    }
  }

  ierr = GlobalSum(&value, &m_normalization, m_grid.com); CHKERRQ(ierr);
  m_normalization *= (scale*scale);
  return 0;
}

PetscErrorCode IPMeanSquareFunctional2S::valueAt(IceModelVec2S &x, double *OUTPUT)  {
  PetscErrorCode   ierr;

  // The value of the objective
  double value = 0;

  ierr = x.begin_access(); CHKERRQ(ierr);

  if (m_weights) {
    ierr = m_weights->begin_access(); CHKERRQ(ierr);
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        double &x_ij = x(i, j);
        value += x_ij*x_ij*(*m_weights)(i, j);
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        double &x_ij = x(i, j);
        value += x_ij*x_ij;
      }
    }
  }
  value /= m_normalization;

  ierr = GlobalSum(&value, OUTPUT, m_grid.com); CHKERRQ(ierr);

  ierr = x.end_access(); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IPMeanSquareFunctional2S::dot(IceModelVec2S &a, IceModelVec2S &b, double *OUTPUT)  {
  PetscErrorCode   ierr;

  // The value of the objective
  double value = 0;

  ierr = a.begin_access(); CHKERRQ(ierr);

  ierr = b.begin_access(); CHKERRQ(ierr);

  if (m_weights) {
    ierr = m_weights->begin_access(); CHKERRQ(ierr);
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        value += (a(i, j)*b(i, j))*(*m_weights)(i, j);
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        value += (a(i, j)*b(i, j));
      }
    }
  }
  value /= m_normalization;

  ierr = GlobalSum(&value, OUTPUT, m_grid.com); CHKERRQ(ierr);

  ierr = a.end_access(); CHKERRQ(ierr);
  ierr = b.end_access(); CHKERRQ(ierr);

  return 0;
}


PetscErrorCode IPMeanSquareFunctional2S::gradientAt(IceModelVec2S &x, IceModelVec2S &gradient)  {
  PetscErrorCode ierr;

  gradient.set(0);

  ierr = x.begin_access(); CHKERRQ(ierr);

  ierr = gradient.begin_access(); CHKERRQ(ierr);

  if (m_weights) {
    ierr = m_weights->begin_access(); CHKERRQ(ierr);
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        gradient(i, j) = 2*x(i, j)*(*m_weights)(i, j) / m_normalization;
      }
    }
    ierr = m_weights->end_access(); CHKERRQ(ierr);
  } else {
    for (int i = m_grid.xs; i < m_grid.xs + m_grid.xm; i++) {
      for (int j = m_grid.ys; j < m_grid.ys + m_grid.ym; j++) {
        gradient(i, j) = 2*x(i, j) / m_normalization;
      }
    }
  }

  ierr = x.end_access(); CHKERRQ(ierr);
  ierr = gradient.end_access(); CHKERRQ(ierr);

  return 0;
}

} // end of namespace pism
