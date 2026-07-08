// Copyright (C) 2026 Andy Aschwanden
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

#include "pism/inverse/functional/IPHuberMisfit2V.hh"
#include "pism/util/Grid.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/pism_utilities.hh"

#include <cmath>

namespace pism {
namespace inverse {

//! Set the normalization constant so that a uniform residual of magnitude
//! `scale` (in the quadratic regime) gives J = 1/2. This matches the
//! IPMeanSquareFunctional2V convention up to the factor of 1/2 inherent in
//! the Huber kernel, so that `inverse.tikhonov.penalty_weight` retains its
//! usual interpretation when residuals are small.
void IPHuberMisfit2V::normalize(double scale) {

  double value = 0;

  if (m_weights != nullptr) {
    array::AccessScope list{m_weights};

    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();

      value += (*m_weights)(i, j);
    }
  } else {
    for (auto p : m_grid->points()) {
      (void) p;
      value += 1;
    }
  }

  m_normalization = GlobalSum(m_grid->com, value);
  m_normalization *= (scale * scale);
}

void IPHuberMisfit2V::valueAt(array::Vector &x, double *OUTPUT) {

  double value = 0;

  array::AccessScope list{&x};
  if (m_weights != nullptr) {
    list.add(*m_weights);
  }

  const double delta = m_delta;
  const double half_delta_sq = 0.5 * delta * delta;

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    Vector2d &x_ij = x(i, j);
    double r_sq = (x_ij.u * x_ij.u) + (x_ij.v * x_ij.v);
    double w = (m_weights != nullptr) ? (*m_weights)(i, j) : 1.0;

    double rho;
    if (r_sq <= delta * delta) {
      // Quadratic regime: rho = (1/2) r^2.
      rho = 0.5 * r_sq;
    } else {
      // Linear regime: rho = delta * (|r| - delta/2).
      rho = delta * std::sqrt(r_sq) - half_delta_sq;
    }
    value += w * rho;
  }
  value /= m_normalization;

  GlobalSum(m_grid->com, &value, OUTPUT, 1);
}

void IPHuberMisfit2V::gradientAt(array::Vector &x, array::Vector &gradient) {
  gradient.set(0);

  array::AccessScope list{&x, &gradient};
  if (m_weights != nullptr) {
    list.add(*m_weights);
  }

  const double delta = m_delta;
  const double delta_sq = delta * delta;

  for (auto p : m_grid->points()) {
    const int i = p.i(), j = p.j();

    Vector2d &x_ij = x(i, j);
    double r_sq = (x_ij.u * x_ij.u) + (x_ij.v * x_ij.v);
    double w = (m_weights != nullptr) ? (*m_weights)(i, j) : 1.0;

    // d(rho)/d(r_k) is r_k itself when |r| <= delta, and delta * r_k / |r|
    // when |r| > delta. Equivalently, multiply by min(1, delta / |r|).
    double factor;
    if (r_sq <= delta_sq) {
      factor = 1.0;
    } else {
      factor = delta / std::sqrt(r_sq);
    }

    gradient(i, j).u = w * factor * x_ij.u / m_normalization;
    gradient(i, j).v = w * factor * x_ij.v / m_normalization;
  }
}

} // end of namespace inverse
} // end of namespace pism
