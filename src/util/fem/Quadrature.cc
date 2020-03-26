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

#include "Quadrature.hh"

namespace pism {
namespace fem {

const std::vector<QuadPoint>& Quadrature::points() const {
  return m_points;
}

const std::vector<double>& Quadrature::weights() const {
  return m_weights;
}

//! Build quadrature points and weights for a tensor product quadrature based on a 1D quadrature
//! rule. Uses the same 1D quadrature in both directions.
/**
   @param[in] n 1D quadrature size (the resulting quadrature has size n*n)
   @param[in] points1 1D quadrature points
   @param[in] weights1 1D quadrature weights
   @param[out] points resulting 2D quadrature points
   @param[out] weights resulting 2D quadrature weights
 */
static void tensor_product_quadrature(unsigned int n,
                                      const double *points1,
                                      const double *weights1,
                                      std::vector<QuadPoint>& points,
                                      std::vector<double>& weights) {
  unsigned int q = 0;
  for (unsigned int j = 0; j < n; ++j) {
    for (unsigned int i = 0; i < n; ++i) {
      points[q].xi = points1[i];
      points[q].eta = points1[j];

      weights[q] = weights1[i] * weights1[j];

      ++q;
    }
  }
}
Gaussian2::Gaussian2(double D) {

  // coordinates and weights of the 2-point 1D Gaussian quadrature
  double A = 1.0 / std::sqrt(3.0);

  m_points  = {{-A, 0.0, 0.0}, {A, 0.0, 0.0}};
  m_weights = {0.5 * D, 0.5 * D};
}

//! Two-by-two Gaussian quadrature on a rectangle.
Q1Quadrature4::Q1Quadrature4() {

  // coordinates and weights of the 2-point 1D Gaussian quadrature
  const double
    A           = 1.0 / sqrt(3.0),
    points2[2]  = {-A, A},
    weights2[2] = {1.0, 1.0};

  m_points.resize(4);
  m_weights.resize(4);
  tensor_product_quadrature(2, points2, weights2, m_points, m_weights);
}

Q1Quadrature9::Q1Quadrature9() {
  const double
    A         = 0.0,
    B         = sqrt(0.6),
    points3[3] = {-B, A, B};

  const double
    w1         = 5.0 / 9.0,
    w2         = 8.0 / 9.0,
    weights3[3] = {w1, w2, w1};

  m_points.resize(9);
  m_weights.resize(9);
  tensor_product_quadrature(3, points3, weights3, m_points, m_weights);
}

Q1Quadrature16::Q1Quadrature16() {
  const double
    A          = sqrt(3.0 / 7.0 - (2.0 / 7.0) * sqrt(6.0 / 5.0)), // smaller magnitude
    B          = sqrt(3.0 / 7.0 + (2.0 / 7.0) * sqrt(6.0 / 5.0)), // larger magnitude
    points4[4] = {-B, -A, A, B};

  // The weights w_i for Gaussian quadrature on the reference element with these
  // quadrature points
  const double
    w1          = (18.0 + sqrt(30.0)) / 36.0, // larger
    w2          = (18.0 - sqrt(30.0)) / 36.0, // smaller
    weights4[4] = {w2, w1, w1, w2};

  m_points.resize(16);
  m_weights.resize(16);
  tensor_product_quadrature(4, points4, weights4, m_points, m_weights);
}


//! @brief N*N-point uniform (*not* Gaussian) quadrature for integrating discontinuous
//! functions.
Q1QuadratureN::Q1QuadratureN(unsigned int N) {

  std::vector<double> xi(N), w(N);
  const double dxi = 2.0 / N;
  for (unsigned int k = 0; k < N; ++k) {
    xi[k] = -1.0 + dxi*(k + 0.5);
    w[k]  = 2.0 / N;
  }

  m_points.resize(N * N);
  m_weights.resize(N * N);
  tensor_product_quadrature(N, xi.data(), w.data(), m_points, m_weights);
}

/*!
 * 3-point Gaussian quadrature on the triangle (0,0)-(1,0)-(0,1)
 */
P1Quadrature3::P1Quadrature3() {

  const double
    one_over_six   = 1.0 / 6.0,
    two_over_three = 2.0 / 3.0;

  m_points = {{two_over_three, one_over_six, 0.0},
              {one_over_six,   two_over_three, 0.0},
              {one_over_six,   one_over_six, 0.0}};

  m_weights = {one_over_six, one_over_six, one_over_six};
}

} // end of namespace fem
} // end of namespace pism
