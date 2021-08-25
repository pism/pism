/* Copyright (C) 2015, 2016, 2017, 2018, 2019, 2021 PISM Authors
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

#include <gsl/gsl_interp.h>
#include <cassert>

#include "interpolation.hh"
#include "error_handling.hh"

namespace pism {

Interpolation::Interpolation(InterpolationType type,
                             const std::vector<double> &input_x,
                             const std::vector<double> &output_x,
                             double period)
  : Interpolation(type, input_x.data(), input_x.size(),
                  output_x.data(), output_x.size(), period) {
  // empty
}

Interpolation::Interpolation(InterpolationType type,
                             const double *input_x, unsigned int input_x_size,
                             const double *output_x, unsigned int output_x_size,
                             double period) {
  switch (type) {
  case LINEAR:
    init_linear(input_x, input_x_size, output_x, output_x_size);
    break;
  case NEAREST:
    init_nearest(input_x, input_x_size, output_x, output_x_size);
    break;
  case PIECEWISE_CONSTANT:
    init_piecewise_constant(input_x, input_x_size, output_x, output_x_size);
    break;
  case LINEAR_PERIODIC:
    init_linear_periodic(input_x, input_x_size, output_x, output_x_size, period);
    break;
  default:
    throw RuntimeError(PISM_ERROR_LOCATION, "invalid interpolation type");
  }
}

/**
 * Compute linear interpolation indexes and weights.
 *
 * @param[in] input_x coordinates of the input grid
 * @param[in] input_x_size number of points in the input grid
 * @param[in] output_x coordinates of the output grid
 * @param[in] output_x_size number of points in the output grid
 */
void Interpolation::init_linear(const double *input_x, unsigned int input_x_size,
                                const double *output_x, unsigned int output_x_size) {

  m_left.resize(output_x_size);
  m_right.resize(output_x_size);
  m_alpha.resize(output_x_size);

  // the trivial case (the code below requires input_x_size >= 2)
  if (input_x_size < 2) {
    for (unsigned int k = 0; k < output_x_size; ++k) {
      m_left[k]  = 0.0;
      m_right[k] = 0.0;
      m_alpha[k] = 0.0;
    }
    return;
  }

  // input grid points have to be stored in the increasing order
  for (unsigned int i = 0; i < input_x_size - 1; ++i) {
    if (input_x[i] >= input_x[i + 1]) {
      throw RuntimeError(PISM_ERROR_LOCATION, "an input grid for linear interpolation has to be "
                         "strictly increasing");
    }
  }

  // compute indexes and weights
  for (unsigned int i = 0; i < output_x_size; ++i) {
    double x = output_x[i];

    unsigned int
      L = gsl_interp_bsearch(input_x, x, 0, input_x_size),
      R = L + 1;

    double alpha = 0.0;
    if (x >= input_x[L] and R < input_x_size) {
      // regular case
      alpha = (x - input_x[L]) / (input_x[R] - input_x[L]);
    } else {
      // extrapolation
      alpha = 0.0;
      R = L;
    }

    assert(L < input_x_size);
    assert(R < input_x_size);
    assert(alpha >= 0.0 and alpha <= 1.0);

    m_left[i]  = L;
    m_right[i] = R;
    m_alpha[i] = alpha;
  }

  init_weights_linear(input_x, input_x_size, output_x, output_x_size);
}

const std::vector<int>& Interpolation::left() const {
  return m_left;
}

const std::vector<int>& Interpolation::right() const {
  return m_right;
}

const std::vector<double>& Interpolation::alpha() const {
  return m_alpha;
}

int Interpolation::left(size_t j) const {
  return m_left[j];
}

int Interpolation::right(size_t j) const {
  return m_right[j];
}

double Interpolation::alpha(size_t j) const {
  return m_alpha[j];
}

std::vector<double> Interpolation::interpolate(const std::vector<double> &input_values) const {
  std::vector<double> result(m_alpha.size());

  interpolate(input_values.data(), result.data());

  return result;
}

void Interpolation::interpolate(const double *input, double *output) const {
  size_t n = m_alpha.size();
  for (size_t k = 0; k < n; ++k) {
    const int
      L = m_left[k],
      R = m_right[k];
    output[k] = input[L] + m_alpha[k] * (input[R] - input[L]);
  }
}

void Interpolation::init_nearest(const double *input_x, unsigned int input_x_size,
                                 const double *output_x, unsigned int output_x_size) {

  init_linear(input_x, input_x_size, output_x, output_x_size);

  for (unsigned int j = 0; j < m_alpha.size(); ++j) {
    m_alpha[j] = m_alpha[j] > 0.5 ? 1.0 : 0.0;
  }
}

/*!
 * Input grid `input_x` corresponds to *left* end-points of intervals.
 */
void Interpolation::init_piecewise_constant(const double *input_x, unsigned int input_x_size,
                                            const double *output_x, unsigned int output_x_size) {

  m_left.resize(output_x_size);
  m_right.resize(output_x_size);
  m_alpha.resize(output_x_size);

  // the trivial case
  if (input_x_size < 2) {
    for (unsigned int i = 0; i < output_x_size; ++i) {
      m_left[i]  = 0;
      m_right[i] = 0;
      m_alpha[i] = 0.0;
    }
    return;
  }

  // input grid points have to be stored in the increasing order
  for (unsigned int i = 0; i < input_x_size - 1; ++i) {
    if (input_x[i] >= input_x[i + 1]) {
      throw RuntimeError(PISM_ERROR_LOCATION, "an input grid for interpolation has to be "
                         "strictly increasing");
    }
  }

  // compute indexes and weights
  for (unsigned int i = 0; i < output_x_size; ++i) {

    size_t L = gsl_interp_bsearch(input_x, output_x[i], 0, input_x_size);

    m_left[i] = L;
    m_right[i] = L;
    m_alpha[i] = 0.0;

    assert(m_left[i] >= 0 and m_left[i] < (int)input_x_size);
    assert(m_right[i] >= 0 and m_right[i] < (int)input_x_size);
    assert(m_alpha[i] >= 0.0 and m_alpha[i] <= 1.0);
  }

  init_weights_piecewise_constant(input_x, input_x_size, output_x, output_x_size);
}

void Interpolation::init_linear_periodic(const double *input_x, unsigned int input_x_size,
                                         const double *output_x, unsigned int output_x_size,
                                         double period) {

  assert(period > 0);

  m_left.resize(output_x_size);
  m_right.resize(output_x_size);
  m_alpha.resize(output_x_size);

  // the trivial case
  if (input_x_size < 2) {
    for (unsigned int i = 0; i < output_x_size; ++i) {
      m_left[i]  = 0;
      m_right[i] = 0;
      m_alpha[i] = 0.0;
    }
    return;
  }

  // input grid points have to be stored in the increasing order
  for (unsigned int i = 0; i < input_x_size - 1; ++i) {
    if (input_x[i] >= input_x[i + 1]) {
      throw RuntimeError(PISM_ERROR_LOCATION, "an input grid for interpolation has to be "
                         "strictly increasing");
    }
  }

  // compute indexes and weights
  for (unsigned int i = 0; i < output_x_size; ++i) {
    double x = output_x[i];

    unsigned int L = 0, R = 0;
    if (x < input_x[0]) {
      L = input_x_size - 1;
      R = 0.0;
    } else {
      L = gsl_interp_bsearch(input_x, x, 0, input_x_size);
      R = L + 1 < input_x_size ? L + 1 : 0;
    }

    double
      x_l = input_x[L],
      x_r = input_x[R],
      alpha = 0.0;
    if (L < R) {
      // regular case
      alpha = (x - x_l) / (x_r - x_l);
    } else {
      double
        x0 = input_x[0],
        dx = (period - x_l) + x0;
      assert(dx > 0);
      if (x > x0) {
        // interval from the last point of the input grid to the period
        alpha = (x - x_l) / dx;
      } else {
        // interval from 0 to the first point of the input grid
        alpha = 1.0 - (x_r - x) / dx;
      }
    }

    assert(L < input_x_size);
    assert(R < input_x_size);
    assert(alpha >= 0.0 and alpha <= 1.0);

    m_left[i]  = L;
    m_right[i] = R;
    m_alpha[i] = alpha;
  }
}

/*!
 * Initialize integration weights (trapezoid rule).
 */
void Interpolation::init_weights_linear(const double *x,
                                        unsigned int x_size,
                                        const double *output_x,
                                        unsigned int output_x_size) {

  if (output_x_size == 1) {
    // only one point requested: this cannot define an interval to integrate over
    m_w = {0.0};
    m_interval_length = 0.0;
    return;
  }

  int
    N = output_x_size - 1,
    al = m_left[0],
    ar = m_right[0],
    bl = m_left[N],
    br = m_right[N];

  double
    a = output_x[0],
    b = output_x[N],
    alpha_a = m_alpha[0],
    alpha_b = m_alpha[N];

  if (a >= b) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION, "invalid interval: (%f, %f)", a, b);
  }

  m_interval_length = b - a;

  m_w.resize(x_size);
  for (unsigned int k = 0; k < x_size; ++k) {
    m_w[k] = 0.0;
  }

  if (al == bl and ar == br) {
    // both end points are in the same interval

    m_w[al] += 0.5 * (2.0 - alpha_a - alpha_b) * (b - a);
    m_w[ar] += 0.5 * (alpha_a + alpha_b) * (b - a);
  } else {
    // first interval
    m_w[al] += 0.5 * (1.0 - alpha_a) * (x[ar] - a);
    m_w[ar] += 0.5 * (1.0 + alpha_a) * (x[ar] - a);

    // intermediate intervals
    for (int k = ar; k < bl; ++k) {
      int
        L = k,
        R = k + 1;
      m_w[L] += 0.5 * (x[R] - x[L]);
      m_w[R] += 0.5 * (x[R] - x[L]);
    }

    // last interval
    m_w[bl] += 0.5 * (2.0 - alpha_b) * (b - x[bl]);
    m_w[br] += 0.5 * alpha_b * (b - x[bl]);
  }
}

/*!
 * Initialize integration weights in the piecewise-constant case
 */
void Interpolation::init_weights_piecewise_constant(const double *x,
                                                    unsigned int x_size,
                                                    const double *output_x,
                                                    unsigned int output_x_size) {

  if (output_x_size == 1) {
    // only one point requested: this cannot define an interval to integrate over
    m_w = {0.0};
    m_interval_length = 0.0;
    return;
  }

  int
    N = output_x_size - 1,
    al = m_left[0],
    ar = m_right[0],
    bl = m_left[N],
    br = m_right[N];

  double
    a = output_x[0],
    b = output_x[N];

  if (a >= b) {
    throw RuntimeError::formatted(PISM_ERROR_LOCATION,
                                  "invalid interval: (%f, %f)", a, b);
  }

  m_interval_length = b - a;

  m_w.resize(x_size);
  for (unsigned int k = 0; k < x_size; ++k) {
    m_w[k] = 0.0;
  }

  if (al == bl and ar == br) {
    // both end points are in the same interval
    m_w[al] += (b - a);
  } else {
    // first interval
    m_w[al] += (x[ar] - a);

    // intermediate intervals
    for (int k = ar; k < bl; ++k) {
      m_w[k] += x[k + 1] - x[k];
    }

    // last interval
    m_w[bl] += b - x[bl];
  }
}

double Interpolation::integral(const double *input) const {
  double result = 0.0;
  for (size_t k = 0; k < m_w.size(); ++k) {
    result += input[k] * m_w[k];
  }
  return result;
}

double Interpolation::integral(const std::vector<double> &input) const {
  return integral(input.data());
}

double Interpolation::interval_length() const {
  return m_interval_length;
}


} // end of namespace pism
