/* Copyright (C) 2015, 2016, 2017, 2018, 2019, 2021, 2023 PISM Authors
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

#include "pism/util/interpolation.hh"
#include "pism/util/error_handling.hh"

namespace pism {

Interpolation::Interpolation(InterpolationType type,
                             const std::vector<double> &input_x,
                             const std::vector<double> &output_x)
  : Interpolation(type, input_x.data(), input_x.size(),
                  output_x.data(), output_x.size()) {
  // empty
}

Interpolation::Interpolation(InterpolationType type,
                             const double *input_x, unsigned int input_x_size,
                             const double *output_x, unsigned int output_x_size) {

  // the trivial case (the code below requires input_x_size >= 2)
  if (input_x_size < 2) {
    m_left.resize(output_x_size);
    m_right.resize(output_x_size);
    m_alpha.resize(output_x_size);

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
      throw RuntimeError(PISM_ERROR_LOCATION,
                         "an input grid for interpolation has to be "
                         "strictly increasing");
    }
  }

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
    // LCOV_EXCL_START
  default:
    throw RuntimeError(PISM_ERROR_LOCATION, "invalid interpolation type");
    // LCOV_EXCL_STOP
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
  assert(input_x_size >= 2);

  m_left.resize(output_x_size);
  m_right.resize(output_x_size);
  m_alpha.resize(output_x_size);

  // compute indexes and weights
  for (unsigned int i = 0; i < output_x_size; ++i) {
    double x = output_x[i];

    // note: use "input_x_size" instead of "input_x_size - 1" to support extrapolation on
    // the right
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
  assert(input_x_size >= 2);

  m_left.resize(output_x_size);
  m_right.resize(output_x_size);
  m_alpha.resize(output_x_size);

  // compute indexes and weights
  for (unsigned int i = 0; i < output_x_size; ++i) {

    // note: use "input_x_size" instead of "input_x_size - 1" to support extrapolation on
    // the right
    size_t L = gsl_interp_bsearch(input_x, output_x[i], 0, input_x_size);

    m_left[i] = L;
    m_right[i] = L;
    m_alpha[i] = 0.0;

    assert(m_left[i] >= 0 and m_left[i] < (int)input_x_size);
    assert(m_right[i] >= 0 and m_right[i] < (int)input_x_size);
    assert(m_alpha[i] >= 0.0 and m_alpha[i] <= 1.0);
  }
}

static std::map<size_t, double> weights_piecewise_constant(const double *x,
                                                           size_t x_size,
                                                           double a,
                                                           double b) {

  size_t
    al = gsl_interp_bsearch(x, a, 0, x_size),
    ar = (a >= x[al] and al + 1 < x_size) ? al + 1 : al,
    bl = gsl_interp_bsearch(x, b, 0, x_size),
    br = (b >= x[bl] and bl + 1 < x_size) ? bl + 1 : bl;

  std::map<size_t, double> result;

  if (al == bl and ar == br) {
    // both end points are in the same interval
    result[al] += (b - a);
  } else {
    // first interval
    result[al] += (x[ar] - a);

    // intermediate intervals
    for (size_t k = ar; k < bl; ++k) {
      result[k] += x[k + 1] - x[k];
    }

    // last interval
    result[bl] += b - x[bl];
  }

  return result;
}

static std::map<size_t, double> weights_piecewise_linear(const double *x,
                                                         size_t x_size,
                                                         double a,
                                                         double b) {

  size_t
    al = gsl_interp_bsearch(x, a, 0, x_size),
    ar = (a >= x[al] and al + 1 < x_size) ? al + 1 : al,
    bl = gsl_interp_bsearch(x, b, 0, x_size),
    br = (b >= x[bl] and bl + 1 < x_size) ? bl + 1 : bl;

  double
    alpha_a = (al == ar) ? 0.0 : (a - x[al]) / (x[ar] - x[al]),
    alpha_b = (bl == br) ? 0.0 : (b - x[bl]) / (x[br] - x[bl]);

  std::map<size_t, double> result;

  if (al == bl and ar == br) {
    // both end points are in the same interval

    result[al] += 0.5 * (2.0 - alpha_a - alpha_b) * (b - a);
    result[ar] += 0.5 * (alpha_a + alpha_b) * (b - a);
  } else {
    // first interval
    result[al] += 0.5 * (1.0 - alpha_a) * (x[ar] - a);
    result[ar] += 0.5 * (1.0 + alpha_a) * (x[ar] - a);

    // intermediate intervals
    for (size_t k = ar; k < bl; ++k) {
      int
        L = k,
        R = k + 1;
      result[L] += 0.5 * (x[R] - x[L]);
      result[R] += 0.5 * (x[R] - x[L]);
    }

    // last interval
    result[bl] += 0.5 * (2.0 - alpha_b) * (b - x[bl]);
    result[br] += 0.5 * alpha_b * (b - x[bl]);
  }

  return result;
}

/*!
 * Compute weights for integrating a piece-wise linear or piece-wise constant function
 * defined on the grid `x` from `a` to `b`.
 *
 * Uses constant extrapolation, both on the left and on the right.
 *
 * In the piece-wise constant case points in `x` are interpreted as left end points of
 * intervals, i.e. the value `data[k]` corresponds to the interval `x[k], x[k + 1]`.
 *
 * To evaluate in the integral compute the dot product of data on the grid `x` with
 * weights computed by this function:
 *
 * ```
 * double result = 0.0;
 * for (const auto &weight : weights) {
 *   size_t k = weight.first;
 *   double w = weight.second;
 *   result += w * data[k];
 * }
 * ```
 */
std::map<size_t, double> integration_weights(const double *x,
                                             size_t x_size,
                                             InterpolationType type,
                                             double a,
                                             double b) {

  if (a >= b) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "invalid integration interval (a >= b)");
  }

  if (type != LINEAR and type != PIECEWISE_CONSTANT) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "unsupported interpolation type");
  }

  auto weights = (type == LINEAR) ?
    weights_piecewise_linear :
    weights_piecewise_constant;

  size_t N  = x_size - 1;
  double t0 = x[0];
  double t1 = x[N];

  // both points are to the left of [t0, t1]
  if (b <= t0) {
    return {{0, b - a}};
  }

  // both points are to the right of [t0, t1]
  if (a >= t1) {
    return {{N, b - a}};
  }

  // a is to the left of [t0, t1]
  if (a < t0) {
    auto W = weights(x, x_size, t0, b);
    W[0] += t0 - a;
    return W;
  }

  // b is to the right of [t0, t1]
  if (b > t1) {
    auto W = weights(x, x_size, a, t1);
    W[N] += b - t1;
    return W;
  }

  // both a and b are inside [t0, t1]
  return weights(x, x_size, a, b);
}

std::map<size_t, double> integration_weights(const std::vector<double> &x,
                                             InterpolationType type,
                                             double a,
                                             double b) {
  return integration_weights(x.data(), x.size(), type, a, b);
}

} // end of namespace pism
