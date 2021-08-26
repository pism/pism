/* Copyright (C) 2015, 2017, 2018, 2019, 2021 PISM Authors
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

#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#include <vector>

namespace pism {

enum InterpolationType {LINEAR, NEAREST, PIECEWISE_CONSTANT, LINEAR_PERIODIC};

/**
 * Class encapsulating linear interpolation indexes and weights.
 *
 * Most library interpolation routines (GSL's, for example) assume that we are working with a fixed
 * function given on an input grid; these libraries support evaluating this function at arbitrary
 * points. In this case an interpolation routine may use values on the input grid to create an
 * interpolant, but has no knowledge of an "output grid."
 *
 * In PISM we have a slightly different situation: we need to transfer several gridded functions
 * between two *fixed* grids, so we can use both grids to compute interpolation weights, but cannot
 * use values on the input grid.
 *
 * If a point on the output grid is outside the interval defined by the input grid this code uses
 * *constant extrapolation*.
 *
 * This class isolates the code computing interpolation weights so that we can use these 1D weights
 * in a 2D or 3D regridding code. This avoids code duplication: in bilinear (trilinear)
 * interpolation X and Y (and Z) grid dimensions are treated the same. This also makes it possible
 * to test the code computing interpolation weights in isolation.
 *
 * We provide getters `left()`, `right()` and `alpha()`.
 *
 * Here `left[i]` is the index in the input grid (`x_input`) such that
 * ~~~ c++
 * input_x[left[i]] < output_x[i]
 * ~~~
 * similar for `right[i]`:
 * ~~~ c++
 * input_x[right[i]] >= output_x[i]
 * ~~~
 * When `output_x[i]` is outside the input interval `left[i]` and `right[i]` are adjusted so that
 * the interpolation formula below corresponds to constant extrapolation.
 *
 * `alpha[i]` is the linear interpolation weight used like this:
 * ~~~ c++
 * const int
 *   L = m_left[k],
 *   R = m_right[k];
 * const double Alpha = m_alpha[k];
 * result[k] = input_values[L] + Alpha * (input_values[R] - input_values[L]);
 * ~~~
 *
 * Piecewise constant 1D interpolation used for time-dependent forcing (fluxes that should
 * be interpreted as piecewise-constant to simplify accounting of mass conservation).
 *
 * Here `input_x` defines left end-points of intervals, For example, [0, 1] defines two
 * intervals: [0, 1) and [1, x_e). Here the value x_e is irrelevant because we use
 * constant extrapolation for points outside the interval covered by input data.
 *
 *
 * Linear interpolation for periodic data (annual temperature cycles, etc).
 *
 * Input data are assumed to be periodic on the interval from 0 to `period`.
 *
 */
class Interpolation {
public:
  Interpolation(InterpolationType type, const std::vector<double> &input_x,
                const std::vector<double> &output_x, double period = 0.0);
  Interpolation(InterpolationType type, const double *input_x, unsigned int input_x_size,
                const double *output_x, unsigned int output_x_size, double period = 0.0);

  const std::vector<int>& left() const;
  const std::vector<int>& right() const;
  const std::vector<double>& alpha() const;

  int left(size_t i) const;
  int right(size_t i) const;
  double alpha(size_t i) const;

  //! Return interpolated values (on the output grid) given `input_values` on the input grid.
  /** This is used for testing. (Regular code calls left(), right(), and alpha().)
   */
  std::vector<double> interpolate(const std::vector<double> &input_values) const;

  /*!
   * Compute interpolated values on the output grid given values on the input grid.
   */
  void interpolate(const double *input, double *output) const;

  double integrate(const double *input) const;
  double integrate(const std::vector<double> &input) const;
  double interval_length() const;
private:
  //! Interpolation indexes
  std::vector<int> m_left, m_right;
  //! Interpolation weights
  std::vector<double> m_alpha;

  //! Integration weights
  std::vector<double> m_w;
  //! Length of the interval defined by `output_x`.
  double m_interval_length;

  void init_linear(const double *input_x, unsigned int input_x_size,
                   const double *output_x, unsigned int output_x_size);
  void init_nearest(const double *input_x, unsigned int input_x_size,
                    const double *output_x, unsigned int output_x_size);
  void init_piecewise_constant(const double *input_x, unsigned int input_x_size,
                               const double *output_x, unsigned int output_x_size);
  void init_linear_periodic(const double *input_x, unsigned int input_x_size,
                            const double *output_x, unsigned int output_x_size,
                            double period);

  void init_weights_linear(const double *x,
                           unsigned int x_size,
                           const double *output_x,
                           unsigned int output_x_size);

  void init_weights_piecewise_constant(const double *x,
                                       unsigned int x_size,
                                       const double *output_x,
                                       unsigned int output_x_size);
};

} // end of namespace pism

#endif /* _INTERPOLATION_H_ */
