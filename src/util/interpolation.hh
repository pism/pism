/* Copyright (C) 2015, 2017, 2018 PISM Authors
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
 */
class Interpolation {
public:
  virtual ~Interpolation();

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
protected:
  std::vector<int> m_left, m_right;
  std::vector<double> m_alpha;
  void init_linear(const double *input_x, unsigned int input_x_size,
                   const double *output_x, unsigned int output_x_size);
  void init_nearest(const double *input_x, unsigned int input_x_size,
                    const double *output_x, unsigned int output_x_size);
  Interpolation();
};

class LinearInterpolation : public Interpolation {
public:
  LinearInterpolation(const std::vector<double> &input_x,
                      const std::vector<double> &output_x);
  LinearInterpolation(const double *input_x, unsigned int input_x_size,
                      const double *output_x, unsigned int output_x_size);
};

/*!
 * Nearest neighbor interpolation.
 *
 * Implemented as a modification of linear interpolation by adjusting interpolation
 * weights.
 */
class NearestNeighbor : public Interpolation {
public:
  NearestNeighbor(const std::vector<double> &input_x,
                      const std::vector<double> &output_x);
  NearestNeighbor(const double *input_x, unsigned int input_x_size,
                  const double *output_x, unsigned int output_x_size);
};

} // end of namespace pism

#endif /* _INTERPOLATION_H_ */
