/* Copyright (C) 2014, 2015, 2021, 2022 PISM Authors
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

#include "ColumnInterpolation.hh"

#include <petscsys.h>
#include <cmath>

namespace pism {

ColumnInterpolation::ColumnInterpolation(const std::vector<double> &new_z_coarse,
                                         const std::vector<double> &new_z_fine)
  : m_z_fine(new_z_fine),
    m_z_coarse(new_z_coarse),
    m_use_linear_interpolation(false),
    m_identical_grids(false) {
  init_interpolation();
}

std::vector<double> ColumnInterpolation::coarse_to_fine(const std::vector<double> &input,
                                                        unsigned int k_max_result) const {
  std::vector<double> result(Mz_fine());
  coarse_to_fine(&input[0], k_max_result, &result[0]);
  return result;
}

void ColumnInterpolation::coarse_to_fine(const double *input, unsigned int k_max_result, double *result) const {
  if (m_identical_grids) {
#if PETSC_VERSION_LT(3, 12, 0)
    PetscMemmove(result, const_cast<double*>(input), Mz_fine()*sizeof(double));
#else
    PetscArraymove(result, input, Mz_fine());
#endif
    return;
  }

  if (m_use_linear_interpolation) {
    coarse_to_fine_linear(input, k_max_result, result);
    return;
  }

  coarse_to_fine_quadratic(input, k_max_result, result);
}

void ColumnInterpolation::coarse_to_fine_linear(const double *input, unsigned int k_max_result,
                                                double *result) const {
  const unsigned int Mzfine = Mz_fine();
  const unsigned int Mzcoarse = Mz_coarse();

  for (unsigned int k = 0; k < Mzfine; ++k) {
    if (k > k_max_result) {
      result[k] = input[m_coarse2fine[k]];
      continue;
    }

    unsigned int m = m_coarse2fine[k];

    // extrapolate (if necessary):
    if (m == Mzcoarse - 1) {
      result[k] = input[Mzcoarse - 1];
      continue;
    }

    const double incr = (m_z_fine[k] - m_z_coarse[m]) / (m_z_coarse[m + 1] - m_z_coarse[m]);
    result[k] = input[m] + incr * (input[m + 1] - input[m]);
  }
}

void ColumnInterpolation::coarse_to_fine_quadratic(const double *input, unsigned int k_max_result,
                                                   double *result) const {
  unsigned int k = 0, m = 0;
  const unsigned int Mz = Mz_coarse();
  for (m = 0; m < Mz - 2 and k <= k_max_result; ++m) {

    const double
      z0      = m_z_coarse[m],
      z1      = m_z_coarse[m + 1],
      dz_inv  = m_constants[3 * m + 0], // = 1.0 / (z1 - z0)
      dz1_inv = m_constants[3 * m + 1], // = 1.0 / (z2 - z0)
      dz2_inv = m_constants[3 * m + 2], // = 1.0 / (z2 - z1)
      f0      = input[m],
      f1      = input[m + 1],
      f2      = input[m + 2];

    const double
      d1 = (f1 - f0) * dz_inv,
      d2 = (f2 - f0) * dz1_inv,
      b  = (d2 - d1) * dz2_inv,
      a  = d1 - b * (z1 - z0),
      c  = f0;

    for (; m_z_fine[k] < z1 and k <= k_max_result; ++k) {
      const double s = m_z_fine[k] - z0;

      result[k] = s * (a + b * s) + c;
    }
  } // m-loop

  // check if we got to the end of the m-loop and use linear
  // interpolation between the remaining 2 coarse levels
  if (m == Mz - 2) {
    const double
      z0 = m_z_coarse[m],
      z1 = m_z_coarse[m + 1],
      f0 = input[m],
      f1 = input[m + 1],
      lambda = (f1 - f0) / (z1 - z0);

    for (; m_z_fine[k] < z1; ++k) {
      result[k] = f0 + lambda * (m_z_fine[k] - z0);
    }
  }

  // fill the rest using constant extrapolation
  const double f0 = input[Mz - 1];
  for (; k <= k_max_result; ++k) {
    result[k] = f0;
  }
}

std::vector<double> ColumnInterpolation::fine_to_coarse(const std::vector<double> &input) const {
  std::vector<double> result(Mz_coarse());
  fine_to_coarse(&input[0], &result[0]);
  return result;
}

void ColumnInterpolation::fine_to_coarse(const double *input, double *result) const {
  if (m_identical_grids) {
#if PETSC_VERSION_LT(3, 12, 0)
    PetscMemmove(result, const_cast<double*>(input), Mz_fine()*sizeof(double));
#else
    PetscArraymove(result, input, Mz_fine());
#endif
    return;
  }

  const unsigned int N = Mz_coarse();

  for (unsigned int k = 0; k < N - 1; ++k) {
    const int m = m_fine2coarse[k];

    const double increment = (m_z_coarse[k] - m_z_fine[m]) / (m_z_fine[m + 1] - m_z_fine[m]);
    result[k] = input[m] + increment * (input[m + 1] - input[m]);
  }

  result[N - 1] = input[m_fine2coarse[N - 1]];
}

unsigned int ColumnInterpolation::Mz_coarse() const {
  return m_z_coarse.size();
}

unsigned int ColumnInterpolation::Mz_fine() const {
  return m_z_fine.size();
}

double ColumnInterpolation::dz_fine() const {
  return m_z_fine[1] - m_z_fine[0];
}

const std::vector<double>& ColumnInterpolation::z_fine() const {
  return m_z_fine;
}

const std::vector<double>& ColumnInterpolation::z_coarse() const {
  return m_z_coarse;
}

/*!
 * Given two 1D grids, `z_input` and `z_output`, for each `z_output`
 * index `k` we find an index `m` so that
 *
 * `z_input[m] < z_output[k] <= z_input[m+1]`
 *
 * In other words, we look for two consecutive points in the input
 * grid that bracket a point in the output grid.
 *
 * This function sets `result[k] = m`. This information is then used
 * to interpolate from the grid defined by `z_input` to the one
 * defined by `z_output`.
 *
 * We use constant extrapolation outside the range defined by `z_input`.
 */
static std::vector<unsigned int> init_interpolation_indexes(const std::vector<double>& z_input,
                                                            const std::vector<double>& z_output) {
  std::vector<unsigned int> result(z_output.size());

  unsigned int m = 0;
  for (unsigned int k = 0; k < z_output.size(); ++k) {

    if (z_output[k] <= z_input.front()) {
      result[k] = 0;
      continue;
    }

    if (z_output[k] >= z_input.back()) {
      result[k] = z_input.size() - 1;
      continue;
    }

    while (z_input[m + 1] < z_output[k]) {
      ++m;
    }

    result[k] = m;
  }

  return result;
}

void ColumnInterpolation::init_interpolation() {

  if (m_z_coarse.size() == m_z_fine.size()) {

    PetscBool identical = PETSC_FALSE;
#if PETSC_VERSION_LT(3, 12, 0)
    size_t bytes = m_z_fine.size() * sizeof(double);
    PetscMemcmp(m_z_coarse.data(), m_z_fine.data(), bytes, &identical);
#else
    PetscArraycmp(m_z_coarse.data(), m_z_fine.data(), m_z_fine.size(), &identical);
#endif

    if (identical == PETSC_TRUE) {
      m_identical_grids = static_cast<bool>(identical);
      return;
    }
  }

  // coarse -> fine
  m_coarse2fine = init_interpolation_indexes(m_z_coarse, m_z_fine);

  // fine -> coarse
  m_fine2coarse = init_interpolation_indexes(m_z_fine, m_z_coarse);

  // decide if we're going to use linear or quadratic interpolation
  double dz_min = m_z_coarse.back();
  double dz_max = 0.0;
  for (unsigned int k = 0; k < Mz_coarse() - 1; ++k) {
    const double dz = m_z_coarse[k + 1] - m_z_coarse[k];
    dz_min = std::min(dz, dz_min);
    dz_max = std::max(dz, dz_max);
  }

  const double eps = 1.0e-8;
  m_use_linear_interpolation = (fabs(dz_max - dz_min) <= eps);

  // initialize quadratic interpolation constants
  if (not m_use_linear_interpolation) {
    const unsigned int N = Mz_coarse() - 2;
    m_constants.resize(3 * N);
    for (unsigned int m = 0; m < N; ++m) {
      const double
        z0 = m_z_coarse[m],
        z1 = m_z_coarse[m + 1],
        z2 = m_z_coarse[m + 2];
      m_constants[3 * m + 0] = 1.0 / (z1 - z0);
      m_constants[3 * m + 1] = 1.0 / (z2 - z0);
      m_constants[3 * m + 2] = 1.0 / (z2 - z1);
    }
  }

}

} // end of namespace pism
