/* Copyright (C) 2014 PISM Authors
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

#include <cmath>

namespace pism {

ColumnInterpolation::ColumnInterpolation(std::vector<double> coarse)
  : m_z_coarse(coarse) {
  init_fine_grid();
  init_interpolation();
}

void ColumnInterpolation::coarse_to_fine(double *coarse, unsigned int ks, double *fine) const {
  if (m_use_linear_interpolation) {
    coarse_to_fine_linear(coarse, ks, fine);
  } else {
    coarse_to_fine_quadratic(coarse, ks, fine);
  }
}

void ColumnInterpolation::coarse_to_fine_linear(double *coarse, unsigned int ks, double *fine) const {
  for (unsigned int k = 0; k < Mz_fine(); k++) {
    if (k > ks) {
      fine[k] = coarse[m_coarse2fine[k]];
      continue;
    }

    unsigned int m = m_coarse2fine[k];

    // extrapolate (if necessary):
    if (m == Mz() - 1) {
      fine[k] = coarse[Mz() - 1];
      continue;
    }

    const double incr = (m_z_fine[k] - m_z_coarse[m]) / (m_z_coarse[m + 1] - m_z_coarse[m]);
    fine[k] = coarse[m] + incr * (coarse[m + 1] - coarse[m]);
  }
}

void ColumnInterpolation::coarse_to_fine_quadratic(double *coarse, unsigned int ks, double *fine) const {
  unsigned int k = 0, m = 0;
  for (m = 0; m < Mz() - 2; m++) {
    if (k > ks) {
      break;
    }

    const double
      z0 = m_z_coarse[m],
      z1 = m_z_coarse[m+1],
      z2 = m_z_coarse[m+2],
      f0 = coarse[m],
      f1 = coarse[m+1],
      f2 = coarse[m+2];

    const double
      d1 = (f1 - f0) / (z1 - z0),
      d2 = (f2 - f0) / (z2 - z0),
      b  = (d2 - d1) / (z2 - z1),
      a  = d1 - b * (z1 - z0),
      c  = f0;

    while (m_z_fine[k] < z1) {
      if (k > ks) {
        break;
      }

      const double s = m_z_fine[k] - z0;

      fine[k] = s * (a + b * s) + c;

      k++;
    }
  } // m-loop

  // check if we got to the end of the m-loop and use linear
  // interpolation between the remaining 2 coarse levels
  if (m == Mz() - 2) {
    const double
      z0 = m_z_coarse[m],
      z1 = m_z_coarse[m+1],
      f0 = coarse[m],
      f1 = coarse[m+1],
      lambda = (f1 - f0) / (z1 - z0);

    while (m_z_fine[k] < z1) {
      fine[k] = f0 + lambda * (m_z_fine[k] - z0);

      k++;
    }
  }

  // fill the rest using constant extrapolation
  const double f0 = coarse[Mz() - 1];
  while (k <= ks) {
    fine[k] = f0;
    k++;
  }
}

void ColumnInterpolation::fine_to_coarse(double *fine, double *coarse) const {
  const unsigned int N = Mz();

  for (unsigned int k = 0; k < N - 1; ++k) {
    const int m = m_fine2coarse[k];

    const double increment = (m_z_coarse[k] - m_z_fine[m]) / (m_z_fine[m + 1] - m_z_fine[m]);
    coarse[k] = fine[m] + increment * (fine[m + 1] - fine[m]);
  }

  coarse[N - 1] = fine[m_fine2coarse[N - 1]];
}

unsigned int ColumnInterpolation::Mz() const {
  return m_z_coarse.size();
}

unsigned int ColumnInterpolation::Mz_fine() const {
  return m_z_fine.size();
}

double ColumnInterpolation::dz_fine() const {
  return m_z_fine[1] - m_z_fine[0];
}

void ColumnInterpolation::init_fine_grid() {
  // Compute dz_fine as the minimum vertical spacing in the coarse
  // grid:
  unsigned int Mz = m_z_coarse.size();
  double Lz = m_z_coarse.back(), dz_fine = Lz;
  for (unsigned int k = 1; k < Mz; ++k) {
    dz_fine = std::min(dz_fine, m_z_coarse[k] - m_z_coarse[k-1]);
  }

  size_t Mz_fine = static_cast<size_t>(ceil(Lz / dz_fine) + 1);
  dz_fine = Lz / (Mz_fine - 1);

  m_z_fine.resize(Mz_fine);
  // compute levels of the fine grid:
  for (unsigned int k = 0; k < Mz_fine; k++) {
    m_z_fine[k] = m_z_coarse[0] + k * dz_fine;
  }
  // Note that it's allowed to go over Lz.
}

void ColumnInterpolation::init_interpolation() {
  unsigned int m = 0;
  double Lz = m_z_coarse.back();

  // coarse -> fine
  m_coarse2fine.resize(Mz_fine());
  m = 0;
  for (unsigned int k = 0; k < Mz_fine(); k++) {
    if (m_z_fine[k] >= Lz) {
      m_coarse2fine[k] = Mz() - 1;
      continue;
    }

    while (m_z_coarse[m + 1] < m_z_fine[k]) {
      m++;
    }

    m_coarse2fine[k] = m;
  }

  // fine -> coarse
  m_fine2coarse.resize(Mz());
  m = 0;
  for (unsigned int k = 0; k < Mz(); k++) {
    while (m < Mz_fine() - 1 &&
           m_z_fine[m + 1] < m_z_coarse[k]) {
      m++;
    }

    m_fine2coarse[k] = m;
  }

  // decide if we're going to use linear or quadratic interpolation
  double dz_min = Lz;
  double dz_max = 0.0;
  for (unsigned int k = 0; k < Mz() - 1; k++) {
    const double dz = m_z_coarse[k+1] - m_z_coarse[k];
    dz_min = std::min(dz, dz_min);
    dz_max = std::max(dz, dz_max);
  }

  if (fabs(dz_max - dz_min) <= 1.0e-8) {
    m_use_linear_interpolation = true;
  } else {
    m_use_linear_interpolation = false;
  }

}

} // end of namespace pism
