/* Copyright (C) 2019, 2023 PISM Authors
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

#include <algorithm>            // std::min, std::max

#include "pism/rheology/grain_size_vostok.hh"

namespace pism {
namespace rheology {

// age in thousands of years
const double grain_size_vostok::m_age[grain_size_vostok::m_N] =
  {0.0000e+00, 5.0000e+01, 1.0000e+02, 1.2500e+02, 1.5000e+02,
   1.5800e+02, 1.6500e+02, 1.7000e+02, 1.8000e+02, 1.8800e+02,
   2.0000e+02, 2.2500e+02, 2.4500e+02, 2.6000e+02, 3.0000e+02,
   3.2000e+02, 3.5000e+02, 4.0000e+02, 5.0000e+02, 6.0000e+02,
   8.0000e+02, 1.0000e+04};

// grain size in meters
const double grain_size_vostok::m_grain_size[grain_size_vostok::m_N] =
  {1.8000e-03, 2.2000e-03, 3.0000e-03, 4.0000e-03, 4.3000e-03,
   3.0000e-03, 3.0000e-03, 4.6000e-03, 3.4000e-03, 3.3000e-03,
   5.9000e-03, 6.2000e-03, 5.4000e-03, 6.8000e-03, 3.5000e-03,
   6.0000e-03, 8.0000e-03, 8.3000e-03, 3.6000e-03, 3.8000e-03,
   9.5000e-03, 1.0000e-02};
grain_size_vostok::grain_size_vostok() {
  m_acc = gsl_interp_accel_alloc();
  m_spline = gsl_spline_alloc(gsl_interp_linear, m_N);
  gsl_spline_init(m_spline, m_age, m_grain_size, m_N);
}

grain_size_vostok::~grain_size_vostok() {
  gsl_spline_free(m_spline);
  gsl_interp_accel_free(m_acc);
}

double grain_size_vostok::operator()(double age) {
  double age_ka = age / 1000.0;

  age_ka = std::max(age_ka, m_age[0]);
  age_ka = std::min(age_ka, m_age[m_N - 1]);

  return gsl_spline_eval(m_spline, age_ka, m_acc);
}

} // end of namespace rheology
} // end of namespace pism
