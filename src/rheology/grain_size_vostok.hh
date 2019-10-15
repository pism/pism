/* Copyright (C) 2019 PISM Authors
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

#ifndef GRAIN_SIZE_VOSTOK_H
#define GRAIN_SIZE_VOSTOK_H

#include <gsl/gsl_spline.h>

namespace pism {
namespace rheology {

//! A relationship between the age of the ice and the grain size from the Vostok core.
/*! A data set is interpolated here. The intention is that the softness of the
  ice has nontrivial dependence on its age, through its grain size, because of
  variable dustiness of the global climate. The grain size is partly determined
  by at which point in the glacial cycle the given ice fell as snow.

  The data is from [@ref DeLaChapelleEtAl98] and [@ref LipenkovEtAl89]. In
  particular, Figure A2 in the former reference was hand-sampled with an
  attempt to include the ``wiggles'' in that figure. Ages of the oldest ice (>=
  300 ka) were estimated in a necessarily ad hoc way. The age value of 10000 ka
  was added simply to give interpolation for very old ice; ages beyond that get
  constant extrapolation. Linear interpolation is done between the samples.
 */
class grain_size_vostok {
public:
  grain_size_vostok();
  ~grain_size_vostok();

  /*!
   * Return grain size in meters given ice age in years.
   */
  double operator()(double a);
private:
  static const int m_N = 22;
  static const double m_age[m_N];
  static const double m_grain_size[m_N];

  gsl_interp_accel* m_acc;
  gsl_spline*       m_spline;

  // disable copy constructor and the assignment operator:
  grain_size_vostok(const grain_size_vostok &other);
  grain_size_vostok& operator=(const grain_size_vostok&);
};

} // end of namespace rheology
} // end of namespace pism

#endif /* GRAIN_SIZE_VOSTOK_H */
