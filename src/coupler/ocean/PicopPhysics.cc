/* Copyright (C) 2018, 2019, 2023, 2025 PISM Authors
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
#include <cmath> // sqrt
#include <cassert>              // assert
#include <algorithm>            // std::min

#include "pism/coupler/ocean/PicoppPhysics.hh"

#include "pism/util/Config.hh"

namespace pism {
namespace ocean {

PicopPhysics::PicopPhysics(const Config &config) {

}



//! equation 8 in the PICO paper.
double PicopPhysics::melt_rate(double pm_point, double Toc) const {
  // in m/s
  return m_gamma_T / (m_nu * m_lambda) * (Toc - pm_point);
}


} // end of namespace ocean
} // end of namespace pism
