/* Copyright (C) 2015, 2016, 2017, 2019, 2021, 2022 PISM Authors
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

#ifndef _SIAFD_REGIONAL_H_
#define _SIAFD_REGIONAL_H_

#include "pism/stressbalance/sia/SIAFD.hh"

namespace pism {
namespace stressbalance {

//! \brief A version of the SIA stress balance with tweaks for outlet glacier
//! simulations.
class SIAFD_Regional : public SIAFD {
public:
  SIAFD_Regional(std::shared_ptr<const Grid> g);
  virtual ~SIAFD_Regional() = default;
  void init();
private:
  void compute_surface_gradient(const Inputs &inputs,
                                array::Staggered1 &h_x,
                                array::Staggered1 &h_y);
  array::Staggered1 m_h_x_no_model;
  array::Staggered1 m_h_y_no_model;
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* _SIAFD_REGIONAL_H_ */
