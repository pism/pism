/* Copyright (C) 2015, 2016, 2017, 2019, 2020, 2021, 2022 PISM Authors
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

#ifndef _SSAFD_REGIONAL_H_
#define _SSAFD_REGIONAL_H_

#include "pism/stressbalance/ssa/SSAFD.hh"

namespace pism {

namespace stressbalance {

//! \brief A version of the SSA stress balance with tweaks for outlet glacier
//! simulations.
class SSAFD_Regional : public SSAFD {
public:
  SSAFD_Regional(std::shared_ptr<const IceGrid> g);
  virtual ~SSAFD_Regional() = default;
  virtual void init();
  virtual void compute_driving_stress(const array::Scalar &ice_thickness,
                                      const array::Scalar1 &surface_elevation,
                                      const array::CellType1 &cell_type,
                                      const array::Scalar1 *no_model_mask,
                                      array::Vector &result) const;

private:
  void update(const Inputs &inputs, bool full_update);

  const array::Scalar1   *m_h_stored;
  const array::Scalar   *m_H_stored;
  const array::Scalar *m_no_model_mask;
};

SSA * SSAFD_RegionalFactory(std::shared_ptr<const IceGrid> grid);

} // end of namespace stressbalance

} // end of namespace pism

#endif /* _SSAFD_REGIONAL_H_ */
