/* Copyright (C) 2020, 2021, 2022 PISM Authors
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

#ifndef PISM_BLATTER_MOD_H
#define PISM_BLATTER_MOD_H

#include <memory>               // std::shared_ptr

#include "Blatter.hh"
#include "pism/stressbalance/SSB_Modifier.hh"

namespace pism {
namespace stressbalance {

/*!
 * The "modifier" post-processing the velocity field computed using `Blatter`.
 */
class BlatterMod : public SSB_Modifier {
public:
  BlatterMod(std::shared_ptr<Blatter> solver);
  virtual ~BlatterMod() = default;

  void init();

  void update(const array::Vector &sliding_velocity,
              const Inputs &inputs,
              bool full_update);
private:
  std::shared_ptr<Blatter> m_solver;

  void transfer(const array::Scalar &ice_thickness);

  void compute_max_diffusivity(const array::Vector &velocity,
                               const array::Scalar &ice_thickness,
                               const array::Scalar1 &surface);
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* PISM_BLATTER_MOD_H */
