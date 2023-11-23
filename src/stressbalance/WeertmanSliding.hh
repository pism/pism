/* Copyright (C) 2018, 2021, 2023 PISM Authors
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

#ifndef WEERTMANSLIDING_H
#define WEERTMANSLIDING_H

#include "pism/stressbalance/ShallowStressBalance.hh"

namespace pism {
namespace stressbalance {

class WeertmanSliding : public ShallowStressBalance {
public:
  WeertmanSliding(std::shared_ptr<const Grid> g);
  virtual ~WeertmanSliding() = default;
  virtual void update(const Inputs &inputs, bool full_update);
protected:
  void init_impl();
};

} // end of namespace stressbalance
} // end of namespace pism

#endif /* WEERTMANSLIDING_H */
