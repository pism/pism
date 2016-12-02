/* Copyright (C) 2016 PISM Authors
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

#ifndef STRESSCALVING_H
#define STRESSCALVING_H

#include "CalvingFrontRetreat.hh"

namespace pism {

namespace stressbalance {
class StressBalance;
}

namespace calving {

/*! @brief An abstract class containing fields used by all stress-based calving methods. */
class StressCalving : public CalvingFrontRetreat {
public:
  StressCalving(IceGrid::ConstPtr g, stressbalance::StressBalance *stress_balance,
                unsigned int stencil_width);
  virtual ~StressCalving();

protected:
  void update_strain_rates() const;

protected:
  mutable IceModelVec2 m_strain_rates;
  const int m_stencil_width;
  stressbalance::StressBalance *m_stress_balance;
};


} // end of namespace calving
} // end of namespace pism


#endif /* STRESSCALVING_H */
