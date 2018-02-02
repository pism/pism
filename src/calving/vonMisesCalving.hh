/* Copyright (C) 2016, 2017, 2018 PISM Authors
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

#ifndef VONMISESCALVING_H
#define VONMISESCALVING_H

#include "pism/util/iceModelVec.hh"
#include "StressCalving.hh"

namespace pism {

class IceModelVec2CellType;

namespace calving {

class vonMisesCalving : public StressCalving {
public:
  vonMisesCalving(IceGrid::ConstPtr g, stressbalance::StressBalance *stress_balance);
  virtual ~vonMisesCalving();

  void init();

  // empty methods that we're required to implement:
protected:
  virtual DiagnosticList diagnostics_impl() const;

  void compute_calving_rate(const IceModelVec2CellType &mask,
                            IceModelVec2S &result) const;

protected:
  double m_sigma_max;
};

} // end of namespace calving
} // end of namespace pism

#endif /* VONMISESCALVING_H */
