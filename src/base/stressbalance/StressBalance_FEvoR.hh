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

#ifndef _STRESSBALANCE_FEVOR_H_
#define _STRESSBALANCE_FEVOR_H_

#include "PISMStressBalance.hh"

namespace pism {

class StressBalance_FEvoR : public StressBalance {
public:
  StressBalance_FEvoR(IceGrid &g, ShallowStressBalance *sb, SSB_Modifier *ssb_mod,
                      const Config &config);
  virtual ~StressBalance_FEvoR();
protected:
  virtual PetscErrorCode compute_volumetric_strain_heating();
};

} // end of namespace pism

#endif /* _STRESSBALANCE_FEVOR_H_ */
