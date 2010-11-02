// Copyright (C) 2010 Constantine Khroulev and Ed Bueler
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 2 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef _SHALLOWSTRESSBALANCE_H_
#define _SHALLOWSTRESSBALANCE_H_

#include "iceModelVec.hh"

//! Shallow stress balance (such as the SSA).
class ShallowStressBalance
{
public:
  ShallowStressBalance(IceGrid &g, IceBasalResistancePlasticLaw &b,
                       const NCConfigVariable &config)
  { grid = g; basal = b; set_bc = false; }
  virtual ~ShallowStressBalance();

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode set_boundary_conditions(IceModelVec2Mask &locations,
                                                 IceModelVec2V &velocities); // done

  virtual PetscErrorCode update(bool fast) = 0;

  //! \brief Get the thickness-advective (SSA) 2D velocity.
  virtual PetscErrorCode get_advective_2d_velocity(IceModelVec2V* &result);
  //! \brief Get the max advective velocity (for the adaptive mass-continuity time-stepping).
  virtual PetscErrorCode get_max_2d_velocity(PetscReal &u_max, PetscReal &v_max);
  //! \brief Get the basal frictional heating (for the adaptive energy time-stepping).
  virtual PetscErrorCode get_basal_frictional_heating(IceModelVec2S* &result);
protected:
  IceGrid &grid;
  IceBasalResistancePlasticLaw &basal;

  IceModelVec2V velocity, *vel_bc;
  IceModelVec2S basal_frictional_heating;
  IceModelVec2Mask bc_locations;
  PetscReal max_u, max_v;
  bool set_bc;
};

class SSB_Trivial : public ShallowStressBalance
{
public:
  SSB_Trivial(IceGrid &g, IceBasalResistancePlasticLaw &b,
              const NCConfigVariable &config)
    : ShallowStressBalance(g, basal, config) {}
  virtual ~SSB_Trivial() {}
  virtual PetscErrorCode update(bool fast)
  {
    if (fast) return 0;

    ierr = velocity.set(0.0); CHKERRQ(ierr);
    max_u = max_v = 0.0;
    ierr = basal_frictional_heating.set(0.0); CHKERRQ(ierr);

    return 0;
  }
};

#endif /* _SHALLOWSTRESSBALANCE_H_ */
