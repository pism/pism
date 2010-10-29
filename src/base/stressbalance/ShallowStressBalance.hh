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
  ShallowStressBalance(IceGrid &g,
                       const NCConfigVariable &config);
  virtual ~ShallowStressBalance();

  PetscErrorCode init();

  PetscErrorCode set_boundary_conditions(IceModelVec2Mask &locations,
                                         IceModelVec2V &velocities);

  PetscErrorCode update(bool fast);

  //! \brief Get the advective (SSA) 2D velocity.
  PetscErrorCode get_advective_2d_velocity(IceModelVec2V* &result);
  //! \brief Get the max advective velocity (for the adaptive time-stepping).
  PetscErrorCode get_max_2d_velocity(PetscReal &result);
  //! \brief Get the basal frictional heating (for the energy time-stepping).
  PetscErrorCode get_basal_frictional_heating(IceModelVec2S* &result);
protected:
  IceGrid &grid;
  const NCConfigVariable &config;

  IceModelVec2V velocity;
  IceModelVec2S basal_frictional_heating;
  PetscReal max_u, max_v;
};

//! PISM's SSA solver implementation
class SSAFD : public ShallowStressBalance
{
public:
  SSAFD(IceGrid &g,
        const NCConfigVariable &config);
  virtual ~SSAFD();
protected:
};

#endif /* _SHALLOWSTRESSBALANCE_H_ */
