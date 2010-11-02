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
  ShallowStressBalance(IceGrid &g, IceBasalResistancePlasticLaw &b, IceFlowLaw &i,
                       const NCConfigVariable &config)
  { grid = g; basal = b; ice = i; set_bc = false; max_u = max_v = 0.0; }
  virtual ~ShallowStressBalance();

  //  initialization and I/O:

  virtual PetscErrorCode init(PISMVars &vars);

  //! \brief Set the initial guess of the vertically-averaged ice velocity.
  virtual PetscErrorCode set_initial_guess(IceModelVec2V &guess)
  { return 0; }

  //! Read the initial guess from file.
  virtual PetscErrorCode read_initial_guess(string filename)
  { return 0; }

  //! \brief Save the initial guess (for restarting).
  virtual PetscErrorCode save_initial_guess(string filename)
  { return 0; }

  virtual PetscErrorCode set_boundary_conditions(IceModelVec2Mask &locations,
                                                 IceModelVec2V &velocities)
  { set_bc = true; vel_bc = &velocities; bc_locations = &locations; return 0; }

  // the "main" routine:

  virtual PetscErrorCode update(bool fast) = 0;

  // interface:

  //! \brief Get the thickness-advective (SSA) 2D velocity.
  virtual PetscErrorCode get_advective_2d_velocity(IceModelVec2V* &result)
  { result = &velocity; return 0; }

  //! \brief Get the max advective velocity (for the adaptive mass-continuity time-stepping).
  virtual PetscErrorCode get_max_2d_velocity(PetscReal &u_max, PetscReal &v_max)
  { u_max = max_u; v_max = max_v; return 0; }

  //! \brief Get the basal frictional heating (for the adaptive energy time-stepping).
  virtual PetscErrorCode get_basal_frictional_heating(IceModelVec2S* &result)
  { result = &basal_frictional_heating; return 0; }

  virtual PetscErrorCode get_D2(IceModelVec2S* &result)
  { result = &D2; return 0; }

  // helpers:

  //! \brief Extends the computational grid (vertically).
  virtual PetscErrorCode extend_the_grid(PetscInt old_Mz)
  { return 0; }
protected:
  IceGrid &grid;
  IceBasalResistancePlasticLaw &basal;
  IceFlowLaw &ice;

  IceModelVec2V velocity, *vel_bc;
  IceModelVec2Mask *bc_locations;
  IceModelVec2S basal_frictional_heating, D2;
  PetscReal max_u, max_v;
  bool set_bc;
};

class SSB_Trivial : public ShallowStressBalance
{
public:
  SSB_Trivial(IceGrid &g, IceBasalResistancePlasticLaw &b, IceFlowLaw &i,
              const NCConfigVariable &config)
    : ShallowStressBalance(g, b, i, config) {}
  virtual ~SSB_Trivial() {}
  virtual PetscErrorCode update(bool fast);
};

#endif /* _SHALLOWSTRESSBALANCE_H_ */
