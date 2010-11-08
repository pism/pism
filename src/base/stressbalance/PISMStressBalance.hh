// Copyright (C) 2010 Constantine Khrushchev and Ed Below
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

#ifndef _PISMSTRESSBALANCE_H_
#define _PISMSTRESSBALANCE_H_

#include "iceModelVec.hh"
#include "ShallowStressBalance.hh"
#include "SSB_Modifier.hh"
#include "enthalpyConverter.hh"

//! The class defining PISM's interface to the shallow stress balance code.
class PISMStressBalance
{
public:
  PISMStressBalance(IceGrid &g, IceFlowLaw &ice, EnthalpyConverter &e,
                    const NCConfigVariable &config); // done
  virtual ~PISMStressBalance();                      // done

  //! \brief Initialize the PISMStressBalance object.
  virtual PetscErrorCode init(PISMVars &vars); // done

  //! \brief Set the initial guess of the vertically-averaged ice velocity.
  virtual PetscErrorCode set_initial_guess(IceModelVec2V &guess); // done

  //! Read the initial guess from file.
  virtual PetscErrorCode read_initial_guess(string filename); // done

  //! \brief Save the initial guess (for restarting).
  virtual PetscErrorCode save_initial_guess(string filename); // done

  //! \brief Set the vertically-averaged ice velocity boundary condition.
  /*!
   * Does not affect the SIA computation.
   */
  virtual PetscErrorCode set_boundary_conditions(IceModelVec2Mask &locations,
                                                 IceModelVec2V &velocities); // done

  virtual PetscErrorCode set_basal_melt_rate(IceModelVec2S &bmr); // done

  //! \brief Update all the fields if fast == false, only update diffusive flux
  //! and max. diffusivity otherwise.
  virtual PetscErrorCode update(bool fast); // done

  //! \brief Get the thickness-advective (SSA) 2D velocity.
  virtual PetscErrorCode get_advective_2d_velocity(IceModelVec2V* &result); // done

  //! \brief Get the diffusive (SIA) vertically-averaged flux on the staggered grid.
  virtual PetscErrorCode get_diffusive_flux(IceModelVec2Stag* &result); // done

  //! \brief Get the max diffusivity (for the adaptive time-stepping).
  virtual PetscErrorCode get_max_diffusivity(PetscReal &D); // done

  //! \brief Get the max advective velocity (for the adaptive time-stepping).
  virtual PetscErrorCode get_max_2d_velocity(PetscReal &u, PetscReal &v); // done

  // for the energy/age time step:

  //! \brief Get the 3D velocity (for the energy/age time-stepping).
  virtual PetscErrorCode get_3d_velocity(IceModelVec3* &u, IceModelVec3* &v, IceModelVec3* &w); // done
  //! \brief Get the max 3D velocity (for the adaptive time-stepping).
  virtual PetscErrorCode get_max_3d_velocity(PetscReal &u, PetscReal &v, PetscReal &w); // done
  //! \brief Get the basal frictional heating (for the energy time-stepping).
  virtual PetscErrorCode get_basal_frictional_heating(IceModelVec2S* &result); // done

  //! \brief Extends the computational grid (vertically).
  virtual PetscErrorCode extend_the_grid(PetscInt old_Mz); // done
protected:
  virtual PetscErrorCode compute_vertical_velocity(IceModelVec3 *u, IceModelVec3 *v,
                                                   IceModelVec2S *bmr, IceModelVec3 &result); // done
  IceGrid &grid;
  IceFlowLaw &ice;
  EnthalpyConverter &EC;
  const NCConfigVariable &config;

  IceModelVec3 w;
  PetscReal w_max;
  IceModelVec2S *basal_melt_rate;

  ShallowStressBalance *stress_balance;
  SSB_Modifier *modifier;
};

#endif /* _PISMSTRESSBALANCE_H_ */
