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

#ifndef _SSB_MODIFIER_H_
#define _SSB_MODIFIER_H_

#include "iceModelVec.hh"

//! Shallow stress balance modifier (such as the non-sliding SIA).
class SSB_Modifier
{
public:
  SSB_Modifier(IceGrid &g, const NCConfigVariable &config)
  { D_max = 0.0; }
  virtual ~SSB_Modifier() {}

  virtual PetscErrorCode init(PISMvars &vars);

  virtual PetscErrorCode update(IceModelVec2V &vel_input,
                                IceModelVec2S &D2_input,
                                bool fast) = 0;

  //! \brief Get the diffusive (SIA) vertically-averaged flux on the staggered grid.
  virtual PetscErrorCode get_diffusive_flux(IceModelVec2Stag* &result)
  { result = &diffusive_flux; return 0; }

  //! \brief Get the max diffusivity (for the adaptive time-stepping).
  virtual PetscErrorCode get_max_diffusivity(PetscReal &D)
  { D = D_max; return 0; }

  virtual PetscErrorCode get_horizontal_3d_velocity(IceModelVec3* &u_result, IceModelVec3* &v_result);
  { u_result = &u; v_result = &v; return 0; }

  virtual PetscErrorCode get_volumetric_strain_heating(IceModelVec3* &result)
  { result = &Sigma; return 0; }
protected:
  virtual PetscErrorCode compute_sigma(IceModelVec3 &Sigma) = 0;

  IceGrid &grid;
  const NCConfigVariable &config;
  PetscReal D_max;

  IceModelVec2Stag diffusive_flux;
  IceModelVec3 u, v, Sigma;
};


//! The trivial Shallow Stress Balance modifier.
class SSB_Trivial : public SSB_Modifier
{
public:
  SSB_Trivial(IceGrid &g, const NCConfigVariable &config);
  virtual ~SSB_Trivial();
  virtual PetscErrorCode update(IceModelVec2V *vel_input, bool fast);
protected:
  virtual PetscErrorCode compute_sigma(IceModelVec3 &Sigma);
};
#endif /* _SSB_MODIFIER_H_ */
