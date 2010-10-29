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

//! Shallow stress balance modifier (such as the SIA).
class SSB_Modifier
{
public:
  SSB_Modifier(IceGrid &g,
               const NCConfigVariable &config);
  virtual ~SSB_Modifier();

  PetscErrorCode init();
  PetscErrorCode update(IceModelVec2V *input, bool fast);

  //! \brief Get the diffusive (SIA) vertically-averaged flux on the staggered grid.
  PetscErrorCode get_diffusive_flux(IceModelVec2Stag* &result) = 0;
  //! \brief Get the max diffusivity (for the adaptive time-stepping).
  PetscErrorCode get_max_diffusivity(PetscReal &D) = 0;
  PetscErrorCode get_horizontal_3d_velocity(IceModelVec3* &u, IceModelVec3* &v) = 0;
protected:
  IceGrid &grid;
  const NCConfigVariable &config;

  IceModelVec3 u, v;
};

//! PISM's SIA implementation.
class SIA
{
public:
  SIA(IceGrid &g,
      const NCConfigVariable &config);
  virtual ~SIA();
protected:
};

#endif /* _SSB_MODIFIER_H_ */
