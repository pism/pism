// Copyright (C) 2010, 2011, 2012 Constantine Khroulev
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

#ifndef _PISMSTRESSBALANCE_DIAGNOSTICS_H_
#define _PISMSTRESSBALANCE_DIAGNOSTICS_H_

#include "PISMStressBalance.hh"
#include "PISMDiagnostic.hh"


//! \brief Computes the vertically-averaged ice velocity.
class PSB_velbar : public PISMDiag<PISMStressBalance>
{
public:
  PSB_velbar(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes cbar, the magnitude of vertically-integrated horizontal
//! velocity of ice and masks out ice-free areas.
class PSB_cbar : public PISMDiag<PISMStressBalance>
{
public:
  PSB_cbar(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes cflx, the magnitude of vertically-integrated horizontal
//! flux of ice.
class PSB_cflx : public PISMDiag<PISMStressBalance>
{
public:
  PSB_cflx(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes cbase, the magnitude of horizontal velocity of ice at base
//! of ice and masks out ice-free areas.
class PSB_cbase : public PISMDiag<PISMStressBalance>
{
public:
  PSB_cbase(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes csurf, the magnitude of horizontal ice velocity at the
//! surface.
class PSB_csurf : public PISMDiag<PISMStressBalance>
{
public:
  PSB_csurf(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes velsurf, the horizontal velocity of ice at ice surface.
class PSB_velsurf : public PISMDiag<PISMStressBalance>
{
public:
  PSB_velsurf(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  PetscErrorCode compute(IceModelVec* &result);
};

//! Computes vertical ice velocity (relative to the geoid).
/*!
  \f[
  w(s) = \tilde w(s) + \frac{\partial b}{\partial t} + U(s) \cdot \nabla b
  \f]
 */
class PSB_wvel : public PISMDiag<PISMStressBalance>
{
public:
  PSB_wvel(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  PetscErrorCode compute(IceModelVec* &result);
};

//! Computes wvelsurf, the vertical velocity of ice at ice surface.
class PSB_wvelsurf : public PISMDiag<PISMStressBalance>
{
public:
  PSB_wvelsurf(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  PetscErrorCode compute(IceModelVec* &result);
};

//! Computes wvelbase, the vertical velocity of ice at the base of ice.
class PSB_wvelbase : public PISMDiag<PISMStressBalance>
{
public:
  PSB_wvelbase(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes horizontal ice velocity at the base of ice.
class PSB_velbase : public PISMDiag<PISMStressBalance>
{
public:
  PSB_velbase(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes basal frictional heating.
class PSB_bfrict : public PISMDiag<PISMStressBalance>
{
public:
  PSB_bfrict(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the x-component of the horizontal ice velocity.
class PSB_uvel : public PISMDiag<PISMStressBalance>
{
public:
  PSB_uvel(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the y-component of the horizontal ice velocity.
class PSB_vvel : public PISMDiag<PISMStressBalance>
{
public:
  PSB_vvel(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes vertical velocity of ice, relative to the bed directly
//! below.
class PSB_wvel_rel : public PISMDiag<PISMStressBalance>
{
public:
  PSB_wvel_rel(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Computes the driving shear stress at the base of ice
//! (diagnostically).
class PSB_taud_mag : public PISMDiag<PISMStressBalance>
{
public:
  PSB_taud_mag(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

//! \brief Reports the volumetric strain heating.
class PSB_strainheat : public PISMDiag<PISMStressBalance>
{
public:
  PSB_strainheat(PISMStressBalance *m, IceGrid &g, PISMVars &my_vars);
  virtual PetscErrorCode compute(IceModelVec* &result);
};

#endif /* _PISMSTRESSBALANCE_DIAGNOSTICS_H_ */
