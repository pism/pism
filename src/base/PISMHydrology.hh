// Copyright (C) 2012 PISM Authors
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

#ifndef _PISMHYDROLOGY_H_
#define _PISMHYDROLOGY_H_

#include "PISMComponent.hh"
class IceModelVec2S;

//! \brief The PISM subglacial hydrology model interface.
// FIXME: perhaps this will be a virtual base class
// FIXME: it is timestepping (PISMComponent_ts?) but will not use PISM's main timestepping
class PISMHydrology : public PISMComponent_Diag
{
public:
  PISMHydrology(IceGrid &g, const NCConfigVariable &conf)
    : PISMComponent_Diag(g, conf) {}
  virtual ~PISMHydrology() {}

/* FIXME:
  A PISM component needs to implement the following I/O methods:
  \li add_vars_to_output(), which adds variable names to the list of fields that need
  to be written.
  \li define_variables(), which defines variables to be written and writes variable metadata.
  \li write_variables(), which writes data itself.
*/

  virtual PetscErrorCode init_steady(IceModelVec2S W0);

  // FIXME: caution.  this updates at each call and thus should not be called repeatedly
  virtual PetscErrorCode update_water_and_pressure(PetscScalar dt);

  virtual PetscErrorCode get_water_thickness(IceModelVec2S &result);

  virtual PetscErrorCode get_water_pressure(IceModelVec2S &result);

protected:
  IceModelVec2S W, P;  // this model's state
  IceModelVec2S *bed, *thickness, *surface; // pointers into IceModel; fields describe ice
  IceModelVec2V *Ubase;  // ice sliding velocity in IceModel
  IceModelVec2S Po, cbase, alph, beta;  // sliding speed, overburden, components of water velocity
};

#endif /* _PISMHYDROLOGY_H_ */

