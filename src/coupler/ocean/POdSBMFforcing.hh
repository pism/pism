// Copyright (C) 2011, 2012 PISM Authors
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

#ifndef _PODSBMFFORCING_H_
#define _PODSBMFFORCING_H_

#include "PScalarForcing.hh"
#include "PISMOcean.hh"
#include "POModifier.hh"

//! \brief Forcing using shelf base mass flux scalar time-dependent offsets.
class POdSBMFforcing : public PScalarForcing<PISMOceanModel,POModifier>
{
public:
  POdSBMFforcing(IceGrid &g, const NCConfigVariable &conf, PISMOceanModel* in);
  virtual ~POdSBMFforcing() {}

  virtual PetscErrorCode init(PISMVars &vars);

  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);
};

#endif /* _PODSBMFFORCING_H_ */
