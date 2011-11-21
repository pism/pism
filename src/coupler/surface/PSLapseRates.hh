// Copyright (C) 2011 PISM Authors
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

#ifndef _PSLAPSERATES_H_
#define _PSLAPSERATES_H_

#include "PLapseRates.hh"
#include "PISMSurface.hh"
#include "PSModifier.hh"

class PSLapseRates : public PLapseRates<PISMSurfaceModel,PSModifier>
{
public:
  PSLapseRates(IceGrid &g, const NCConfigVariable &conf, PISMSurfaceModel* in)
    : PLapseRates<PISMSurfaceModel,PSModifier>(g, conf, in)
  {
    smb_lapse_rate = 0;
    option_prefix = "-surface";
  }

  virtual ~PSLapseRates() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode ice_surface_mass_flux(IceModelVec2S &result);
  virtual PetscErrorCode ice_surface_temperature(IceModelVec2S &result);
protected:
  PetscReal smb_lapse_rate;
};

#endif /* _PSLAPSERATES_H_ */
