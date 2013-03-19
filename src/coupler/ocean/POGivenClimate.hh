// Copyright (C) 2011, 2013 Constantine Khroulev
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

#ifndef _PODIRECTFORCING_H_
#define _PODIRECTFORCING_H_

#include "PGivenClimate.hh"
#include "POModifier.hh"

class POGiven : public PGivenClimate<POModifier,PISMOceanModel>
{
public:
  POGiven(IceGrid &g, const NCConfigVariable &conf);
  virtual ~POGiven();

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode update(PetscReal my_t, PetscReal my_dt);

  virtual PetscErrorCode sea_level_elevation(PetscReal &result);

  virtual PetscErrorCode shelf_base_temperature(IceModelVec2S &result);
  virtual PetscErrorCode shelf_base_mass_flux(IceModelVec2S &result);
protected:
  IceModelVec2T *shelfbtemp, *shelfbmassflux;
private:
  PetscErrorCode allocate_POGiven();
};


#endif /* _PODIRECTFORCING_H_ */
