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

#ifndef _PODSLFORCING_H_
#define _PODSLFORCING_H_

#include "PScalarForcing.hh"
#include "POModifier.hh"

class POdSLforcing : public PScalarForcing<PISMOceanModel,POModifier>
{
public:
  POdSLforcing(IceGrid &g, const NCConfigVariable &conf, PISMOceanModel* in);
  virtual ~POdSLforcing() {}

  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode sea_level_elevation(PetscReal &result);
};

#endif /* _PODSLFORCING_H_ */
