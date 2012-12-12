// Copyright (C) 2004--2012 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef _PISMYIELDSTRESS_H_
#define _PISMYIELDSTRESS_H_

#include "PISMComponent.hh"
class IceModelVec2S;

//! \brief The PISM basal yield stress model interface (virtual base class)
class PISMYieldStress : public PISMComponent_TS
{
public:
  PISMYieldStress(IceGrid &g, const NCConfigVariable &conf)
    : PISMComponent_TS(g, conf) {}
  virtual ~PISMYieldStress() {}

  virtual PetscErrorCode basal_material_yield_stress(IceModelVec2S &result) = 0;
};

#endif /* _PISMYIELDSTRESS_H_ */
