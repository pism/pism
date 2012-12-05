// Copyright (C) 2010, 2011, 2012 Ed Bueler, Daniella DellaGiustina, Constantine Khroulev, and Andy Aschwanden
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

#ifndef _REGIONAL_H_
#define _REGIONAL_H_

#include <petsc.h>
#include "IceGrid.hh"
#include "iceModel.hh"
#include "SIAFD.hh"
#include "SSAFD.hh"
#include "PISMMohrCoulombYieldStress.hh"
#include "PISMHydrology.hh"


//! \brief A version of the SIA stress balance with tweaks for outlet glacier
//! simulations.
class SIAFD_Regional : public SIAFD
{
public:
  SIAFD_Regional(IceGrid &g, EnthalpyConverter &e, const NCConfigVariable &c)
    : SIAFD(g, e, c) {}
  virtual ~SIAFD_Regional() {}
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode compute_surface_gradient(IceModelVec2Stag &h_x, IceModelVec2Stag &h_y);
protected:
  IceModelVec2Int *no_model_mask;
  IceModelVec2S   *usurfstore;   
};

//! \brief A version of the SSA stress balance with tweaks for outlet glacier
//! simulations.
class SSAFD_Regional : public SSAFD
{
public:
  SSAFD_Regional(IceGrid &g, IceBasalResistancePlasticLaw &b, EnthalpyConverter &e,
                 const NCConfigVariable &c)
    : SSAFD(g, b, e, c) {}
  virtual ~SSAFD_Regional() {}
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode compute_driving_stress(IceModelVec2V &taud);
protected:
  IceModelVec2Int *no_model_mask;    
  IceModelVec2S   *usurfstore, *thkstore;
};

class PISMRegionalDefaultYieldStress : public PISMMohrCoulombYieldStress
{
public:
  PISMRegionalDefaultYieldStress(IceGrid &g, const NCConfigVariable &conf, PISMHydrology *hydro)
    : PISMMohrCoulombYieldStress(g, conf, hydro) {}
  virtual ~PISMRegionalDefaultYieldStress() {}
  virtual PetscErrorCode init(PISMVars &vars);
  virtual PetscErrorCode basal_material_yield_stress(IceModelVec2S &result);
protected:
  IceModelVec2Int *no_model_mask;
};

#endif /* _REGIONAL_H_ */
