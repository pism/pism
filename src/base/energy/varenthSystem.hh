// Copyright (C) 2011 Ed Bueler
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

#ifndef __varenthSystem_hh
#define __varenthSystem_hh

#include "enthSystem.hh"

//! Replacement column solver for enthalpy method, using variable conductivity for cold ice.
/*!
Like base class enthSystemCtx.  Solves a tridiagonal linear system for conservation
of energy in vertical column, for the enthalpy method [\ref AschwandenBuelerKhroulevBlatter].
Allows enthalpy-dependent conductivity and/or heat capacity in cold
ice.  Everything is the same except that the assemble_R() method is
based on formulas (4.37) and (4.39) in [\ref GreveBlatter2009].
 */
class varenthSystemCtx : public enthSystemCtx {

public:
  varenthSystemCtx(const NCConfigVariable &config, IceModelVec3 &my_Enth3,
                    int my_Mz, string my_prefix, EnthalpyConverter *EC);
  virtual ~varenthSystemCtx() {}

  PetscErrorCode initThisColumn(bool my_ismarginal,
                                PetscScalar my_lambda,
                                PetscReal ice_thickness);  
protected:
  PetscScalar k_from_T(PetscScalar T);
  virtual PetscErrorCode assemble_R();
  EnthalpyConverter *EC;  // conductivity has known dependence on T, not enthalpy
  PetscReal ice_thickness;
  bool k_depends_on_T;
};

#endif   //  ifndef __varenthSystem_hh

