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

#ifndef __varkenthSystem_hh
#define __varkenthSystem_hh

#include "enthSystem.hh"

//! Replacement column solver for enthalpy using variable conductivity for cold ice.
/*!
Like enthSystemCtx, just does a tridiagonal linear system for conservation
of energy in vertical column, for enthalpy.

Has additional enthalpy-dependent conductivity in cold ice. Also supports
enthalpy- (temperature) dependent heat capacity if
use_linear_in_temperature_heat_capacity is set. Everything is the same except
the assemble_R() method is based on formula (4.37) in [\ref GreveBlatter2009].

This implementation does stupid code duplication.  If we use this and think
it is worth keeping then FIXME: it should be made configurable and this code
duplication should be removed.

This is to address R. Greve's concerns about the submitted paper
[\ref AschwandenBuelerKhroulevBlatter].
 */
class varkenthSystemCtx : public enthSystemCtx {

public:
  varkenthSystemCtx(const NCConfigVariable &config, IceModelVec3 &my_Enth3,
                    int my_Mz, string my_prefix, EnthalpyConverter *EC);
  ~varkenthSystemCtx() {}

  PetscErrorCode initAllColumns(const PetscScalar my_dx, const PetscScalar my_dy, 
                                const PetscScalar my_dtTemp, const PetscScalar my_dzEQ);

  PetscErrorCode viewConstants(PetscViewer viewer, bool show_col_dependent);

protected:
  PetscScalar       getvark(PetscScalar T);
  PetscScalar       getvarc(PetscScalar T);
  PetscErrorCode    assemble_R();
  EnthalpyConverter *EC;  // needed to get temperature from enthalpy because that is
                          //   what conductivity depends on
  bool use_variable_c;
};

#endif   //  ifndef __varkenthSystem_hh

