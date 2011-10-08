// Copyright (C) 2009--2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#ifndef __flowlaw_factory
#define __flowlaw_factory

#include "flowlaws.hh"

#define ICE_CUSTOM  "custom"        /* Plain isothermal Glen with customizable parameters */
#define ICE_PB      "pb"            /* Paterson-Budd (ThermoGlenIce) */
#define ICE_GPBLD   "gpbld"         /* Paterson-Budd-Lliboutry-Duval (PolyThermalGPBLDIce) */
#define ICE_HOOKE   "hooke"         /* Hooke (ThermoGlenIceHooke) */
#define ICE_ARR     "arr"           /* Temperature dependent Arrhenius (either warm or cold) */
#define ICE_HYBRID  "hybrid"        /* Goldsby-Kohlstedt for SIA, PB for SSA */
#define ICE_ARRWARM "arrwarm"       /* Temperature dependent Arrhenius (should be refactored into ICE_ARR) */

class IceFlowLawFactory {
public:
  IceFlowLawFactory(MPI_Comm, const char prefix[], const NCConfigVariable &conf,
                    EnthalpyConverter *my_EC);
  ~IceFlowLawFactory();
  PetscErrorCode setType(const char[]);
  PetscErrorCode setFromOptions();
  PetscErrorCode registerType(const char[],
		  PetscErrorCode(*)(MPI_Comm,const char[], const NCConfigVariable &,
                                    EnthalpyConverter*, IceFlowLaw **));
  PetscErrorCode create(IceFlowLaw **);
private:
  PetscErrorCode registerAll();
private:
  MPI_Comm comm;
  char prefix[256], type_name[256];
  PetscFList type_list;
  const NCConfigVariable &config;
  EnthalpyConverter *EC;
};


#endif  // __flowlaw_factory
