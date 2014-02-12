// Copyright (C) 2009--2014 Jed Brown, Ed Bueler and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
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
#include <map>
#include <string>

#define ICE_ISOTHERMAL_GLEN  "isothermal_glen" /* Plain isothermal Glen */
#define ICE_PB      "pb"            /* Paterson-Budd (ThermoGlenIce) */
#define ICE_GPBLD   "gpbld"         /* Paterson-Budd-Lliboutry-Duval (PolyThermalGPBLDIce) */
#define ICE_HOOKE   "hooke"         /* Hooke (ThermoGlenIceHooke) */
#define ICE_ARR     "arr"           /* Temperature dependent Arrhenius (either warm or cold) */
#define ICE_GOLDSBY_KOHLSTEDT "gk"  /* Goldsby-Kohlstedt for SIA */
#define ICE_ARRWARM "arrwarm"       /* Temperature dependent Arrhenius (should be refactored into ICE_ARR) */

typedef PetscErrorCode(*IceFlowLawCreator)(MPI_Comm, const char[],
                                           const PISMConfig &, EnthalpyConverter*, IceFlowLaw **);

class IceFlowLawFactory {
public:
  IceFlowLawFactory(MPI_Comm, const char prefix[],
                    const PISMConfig &conf,
                    EnthalpyConverter *my_EC);
  ~IceFlowLawFactory();
  PetscErrorCode setType(std::string name);
  PetscErrorCode setFromOptions();
  PetscErrorCode registerType(std::string name, IceFlowLawCreator);
  PetscErrorCode removeType(std::string name);
  PetscErrorCode create(IceFlowLaw **);
private:
  PetscErrorCode registerAll();
private:
  MPI_Comm com;
  char prefix[256];
  std::string type_name;
  std::map<std::string, IceFlowLawCreator> flow_laws;
  const PISMConfig &config;
  EnthalpyConverter *EC;
};


#endif  // __flowlaw_factory
