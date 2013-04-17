// Copyright (C) 2009--2013 Jed Brown, Ed Bueler and Constantine Khroulev
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
#include <map>
#include <string>

// use namespace std BUT remove trivial namespace browser from doxygen-erated HTML source browser
/// @cond NAMESPACE_BROWSER
using namespace std;
/// @endcond

#define ICE_ISOTHERMAL_GLEN  "isothermal_glen" /* Plain isothermal Glen */
#define ICE_PB      "pb"            /* Paterson-Budd (ThermoGlenIce) */
#define ICE_GPBLD   "gpbld"         /* Paterson-Budd-Lliboutry-Duval (PolyThermalGPBLDIce) */
#define ICE_HOOKE   "hooke"         /* Hooke (ThermoGlenIceHooke) */
#define ICE_ARR     "arr"           /* Temperature dependent Arrhenius (either warm or cold) */
#define ICE_GOLDSBY_KOHLSTEDT "gk"  /* Goldsby-Kohlstedt for SIA */
#define ICE_ARRWARM "arrwarm"       /* Temperature dependent Arrhenius (should be refactored into ICE_ARR) */

typedef PetscErrorCode(*IceFlowLawCreator)(MPI_Comm, const char[], PISMUnitSystem,
                                           const NCConfigVariable &, EnthalpyConverter*, IceFlowLaw **);

class IceFlowLawFactory {
public:
  IceFlowLawFactory(MPI_Comm, const char prefix[], PISMUnitSystem unit_system,
                    const NCConfigVariable &conf,
                    EnthalpyConverter *my_EC);
  ~IceFlowLawFactory();
  PetscErrorCode setType(string name);
  PetscErrorCode setFromOptions();
  PetscErrorCode registerType(string name, IceFlowLawCreator);
  PetscErrorCode removeType(string name);
  PetscErrorCode create(IceFlowLaw **);
private:
  PetscErrorCode registerAll();
private:
  MPI_Comm com;
  char prefix[256];
  string type_name;
  map<string, IceFlowLawCreator> flow_laws;
  const NCConfigVariable &config;
  EnthalpyConverter *EC;
  PISMUnitSystem m_unit_system;
};


#endif  // __flowlaw_factory
