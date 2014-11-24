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

namespace pism {

#define ICE_ISOTHERMAL_GLEN  "isothermal_glen" /* Plain isothermal Glen */
#define ICE_PB      "pb"            /* Paterson-Budd (ThermoGlenIce) */
#define ICE_GPBLD   "gpbld"         /* Paterson-Budd-Lliboutry-Duval (PolyThermalGPBLDIce) */
#define ICE_HOOKE   "hooke"         /* Hooke (ThermoGlenIceHooke) */
#define ICE_ARR     "arr"           /* Temperature dependent Arrhenius (either warm or cold) */
#define ICE_GOLDSBY_KOHLSTEDT "gk"  /* Goldsby-Kohlstedt for SIA */
#define ICE_ARRWARM "arrwarm"       /* Temperature dependent Arrhenius (should be refactored into ICE_ARR) */

typedef IceFlowLaw*(*IceFlowLawCreator)(MPI_Comm, const std::string &,
                                        const Config &, EnthalpyConverter*);

class IceFlowLawFactory {
public:
  IceFlowLawFactory(MPI_Comm, const std::string &prefix,
                    const Config &conf,
                    EnthalpyConverter *my_EC);
  ~IceFlowLawFactory();
  void setType(const std::string &name);
  void setFromOptions();
  void registerType(const std::string &name, IceFlowLawCreator);
  void removeType(const std::string &name);
  IceFlowLaw* create();
private:
  void registerAll();
private:
  MPI_Comm com;
  std::string type_name, prefix;
  std::map<std::string, IceFlowLawCreator> flow_laws;
  const Config &config;
  EnthalpyConverter *EC;
};


} // end of namespace pism

#endif  // __flowlaw_factory
