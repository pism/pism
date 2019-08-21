// Copyright (C) 2009--2015, 2017, 2018 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include <map>
#include <string>
#include <memory>

#include "FlowLaw.hh"
#include "pism/util/ConfigInterface.hh"

namespace pism {
namespace rheology {

#define ICE_ISOTHERMAL_GLEN  "isothermal_glen" /* Plain isothermal Glen */
#define ICE_PB      "pb"            /* Paterson-Budd (PatersonBudd) */
#define ICE_GPBLD   "gpbld"         /* Paterson-Budd-Lliboutry-Duval (GPBLD) */
#define ICE_HOOKE   "hooke"         /* Hooke (Hooke) */
#define ICE_ARR     "arr"           /* Temperature dependent Arrhenius (either warm or cold) */
#define ICE_GOLDSBY_KOHLSTEDT "gk"  /* Goldsby-Kohlstedt for SIA */
#define ICE_ARRWARM "arrwarm"       /* Temperature dependent Arrhenius (should be refactored into ICE_ARR) */

typedef FlowLaw*(*FlowLawCreator)(const std::string &,
                                  const Config &, EnthalpyConverter::Ptr);

class FlowLawFactory {
public:
  FlowLawFactory(const std::string &prefix,
                 Config::ConstPtr conf,
                 EnthalpyConverter::Ptr my_EC);
  ~FlowLawFactory();
  void set_default(const std::string &name);
  void add(const std::string &name, FlowLawCreator);
  void remove(const std::string &name);
  std::shared_ptr<FlowLaw> create();
private:
  std::string m_type_name, m_prefix;
  std::map<std::string, FlowLawCreator> m_flow_laws;
  const Config::ConstPtr m_config;
  EnthalpyConverter::Ptr m_EC;
};


} // end of namespace rheology
} // end of namespace pism

#endif  // __flowlaw_factory
