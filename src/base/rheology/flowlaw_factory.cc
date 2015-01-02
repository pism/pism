// Copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "flowlaw_factory.hh"
#include "pism_const.hh"
#include "pism_options.hh"
#include "PISMUnits.hh"
#include <cassert>
#include <stdexcept>

#include "error_handling.hh"

namespace pism {

IceFlowLawFactory::IceFlowLawFactory(MPI_Comm c,
                                     const std::string &pre,
                                     const Config &conf,
                                     EnthalpyConverter *my_EC)
  : com(c), config(conf), EC(my_EC) {

  prefix = pre;

  assert(prefix.empty() == false);

  registerAll();

  setType(ICE_PB);    
}

IceFlowLawFactory::~IceFlowLawFactory()
{
}

void IceFlowLawFactory::registerType(const std::string &name, IceFlowLawCreator icreate)
{
  flow_laws[name] = icreate;
}

void IceFlowLawFactory::removeType(const std::string &name) {
  flow_laws.erase(name);
}


IceFlowLaw* create_isothermal_glen(MPI_Comm com,const std::string &pre,
                                             const Config &config, EnthalpyConverter *EC) {
  return new (IsothermalGlenIce)(com, pre, config, EC);
}

IceFlowLaw* create_pb(MPI_Comm com,const std::string &pre,
                                const Config &config, EnthalpyConverter *EC) {
  return new (ThermoGlenIce)(com, pre, config, EC);
}

IceFlowLaw* create_gpbld(MPI_Comm com,const std::string &pre,
                                   const Config &config, EnthalpyConverter *EC) {
  return new (GPBLDIce)(com, pre, config, EC);
}

IceFlowLaw* create_hooke(MPI_Comm com,const std::string &pre,
                                   const Config &config, EnthalpyConverter *EC) {
  return new (HookeIce)(com, pre, config, EC);
}

IceFlowLaw* create_arr(MPI_Comm com,const std::string &pre,
                                 const Config &config, EnthalpyConverter *EC) {
  return new (ThermoGlenArrIce)(com, pre, config, EC);
}

IceFlowLaw* create_arrwarm(MPI_Comm com,const std::string &pre,
                                     const Config &config, EnthalpyConverter *EC) {
  return new (ThermoGlenArrIceWarm)(com, pre, config, EC);
}

IceFlowLaw* create_goldsby_kohlstedt(MPI_Comm com,const std::string &pre,
                                               const Config &config, EnthalpyConverter *EC) {
  return new (GoldsbyKohlstedtIce)(com, pre, config, EC);
}

void IceFlowLawFactory::registerAll()
{
  flow_laws.clear();
  registerType(ICE_ISOTHERMAL_GLEN, &create_isothermal_glen);
  registerType(ICE_PB, &create_pb);
  registerType(ICE_GPBLD, &create_gpbld);
  registerType(ICE_HOOKE, &create_hooke);
  registerType(ICE_ARR, &create_arr);
  registerType(ICE_ARRWARM, &create_arrwarm);
  registerType(ICE_GOLDSBY_KOHLSTEDT, &create_goldsby_kohlstedt);

}

void IceFlowLawFactory::setType(const std::string &type)
{
  IceFlowLawCreator r = flow_laws[type];
  if (not r) {
    throw RuntimeError::formatted("Selected ice type \"%s\" is not available.\n",
                                  type.c_str());
  }

  type_name = type;

}

void IceFlowLawFactory::setFromOptions()
{
  {
    // build the list of choices
    std::map<std::string,IceFlowLawCreator>::iterator j = flow_laws.begin();
    std::vector<std::string> choices;
    while (j != flow_laws.end()) {
      choices.push_back(j->first);
      ++j;
    }

    options::Keyword type("-" + prefix + "flow_law", "flow law type",
                          join(choices, ","), type_name);

    if (type.is_set()) {
      setType(type);
    }
  }
}

IceFlowLaw* IceFlowLawFactory::create()
{
  // find the function that can create selected ice type:
  IceFlowLawCreator r = flow_laws[type_name];
  if (r == NULL) {
    throw RuntimeError::formatted("Selected ice type %s is not available,\n"
                                  "but we shouldn't be able to get here anyway",
                                  type_name.c_str());
  }

  // create an IceFlowLaw instance:
  return (*r)(com, prefix, config, EC);
}

} // end of namespace pism
