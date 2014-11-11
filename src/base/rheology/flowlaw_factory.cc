// Copyright (C) 2009, 2010, 2011, 2012, 2013, 2014 Jed Brown, Ed Bueler and Constantine Khroulev
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

  if (registerAll()) {
    throw std::runtime_error("IceFlowLawFactory constructor: failed to register flow law types");
  }

  if (setType(ICE_PB)) {        // Set's a default type
    throw std::runtime_error("IceFlowLawFactory constructor: failed to set the default flow law type");
  }
}

IceFlowLawFactory::~IceFlowLawFactory()
{
}

PetscErrorCode IceFlowLawFactory::registerType(const std::string &name, IceFlowLawCreator icreate)
{
  flow_laws[name] = icreate;
  return 0;
}

PetscErrorCode IceFlowLawFactory::removeType(const std::string &name) {
  flow_laws.erase(name);
  return 0;
}


static PetscErrorCode create_isothermal_glen(MPI_Comm com,const std::string &pre,
                                             const Config &config, EnthalpyConverter *EC, IceFlowLaw **i) {
  *i = new (IsothermalGlenIce)(com, pre, config, EC);
  return 0;
}

static PetscErrorCode create_pb(MPI_Comm com,const std::string &pre,
                                const Config &config, EnthalpyConverter *EC, IceFlowLaw **i) {
  *i = new (ThermoGlenIce)(com, pre, config, EC);
  return 0;
}

static PetscErrorCode create_gpbld(MPI_Comm com,const std::string &pre,
                                   const Config &config, EnthalpyConverter *EC, IceFlowLaw **i) {
  *i = new (GPBLDIce)(com, pre, config, EC);
  return 0;
}

static PetscErrorCode create_hooke(MPI_Comm com,const std::string &pre,
                                   const Config &config, EnthalpyConverter *EC, IceFlowLaw **i) {
  *i = new (HookeIce)(com, pre, config, EC);
  return 0;
}

static PetscErrorCode create_arr(MPI_Comm com,const std::string &pre,
                                 const Config &config, EnthalpyConverter *EC, IceFlowLaw **i) {
  *i = new (ThermoGlenArrIce)(com, pre, config, EC);
  return 0;
}

static PetscErrorCode create_arrwarm(MPI_Comm com,const std::string &pre,
                                     const Config &config, EnthalpyConverter *EC, IceFlowLaw **i) {
  *i = new (ThermoGlenArrIceWarm)(com, pre, config, EC);
  return 0;
}

static PetscErrorCode create_goldsby_kohlstedt(MPI_Comm com,const std::string &pre,
                                               const Config &config, EnthalpyConverter *EC, IceFlowLaw **i) {
  *i = new (GoldsbyKohlstedtIce)(com, pre, config, EC);
  return 0;
}

PetscErrorCode IceFlowLawFactory::registerAll()
{
  PetscErrorCode ierr;

  flow_laws.clear();
  ierr = registerType(ICE_ISOTHERMAL_GLEN, &create_isothermal_glen); CHKERRQ(ierr);
  ierr = registerType(ICE_PB, &create_pb);     CHKERRQ(ierr);
  ierr = registerType(ICE_GPBLD, &create_gpbld);  CHKERRQ(ierr);
  ierr = registerType(ICE_HOOKE, &create_hooke);  CHKERRQ(ierr);
  ierr = registerType(ICE_ARR, &create_arr);    CHKERRQ(ierr);
  ierr = registerType(ICE_ARRWARM, &create_arrwarm); CHKERRQ(ierr);
  ierr = registerType(ICE_GOLDSBY_KOHLSTEDT, &create_goldsby_kohlstedt); CHKERRQ(ierr);

  return 0;
}

PetscErrorCode IceFlowLawFactory::setType(const std::string &type)
{
  IceFlowLawCreator r = flow_laws[type];
  if (not r) {
    throw RuntimeError::formatted("Selected ice type \"%s\" is not available.\n",
                                  type.c_str());
  }

  type_name = type;

  return 0;
}

PetscErrorCode IceFlowLawFactory::setFromOptions()
{
  PetscErrorCode ierr;
  bool flag;
  std::string my_type_name;

  ierr = PetscOptionsBegin(com, prefix.c_str(), "IceFlowLawFactory options", "IceFlowLaw");
  PISM_PETSC_CHK(ierr, "PetscOptionsBegin");
  {

    // build the list of choices
    std::map<std::string,IceFlowLawCreator>::iterator j = flow_laws.begin();
    std::set<std::string> choices;
    while (j != flow_laws.end()) {
      choices.insert(j->first);
      ++j;
    }

    ierr = OptionsList("-flow_law", "flow law type", choices,
                       type_name, my_type_name, flag); CHKERRQ(ierr);

    if (flag) {
      ierr = setType(my_type_name); CHKERRQ(ierr);
    }

  }
  ierr = PetscOptionsEnd();
  PISM_PETSC_CHK(ierr, "PetscOptionsEnd");

  return 0;
}

PetscErrorCode IceFlowLawFactory::create(IceFlowLaw **inice)
{
  PetscErrorCode ierr;
  IceFlowLawCreator r;
  IceFlowLaw *ice;

  PetscFunctionBegin;
#if PETSC_VERSION_LT(3,4,0)
  PetscValidPointer(inice,3);
#endif
  *inice = 0;

  // find the function that can create selected ice type:
  r = flow_laws[type_name];
  if (r == NULL) {
    throw RuntimeError::formatted("Selected ice type %s is not available,\n"
                                  "but we shouldn't be able to get here anyway",
                                  type_name.c_str());
  }

  // create an IceFlowLaw instance:
  ierr = (*r)(com, prefix, config, EC, &ice); CHKERRQ(ierr);
  *inice = ice;

  PetscFunctionReturn(0);
}

} // end of namespace pism
