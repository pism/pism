// Copyright (C) 2009, 2010, 2011 Jed Brown, Ed Bueler and Constantine Khroulev
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

#include "flowlaw_factory.hh"
#include "pism_const.hh"

#undef ALEN
#define ALEN(a) (sizeof(a)/sizeof(a)[0])

IceFlowLawFactory::IceFlowLawFactory(MPI_Comm c,
                                     const char pre[], const NCConfigVariable &conf,
                                     EnthalpyConverter *my_EC)
  : comm(c), config(conf), EC(my_EC) {
  prefix[0] = 0;
  if (pre) {
    PetscStrncpy(prefix,pre,sizeof(prefix));
  }
  if (registerAll()) {
    PetscPrintf(comm,"IceFlowLawFactory::registerAll returned an error but we're in a constructor\n");
    PISMEnd();
  }
  if (setType(ICE_PB)) {       // Set's a default type
    PetscPrintf(comm,"IceFlowLawFactory::setType(\"%s\") returned an error, but we're in a constructor\n",ICE_PB);
    PISMEnd();
  }
}

IceFlowLawFactory::~IceFlowLawFactory()
{
  PetscFListDestroy(&type_list);
}

#undef __FUNCT__
#define __FUNCT__ "IceFlowLawFactory::registerType"
PetscErrorCode IceFlowLawFactory::registerType(const char tname[],
               PetscErrorCode(*icreate)(MPI_Comm,const char[],const NCConfigVariable &,
                                        EnthalpyConverter*, IceFlowLaw**))
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFListAdd(&type_list,tname,NULL,(void(*)(void))icreate);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


static PetscErrorCode create_isothermal_glen(MPI_Comm comm,const char pre[],
                                             const NCConfigVariable &config, EnthalpyConverter *EC, IceFlowLaw **i) {
  *i = new (IsothermalGlenIce)(comm, pre, config, EC);  return 0;
}
static PetscErrorCode create_pb(MPI_Comm comm,const char pre[],
                                const NCConfigVariable &config, EnthalpyConverter *EC, IceFlowLaw **i) {
  *i = new (ThermoGlenIce)(comm, pre, config, EC);  return 0;
}
static PetscErrorCode create_gpbld(MPI_Comm comm,const char pre[],
                                   const NCConfigVariable &config, EnthalpyConverter *EC, IceFlowLaw **i) {
  *i = new (GPBLDIce)(comm, pre, config, EC);  return 0;
}
static PetscErrorCode create_hooke(MPI_Comm comm,const char pre[],
                                   const NCConfigVariable &config, EnthalpyConverter *EC, IceFlowLaw **i) {
  *i = new (HookeIce)(comm, pre, config, EC);  return 0;
}
static PetscErrorCode create_arr(MPI_Comm comm,const char pre[],
                                 const NCConfigVariable &config, EnthalpyConverter *EC, IceFlowLaw **i) {
  *i = new (ThermoGlenArrIce)(comm, pre, config, EC);  return 0;
}
static PetscErrorCode create_arrwarm(MPI_Comm comm,const char pre[],
                                     const NCConfigVariable &config, EnthalpyConverter *EC, IceFlowLaw **i) {
  *i = new (ThermoGlenArrIceWarm)(comm, pre, config, EC);  return 0;
}
static PetscErrorCode create_hybrid(MPI_Comm comm,const char pre[],
                                    const NCConfigVariable &config, EnthalpyConverter *EC, IceFlowLaw **i) {
  *i = new (HybridIce)(comm, pre, config, EC);  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "IceFlowLawFactory::registerAll"
PetscErrorCode IceFlowLawFactory::registerAll()
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscMemzero(&type_list,sizeof(type_list));CHKERRQ(ierr);
  ierr = registerType(ICE_ISOTHERMAL_GLEN, &create_isothermal_glen); CHKERRQ(ierr);
  ierr = registerType(ICE_PB,     &create_pb);     CHKERRQ(ierr);
  ierr = registerType(ICE_GPBLD,  &create_gpbld);  CHKERRQ(ierr);
  ierr = registerType(ICE_HOOKE,  &create_hooke);  CHKERRQ(ierr);
  ierr = registerType(ICE_ARR,    &create_arr);    CHKERRQ(ierr);
  ierr = registerType(ICE_ARRWARM,&create_arrwarm);CHKERRQ(ierr);
  ierr = registerType(ICE_HYBRID, &create_hybrid); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "IceFlowLawFactory::setType"
PetscErrorCode IceFlowLawFactory::setType(const char type[])
{
  void (*r)(void);
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscFListFind(type_list, comm, type, PETSC_FALSE, (void(**)(void))&r);CHKERRQ(ierr);
  if (!r) {
    ierr = PetscPrintf(comm, "PISM ERROR: Selected ice type \"%s\" is not available.\n",type); CHKERRQ(ierr);
    PISMEnd();
  }
  ierr = PetscStrncpy(type_name,type,sizeof(type_name));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IceFlowLawFactory::setFromOptions"
PetscErrorCode IceFlowLawFactory::setFromOptions()
{
  PetscErrorCode ierr;
  PetscBool flg;
  char my_type_name[256];

  PetscFunctionBegin;
  // These options will choose Goldsby-Kohlstedt ice by default (see IceModel::setFromOptions()) but if a derived class
  // uses a different initialization procedure, we'll recognize them here as well.  A better long-term solution would be
  // to separate tracking of grain size from a particular flow law (since in principle they are unrelated) but since
  // HYBRID is the only one that currently uses grain size, this solution is acceptable.
  ierr = PetscOptionsHasName(prefix, "-gk_age", &flg); CHKERRQ(ierr);
  if (flg) {
    ierr = setType(ICE_HYBRID);CHKERRQ(ierr);
  }
  // -gk 0 does not make sense, so using PetscOptionsHasName is OK.
  ierr = PetscOptionsHasName(prefix, "-gk", &flg); CHKERRQ(ierr);
  if (flg) {
    ierr = setType(ICE_HYBRID);CHKERRQ(ierr);
  }
  ierr = PetscOptionsBegin(comm, prefix, "IceFlowLawFactory options", "IceFlowLaw");CHKERRQ(ierr);
  {
    ierr = PetscOptionsList("-flow_law","Flow law type","IceFlowLawFactory::setType",
                            type_list,type_name,my_type_name,sizeof(my_type_name),&flg);CHKERRQ(ierr);
    if (flg) {ierr = setType(my_type_name);CHKERRQ(ierr);}
  }
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "IceFlowLawFactory::create"
PetscErrorCode IceFlowLawFactory::create(IceFlowLaw **inice)
{
  PetscErrorCode ierr,(*r)(MPI_Comm,const char[],const NCConfigVariable &,EnthalpyConverter*,IceFlowLaw**);
  IceFlowLaw *ice;

  PetscFunctionBegin;
  PetscValidPointer(inice,3);
  *inice = 0;
  // find the function that can create selected ice type:
  ierr = PetscFListFind(type_list, comm, type_name, PETSC_FALSE, (void(**)(void))&r);CHKERRQ(ierr);
  if (!r) SETERRQ1(comm, 1,"Selected Ice type %s not available, but we shouldn't be able to get here anyway",type_name);
  // create an IceFlowLaw instance:
  ierr = (*r)(comm, prefix, config, EC, &ice);CHKERRQ(ierr);
  *inice = ice;
  PetscFunctionReturn(0);
}
