%{
#include "energy/EnthalpyModel.hh"
#include "energy/TemperatureModel.hh"
#include "regional/EnthalpyModel_Regional.hh"
#include "energy/enthSystem.hh"
using pism::MaskValue;
#include "energy/tempSystem.hh"
%}

%rename(EnergyModelInputs) pism::energy::Inputs;

/* wrap the enthalpy solver to make testing easier */
%include "util/ColumnSystem.hh"

%rename(get_lambda) pism::energy::enthSystemCtx::lambda;
%include "energy/enthSystem.hh"

%rename(get_lambda) pism::energy::tempSystemCtx::lambda;
%include "energy/tempSystem.hh"

%shared_ptr(pism::energy::EnergyModel)
%include "energy/EnergyModel.hh"

%shared_ptr(pism::energy::EnthalpyModel)
%shared_ptr(pism::energy::DummyEnergyModel)
%include "energy/EnthalpyModel.hh"

%shared_ptr(pism::energy::TemperatureModel)
%include "energy/TemperatureModel.hh"

%shared_ptr(pism::energy::EnthalpyModel_Regional)
%include "regional/EnthalpyModel_Regional.hh"
