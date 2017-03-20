%{
#include "base/energy/EnthalpyModel.hh"
#include "base/energy/TemperatureModel.hh"
#include "regional/EnthalpyModel_Regional.hh"
#include "base/energy/enthSystem.hh"
using pism::MaskValue;
#include "base/energy/tempSystem.hh"
%}

/* wrap the enthalpy solver to make testing easier */
%include "base/columnSystem.hh"
%rename(get_lambda) pism::energy::enthSystemCtx::lambda;
%include "base/energy/enthSystem.hh"
%rename(get_lambda) pism::energy::tempSystemCtx::lambda;
%include "base/energy/tempSystem.hh"

%shared_ptr(pism::energy::EnergyModel)
%include "base/energy/EnergyModel.hh"
%shared_ptr(pism::energy::EnthalpyModel)
%shared_ptr(pism::energy::DummyEnergyModel)
%include "base/energy/EnthalpyModel.hh"
%shared_ptr(pism::energy::TemperatureModel)
%include "base/energy/TemperatureModel.hh"
%shared_ptr(pism::energy::EnthalpyModel_Regional)
%include "regional/EnthalpyModel_Regional.hh"
