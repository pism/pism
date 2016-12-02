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

%include "base/energy/EnergyModel.hh"
%include "base/energy/EnthalpyModel.hh"
%include "base/energy/TemperatureModel.hh"
%include "regional/EnthalpyModel_Regional.hh"
