%{
#include "energy/EnthalpyModel.hh"
#include "energy/TemperatureModel.hh"
#include "regional/EnthalpyModel_Regional.hh"
#include "energy/enthSystem.hh"

#include "energy/tempSystem.hh"
#include "energy/BedrockColumn.hh"
#include "energy/utilities.hh"
%}

%rename(EnergyModelInputs) pism::energy::Inputs;

/* wrap the enthalpy solver to make testing easier */
%ignore pism::TridiagonalSystem::solve(unsigned int, double *);
%include "util/ColumnSystem.hh"

%rename(get_lambda) pism::energy::enthSystemCtx::lambda;
%include "energy/enthSystem.hh"

%rename(get_lambda) pism::energy::tempSystemCtx::lambda;
%include "energy/tempSystem.hh"

%shared_ptr(pism::energy::EnergyModel)
%include "energy/EnergyModel.hh"

%shared_ptr(pism::energy::EnthalpyModel)
%shared_ptr(pism::energy::DummyEnergyModel)
%feature("notabstract") pism::energy::EnthalpyModel;
%feature("notabstract") pism::energy::DummyEnergyModel;
%include "energy/EnthalpyModel.hh"

%shared_ptr(pism::energy::TemperatureModel)
%feature("notabstract") pism::energy::TemperatureModel;
%include "energy/TemperatureModel.hh"

%shared_ptr(pism::energy::EnthalpyModel_Regional)
%feature("notabstract") pism::energy::EnthalpyModel_Regional;
%include "regional/EnthalpyModel_Regional.hh"

%ignore pism::energy::BedrockColumn::solve(double, double, double, const double *, double *);
%include "energy/BedrockColumn.hh"

%include "energy/utilities.hh"
