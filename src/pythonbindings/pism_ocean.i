%{
#include "coupler/ocean/POConstant.hh"
%}

%template(OceanDiagnostic) pism::Diag<pism::ocean::OceanModel>;
%include "coupler/PISMOcean.hh"

%rename(OceanConstant) pism::ocean::Constant;
%include "coupler/ocean/POConstant.hh"

