%{
#include "coupler/ocean/POConstant.hh"
%}

%template(OceanDiagnostic) pism::Diag<pism::ocean::OceanModel>;
%shared_ptr(pism::ocean::OceanModel)
%include "coupler/PISMOcean.hh"

%shared_ptr(pism::ocean::Constant)
%rename(OceanConstant) pism::ocean::Constant;
%include "coupler/ocean/POConstant.hh"
