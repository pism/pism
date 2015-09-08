%{
#include "coupler/ocean/POConstant.hh"
%}

%include "coupler/PISMOcean.hh"

%rename(OceanConstant) pism::ocean::Constant;
%include "coupler/ocean/POConstant.hh"

