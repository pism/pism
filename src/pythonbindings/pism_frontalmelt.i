%{
#include "coupler/frontalmelt/Constant.hh"
#include "coupler/frontalmelt/GivenClimate.hh"
%}

%shared_ptr(pism::frontalmelt::FrontalMeltModel)
%include "coupler/FrontalMeltModel.hh"

%shared_ptr(pism::frontalmelt::CompleteFrontalMeltModel)
%include "coupler/frontalmelt/CompleteFrontalMeltModel.hh"

%shared_ptr(pism::frontalmelt::Constant)
%rename(FrontalMeltConstant) pism::frontalmelt::Constant;
%include "coupler/frontalmelt/Constant.hh"

%shared_ptr(pism::frontalmelt::Given)
%rename(FrontalMeltGiven) pism::frontalmelt::Given;
%include "coupler/frontalmelt/GivenClimate.hh"
