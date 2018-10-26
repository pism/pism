%{
#include "coupler/frontal_melt/Constant.hh"
#include "coupler/frontal_melt/GivenClimate.hh"
%}

%shared_ptr(pism::frontalmelt::FrontalMeltModel)
%include "coupler/frontal_melt/Model.hh"

%shared_ptr(pism::frontalmelt::CompleteFrontalMeltModel)
%include "coupler/frontal_melt/CompleteFrontalMeltModel.hh"

%shared_ptr(pism::frontalmelt::Constant)
%rename(FrontalMeltConstant) pism::ocean::Constant;
%include "coupler/frontal_melt/Constant.hh"
