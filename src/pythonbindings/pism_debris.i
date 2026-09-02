%{
#include "coupler/DebrisModel.hh"
#include "coupler/debris/GivenDebris.hh"
#include "coupler/debris/IceMeltEnhancement.hh"
%}

%shared_ptr(pism::debris::DebrisModel)
%include "coupler/DebrisModel.hh"

%shared_ptr(pism::debris::Given)
%rename(DebrisGiven) pism::debris::Given;
%include "coupler/debris/GivenDebris.hh"

%shared_ptr(pism::debris::IceMeltEnhancement)
%include "coupler/debris/IceMeltEnhancement.hh"
