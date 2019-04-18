%{
#include "coupler/frontalmelt/Constant.hh"
#include "coupler/frontalmelt/DischargeRouting.hh"
#include "coupler/frontalmelt/DischargeGiven.hh"
#include "coupler/frontalmelt/Given.hh"

#include "coupler/frontalmelt/FrontalMeltPhysics.hh"
%}

%shared_ptr(pism::frontalmelt::FrontalMelt)
%include "coupler/FrontalMelt.hh"

%shared_ptr(pism::frontalmelt::Constant)
%rename(FrontalMeltConstant) pism::frontalmelt::Constant;
%include "coupler/frontalmelt/Constant.hh"

%shared_ptr(pism::frontalmelt::DischargeRouting)
%rename(FrontalMeltDischargeRouting) pism::frontalmelt::DischargeRouting;
%include "coupler/frontalmelt/DischargeRouting.hh"

%shared_ptr(pism::frontalmelt::Given)
%rename(FrontalMeltGiven) pism::frontalmelt::Given;
%include "coupler/frontalmelt/Given.hh"

%shared_ptr(pism::frontalmelt::DischargeGiven)
%rename(FrontalMeltDischargeGiven) pism::frontalmelt::DischargeGiven;
%include "coupler/frontalmelt/DischargeGiven.hh"

%shared_ptr(pism::frontalmelt::FrontalMeltPhysics)
%include "coupler/frontalmelt/FrontalMeltPhysics.hh"
