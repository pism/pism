%{
#include "coupler/DebrisModel.hh"
#include "coupler/debris/GivenDebris.hh"
#include "coupler/debris/IceMeltEnhancement.hh"
#include "coupler/debris/DebrisInput.hh"
#include "coupler/debris/EnglacialTransport.hh"
#include "coupler/debris/SupraglacialTransport.hh"
#include "coupler/debris/GravitationalTransport.hh"
#include "coupler/debris/TerminusRemoval.hh"
#include "coupler/debris/DebrisTransport.hh"
#include "coupler/debris/Factory.hh"
%}

%rename(DebrisInputs) pism::debris::Inputs;
%shared_ptr(pism::debris::DebrisModel)
%include "coupler/DebrisModel.hh"

%shared_ptr(pism::debris::Given)
%rename(DebrisGiven) pism::debris::Given;
%include "coupler/debris/GivenDebris.hh"

%shared_ptr(pism::debris::IceMeltEnhancement)
%include "coupler/debris/IceMeltEnhancement.hh"

%shared_ptr(pism::debris::DebrisInput)
%include "coupler/debris/DebrisInput.hh"

%shared_ptr(pism::debris::EnglacialTransport)
%include "coupler/debris/EnglacialTransport.hh"

%shared_ptr(pism::debris::SupraglacialTransport)
%include "coupler/debris/SupraglacialTransport.hh"

%shared_ptr(pism::debris::GravitationalTransport)
%include "coupler/debris/GravitationalTransport.hh"

%shared_ptr(pism::debris::TerminusRemoval)
%include "coupler/debris/TerminusRemoval.hh"

%shared_ptr(pism::debris::DebrisTransport)
%include "coupler/debris/DebrisTransport.hh"

%shared_ptr(pism::debris::Factory)
%rename(DebrisFactory) pism::debris::Factory;
%include "coupler/debris/Factory.hh"
