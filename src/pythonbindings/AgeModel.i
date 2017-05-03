%{
#include "base/age/AgeModel.hh"
#include "base/age/AgeColumnSystem.hh"
%}

%shared_ptr(pism::AgeModel)
%include "base/age/AgeModel.hh"
%include "base/age/AgeColumnSystem.hh"
