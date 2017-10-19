%{
#include "age/AgeModel.hh"
#include "age/AgeColumnSystem.hh"
%}

%shared_ptr(pism::AgeModel)
%include "age/AgeModel.hh"
%include "age/AgeColumnSystem.hh"
