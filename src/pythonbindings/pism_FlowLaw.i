%{
#include "rheology/FlowLaw.hh"
#include "rheology/GPBLD.hh"
#include "rheology/FlowLawFactory.hh"
%}

%shared_ptr(pism::rheology::FlowLaw)
%shared_ptr(pism::rheology::GPBLD)

%include "rheology/FlowLaw.hh"
%include "rheology/GPBLD.hh"

%include "rheology/FlowLawFactory.hh"
