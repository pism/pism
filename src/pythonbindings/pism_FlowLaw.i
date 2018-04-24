%{
#include "rheology/FlowLaw.hh"
#include "rheology/GPBLD3.hh"
#include "rheology/GPBLD.hh"
#include "rheology/FlowLawFactory.hh"
%}

%shared_ptr(pism::rheology::FlowLaw)
%shared_ptr(pism::rheology::GPBLD)
%shared_ptr(pism::rheology::GPBLD3)

%include "rheology/FlowLaw.hh"
%include "rheology/GPBLD.hh"
%include "rheology/GPBLD3.hh"

%include "rheology/FlowLawFactory.hh"
