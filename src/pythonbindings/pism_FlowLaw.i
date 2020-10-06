%{
#include "rheology/FlowLaw.hh"
#include "rheology/GPBLD.hh"
#include "rheology/FlowLawFactory.hh"
#include "rheology/PatersonBudd.hh"
#include "rheology/PatersonBuddCold.hh"
#include "rheology/PatersonBuddWarm.hh"
#include "rheology/grain_size_vostok.hh"
%}

%shared_ptr(pism::rheology::FlowLaw)
%shared_ptr(pism::rheology::GPBLD)
%shared_ptr(pism::rheology::PatersonBudd)
%shared_ptr(pism::rheology::PatersonBuddCold)
%shared_ptr(pism::rheology::PatersonBuddWarm)

%include "rheology/FlowLaw.hh"
%include "rheology/GPBLD.hh"
%include "rheology/PatersonBudd.hh"
%include "rheology/PatersonBuddCold.hh"
%include "rheology/PatersonBuddWarm.hh"

%include "rheology/FlowLawFactory.hh"

%include "rheology/grain_size_vostok.hh"
