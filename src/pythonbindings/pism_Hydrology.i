%{
#include "hydrology/Hydrology.hh"
%}

%rename(DistributedHydrology) pism::hydrology::Distributed;
%rename(RoutingHydrology) pism::hydrology::Routing;
%rename(NullTransportHydrology) pism::hydrology::NullTransport;

%shared_ptr(pism::hydrology::Hydrology)
%shared_ptr(pism::hydrology::NullTransport)
%shared_ptr(pism::hydrology::Routing)
%shared_ptr(pism::hydrology::Distributed)
%include "hydrology/Hydrology.hh"
