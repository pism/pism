%{
#include "base/hydrology/PISMHydrology.hh"
%}

%rename(DistributedHydrology) pism::hydrology::Distributed;
%rename(RoutingHydrology) pism::hydrology::Routing;
%rename(NullTransportHydrology) pism::hydrology::NullTransport;

%include "base/hydrology/PISMHydrology.hh"
