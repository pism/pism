%{
#include "hydrology/Hydrology.hh"
#include "hydrology/NullTransport.hh"
#include "hydrology/Routing.hh"
#include "hydrology/SteadyState.hh"
#include "hydrology/Distributed.hh"
#include "hydrology/EmptyingProblem.hh"
%}

%rename(DistributedHydrology) pism::hydrology::Distributed;
%rename(RoutingHydrology) pism::hydrology::Routing;
%rename(RoutingSteadyHydrology) pism::hydrology::RoutingSteady;
%rename(NullTransportHydrology) pism::hydrology::NullTransport;

%rename(HydrologyInputs) pism::hydrology::Inputs;

%shared_ptr(pism::hydrology::Hydrology)
%shared_ptr(pism::hydrology::NullTransport)
%shared_ptr(pism::hydrology::Routing)
%shared_ptr(pism::hydrology::SteadyState)
%shared_ptr(pism::hydrology::Distributed)
%include "hydrology/Hydrology.hh"
%include "hydrology/NullTransport.hh"
%include "hydrology/Routing.hh"
%include "hydrology/SteadyState.hh"
%include "hydrology/Distributed.hh"

%shared_ptr(pism::hydrology::EmptyingProblem)
%include "hydrology/EmptyingProblem.hh"
