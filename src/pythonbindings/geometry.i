%{
#include "base/Geometry.hh"
#include "base/GeometryEvolution.hh"
%}


%shared_ptr(pism::Geometry)
%include "base/Geometry.hh"

%shared_ptr(pism::GeometryEvolution)
%shared_ptr(pism::RegionalGeometryEvolution)
%include "base/GeometryEvolution.hh"
