%{
#include "base/Geometry.hh"
#include "base/GeometryEvolution.hh"
%}


%shared_ptr(pism::Geometry)
%include "base/Geometry.hh"

%shared_ptr(pism::GeometryEvolution)
%include "base/GeometryEvolution.hh"
