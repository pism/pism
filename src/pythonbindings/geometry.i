%{
#include "base/Geometry.hh"
#include "base/GeometryEvolution.hh"
%}


%shared_ptr(pism::Geometry)

// Treat data members of Geometry as read-only.
%feature("immutable", "1");
%include "base/Geometry.hh"
%feature("immutable", "0");

%shared_ptr(pism::GeometryEvolution)
%shared_ptr(pism::RegionalGeometryEvolution)
%include "base/GeometryEvolution.hh"
