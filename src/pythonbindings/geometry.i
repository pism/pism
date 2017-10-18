%{
#include "base/geometry/Geometry.hh"
#include "base/geometry/GeometryEvolution.hh"
%}


%shared_ptr(pism::Geometry)

// Treat data members of Geometry as read-only.
%feature("immutable", "1");
%include "base/geometry/Geometry.hh"
%feature("immutable", "0");

%shared_ptr(pism::GeometryEvolution)
%shared_ptr(pism::RegionalGeometryEvolution)
%include "base/geometry/GeometryEvolution.hh"
