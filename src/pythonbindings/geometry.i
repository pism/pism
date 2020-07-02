%{
#include "geometry/Geometry.hh"
#include "geometry/GeometryEvolution.hh"
%}


%shared_ptr(pism::Geometry)

// Treat data members of Geometry as read-only.
%feature("immutable", "1");
%include "geometry/Geometry.hh"
%feature("immutable", "0");

%shared_ptr(pism::GeometryEvolution)
%shared_ptr(pism::RegionalGeometryEvolution)
%include "geometry/GeometryEvolution.hh"

pism_class(pism::MPDATA2, "pism/geometry/MPDATA2.hh")
