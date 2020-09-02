%{
#include "earth/LingleClark.hh"
#include "earth/Given.hh"
#include "earth/greens.hh"
%}

%shared_ptr(pism::bed::BedDef)
%shared_ptr(pism::bed::Null)
%shared_ptr(pism::bed::PointwiseIsostasy)
%include "earth/BedDef.hh"
%shared_ptr(pism::bed::LingleClark)
%include "earth/LingleClark.hh"

%shared_ptr(pism::bed::Given)
%rename(GivenTopography) pism::bed::Given;
%include "earth/Given.hh"

%ignore pism::bed::ge_params;
%ignore pism::bed::ge_integrand;

%ignore pism::bed::vd_params;
%ignore pism::bed::vd_integrand;

%include "earth/greens.hh"
