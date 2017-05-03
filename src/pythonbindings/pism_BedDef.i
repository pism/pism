%{
#include "earth/PBLingleClark.hh"
#include "earth/greens.hh"
%}

%shared_ptr(pism::bed::BedDef)
%shared_ptr(pism::bed::PBNull)
%shared_ptr(pism::bed::PBPointwiseIsostasy)
%include "earth/PISMBedDef.hh"
%shared_ptr(pism::bed::PBLingleClark)
%include "earth/PBLingleClark.hh"

%ignore pism::bed::ge_params;
%ignore pism::bed::ge_integrand;

%ignore pism::bed::vd_params;
%ignore pism::bed::vd_integrand;

%include "earth/greens.hh"
