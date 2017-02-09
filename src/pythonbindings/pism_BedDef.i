%{
#include "earth/PBLingleClark.hh"
#include "earth/greens.hh"
%}

%include "earth/PISMBedDef.hh"
%include "earth/PBLingleClark.hh"

%ignore pism::bed::ge_params;
%ignore pism::bed::ge_integrand;

%ignore pism::bed::vd_params;
%ignore pism::bed::vd_integrand;

%include "earth/greens.hh"
