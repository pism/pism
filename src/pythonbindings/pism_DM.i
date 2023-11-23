%{
#include "util/petscwrappers/DM.hh"
%}

/* This is needed to wrap IceGrid::get_dm() */
%shared_ptr(pism::petsc::DM)
%shared_ptr(pism::Wrapper< ::DM >)

%include "util/Wrapper.hh"
%template(_DMWrapper) pism::Wrapper< ::DM >;

%include "util/petscwrappers/DM.hh"
