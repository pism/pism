%{
#include "util/petscwrappers/DM.hh"
%}

/* This is needed to wrap IceGrid::get_dm() */
%shared_ptr(pism::petsc::DM)
%shared_ptr(pism::petsc::Wrapper< ::DM >)

%include "util/petscwrappers/Wrapper.hh"
%template(_DMWrapper) pism::petsc::Wrapper< ::DM >;

%include "util/petscwrappers/DM.hh"
