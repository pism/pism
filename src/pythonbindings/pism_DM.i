%{
#include "DM.hh"
%}

/* This is needed to wrap IceGrid::get_dm() */
%shared_ptr(pism::PISMDM)
%shared_ptr(pism::petsc::Wrapper<DM>)

%include "Wrapper.hh"
%template(DMWrapper) pism::petsc::Wrapper<DM>;

%include "DM.hh"
