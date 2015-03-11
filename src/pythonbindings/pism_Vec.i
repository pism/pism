%{
#include "Vec.hh"
%}

/* This may be needed to use PISM's scatters to/from processor 0. */
%shared_ptr(pism::petsc::Vec)
%shared_ptr(pism::petsc::TemporaryGlobalVec)
%shared_ptr(pism::petsc::Wrapper< ::Vec >)

%include "Wrapper.hh"
%template(VecWrapper) pism::petsc::Wrapper< ::Vec >;

%include "Vec.hh"
