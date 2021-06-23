%{
#include "util/petscwrappers/Vec.hh"
%}

/* This may be needed to use PISM's scatters to/from processor 0. */
%shared_ptr(pism::petsc::Vec)
%shared_ptr(pism::petsc::TemporaryGlobalVec)
%shared_ptr(pism::Wrapper< ::Vec >)

%include "util/Wrapper.hh"
%template(_VecWrapper) pism::Wrapper< ::Vec >;

%include "util/petscwrappers/Vec.hh"
