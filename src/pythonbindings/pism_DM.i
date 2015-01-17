%{
#include "DM.hh"
%}

/* This is needed to wrap IceGrid::get_dm() */
%shared_ptr(pism::PISMDM)
%ignore pism::PISMDM::operator DM;

%include "DM.hh"
