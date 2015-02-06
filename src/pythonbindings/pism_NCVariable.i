/* NCVariable and NCSpatialVariable don't have default constructors. */
/* This should go before IceModelVec so that IceModelVec::metadata()
   is wrapped properly. */
%feature("valuewrapper") pism::NCVariable;
%feature("valuewrapper") pism::NCSpatialVariable;

%include "NCVariable.hh"
