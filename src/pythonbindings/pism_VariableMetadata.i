/* VariableMetadata and SpatialVariableMetadata don't have default constructors. */
/* This should go before IceModelVec so that IceModelVec::metadata()
   is wrapped properly. */
%feature("valuewrapper") pism::VariableMetadata;
%feature("valuewrapper") pism::SpatialVariableMetadata;
%feature("valuewrapper") pism::TimeseriesMetadata;

%include "util/VariableMetadata.hh"
