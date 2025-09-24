/* VariableMetadata does not have the default constructor. */
/* This should go before array::Array so that array::Array::metadata()
   is wrapped properly. */
%feature("valuewrapper") pism::VariableMetadata;

%ignore pism::Attribute::operator=;
%ignore pism::VariableMetadata::operator[];

%include "util/VariableMetadata.hh"
