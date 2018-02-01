%{
#include "util/pism_options.hh"
%}

%include "util/options.hh"
// instantiate templates used by option processing classes below
%template(_OptionStdString) pism::options::Option<std::string>;
%template(_OptionStdVectorStdString) pism::options::Option<std::vector<std::string> >;
%template(_OptionStdSetStdString) pism::options::Option<std::set<std::string> >;
%template(_OptionInt) pism::options::Option<int>;
%template(_OptionStdVectorInt) pism::options::Option<std::vector<int> >;
%template(_OptionDouble) pism::options::Option<double>;
%template(_OptionStdVectorDouble) pism::options::Option<std::vector<double> >;

// rename classes in pism::options (SWIG flattens namespaces)
%rename(OptionString) pism::options::String;
%rename(OptionStringList) pism::options::StringList;
%rename(OptionStringSet) pism::options::StringSet;
%rename(OptionKeyword) pism::options::Keyword;
%rename(OptionInteger) pism::options::Integer;
%rename(OptionIntegerList) pism::options::IntegerList;
%rename(OptionReal) pism::options::Real;
%rename(OptionRealList) pism::options::RealList;
%rename(OptionBool) pism::options::Bool;

%ignore pism::options::StringList::operator[];
%ignore pism::options::IntegerList::operator[];
%ignore pism::options::RealList::operator[];

%include "util/pism_options.hh"
