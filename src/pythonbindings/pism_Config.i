%{
#include "util/Config.hh"
%}

%extend pism::Config
{
  %pythoncode "Config.py"
}
