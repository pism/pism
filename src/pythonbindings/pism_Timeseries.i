%{
#include "Timeseries.hh"
using pism::Time;
using pism::NCTimeseries;
%}

%extend pism::Timeseries
{
    %ignore operator[];
    double getitem(unsigned int i)
    {
        return (*$self)[i];
    }
    
    %pythoncode {
    def __getitem__(self,*args):
        return self.getitem(*args)
    }
};

%include "Timeseries.hh"
