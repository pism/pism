%{
#include "util/Timeseries.hh"
using pism::Time;
using pism::TimeseriesMetadata;
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

%include "util/Timeseries.hh"
