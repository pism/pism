// Support for nc_types (e.g. NC_BYTE, etc).  In NetCDF3, an nc_type is an enum, and in
// NetCDF4 it is typedef'ed to be an int. The enums pose a small problem in C++ because
// you can't assign an arbitrary integer to an enum without a cast, and you can't assume
// even in C that an enum is an int, so you have to be careful about pointers to enums
// versus pointers to ints.  Moreover, I don't know how to grab the definitions from
// netcdf.h here without wrapping everything in the file.
//
// So: we assume that nc_type is an enum.  On input, we force the python input to be an int,
// use pointers to the int variable where needed, and then do a static cast to shove the int
// into an nc_type.  This procedure works correctly if nc_type is an int instead of an enum.
// As for the allowed values, we copy the defines from (NetCDF4) netcdf.h.  No typechecking
// is being done to ensure that a python int on input is a valid nc_type, which isn't good.
// In particular, the allowed values are different in NetCDF4 vs. NetCDF3 (there are more of them.)
// A constraint check to the minimal set of NetCDF3 types would be the right thing to do. (FIXME)
%typemap(in) IO_Type (int tmp){
    SWIG_AsVal(int)($input,&tmp);
    $1 = static_cast<pism::IO_Type>(tmp);
}
%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) IO_Type {
    $1 = PyInt_Check($input);
}

// same for IO_Mode
%typemap(in) IO_Mode (int tmp){
    SWIG_AsVal(int)($input,&tmp);
    $1 = static_cast<pism::IO_Mode>(tmp);
}
%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) IO_Mode {
    $1 = PyInt_Check($input);
}

%include "IO_Flags.hh"
%include "PIO.hh"               // include before NCVariable

%{
#include "PIO.hh"
%}
