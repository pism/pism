%header%{
#include "pism_python_signal.hh"
#include "pism_const.hh"
%}

%exception {
  try {
    {
      SigInstaller handler(SIGINT, pism_sigint_handler); 
      $action
        }
    int sig = pism_check_signal();
    if (sig == SIGINT) {
      PyErr_SetString(PyExc_KeyboardInterrupt, "");
      return NULL;
    } else if (sig) {
      SWIG_exception(SWIG_RuntimeError, "Caught an unknown signal.");
      return NULL;      
    } 
  }
  catch(Swig::DirectorException &e) {
    SWIG_exception(SWIG_RuntimeError, e.getMessage());
  }
  catch (...) {
    SWIG_exception(SWIG_RuntimeError, "Caught an unexpected C++ exception");
  }
 }
