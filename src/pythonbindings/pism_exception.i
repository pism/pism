%header%{
#include "pism_python.hh"
#include "util/pism_utilities.hh"
%}

%exception {
  try {
    {
      pism::python::SigInstaller handler(SIGINT, pism::python::sigint_handler);
      $action
    }
    int sig = pism::python::check_signal();
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

/* FIXME: this completely overrides the %exception block above. I need
   to fix this. -- CK */
%include exception.i
%exception {
  try {
    $action
  } catch (pism::RuntimeError &e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  } catch (std::exception &e) {
    SWIG_exception(SWIG_UnknownError, e.what());
  } catch (...) {
    SWIG_exception(SWIG_UnknownError, "unknown C++ exception");
  }
}

