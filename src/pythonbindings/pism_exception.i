%header%{
#include "pism_python_signal.hh"
#include "pism_const.hh"
%}

%exception {
   try
   {
      {
        SigInstaller h(SIGINT,pism_sigint_handler); 
        $action
      }
      int sig = pism_check_signal();
      if( sig == SIGINT )
      {
        PyErr_SetString(PyExc_KeyboardInterrupt, "");
        return NULL;
      }
      else if(sig)
      {
        PyErr_SetString(PyExc_Exception, "Caught an unknown signal.");
        return NULL;      
      } 
   }
   catch(Swig::DirectorException &e) {
     PyErr_SetString(PyExc_RuntimeError,e.getMessage());
   }
   catch (...) {
     pism::verbPrintf(1,PETSC_COMM_WORLD,"Caught a C++ exception; PISM memory may be corrupt. Aborting.");
     pism::PISMEnd();
   }
}
