%{
#include "stressbalance/ssa/SSAFEM.hh"
#include "stressbalance/ssa/SSAFD.hh"
%}

%include "stressbalance/ShallowStressBalance.hh"
// The template used in SSA.hh needs to be instantiated in SWIG before
// it is used.
%template(Diag_SSA) pism::Diag<pism::SSA>;
%include "stressbalance/ssa/SSA.hh"
%ignore pism::SSAFEFunction;
%ignore pism::SSAFEJacobian;
%include "stressbalance/ssa/SSAFEM.hh"
%template(Diag_SSAFD) pism::Diag<pism::SSAFD>;
%include "stressbalance/ssa/SSAFD.hh"

/* Wrap C code implementing exact solutions for SSA verification
 * tests.
 */

// Tell SWIG that input arguments of type double * are to be treated as return values,
// and that int return values are to be error checked as per a PetscErrorCode.
%apply double *OUTPUT  {double *};
%typemap(out,noblock=1) int {
PyPetsc_ChkErrQ($1); %set_output(VOID_Object);
}
%include "exactTestsIJ.h"
// FIXME! I don't know how to undo the output typemap.
// %typemap(out,noblock=1) int = PREVIOUS;
%clear double *;
