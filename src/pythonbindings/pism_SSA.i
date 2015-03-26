%{
#include "stressbalance/ssa/SSAFEM.hh"
#include "stressbalance/ssa/SSAFD.hh"
#include "exactTestsIJ.h"
%}

%include "stressbalance/ShallowStressBalance.hh"
%include "stressbalance/ssa/SSA.hh"
%include "stressbalance/ssa/SSAFD.hh"
%include "stressbalance/ssa/SSAFEM.hh"

/* Wrap C code implementing exact solutions for SSA verification
 * tests.
 */
%include "exactTestsIJ.h"
