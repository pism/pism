%{
#include "base/stressbalance/ssa/SSAFEM.hh"
#include "base/stressbalance/ssa/SSAFD.hh"
#include "verif/tests/exactTestsIJ.h"
%}

%include "base/stressbalance/ShallowStressBalance.hh"
%include "base/stressbalance/ssa/SSA.hh"
%include "base/stressbalance/ssa/SSAFD.hh"
%include "base/stressbalance/ssa/SSAFEM.hh"

/* Wrap C code implementing exact solutions for SSA verification
 * tests.
 */
%include "verif/tests/exactTestsIJ.h"
