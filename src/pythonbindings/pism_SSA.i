%{
#include "base/stressbalance/ssa/node_types.hh"
#include "base/stressbalance/ssa/SSAFEM.hh"
#include "base/stressbalance/ssa/SSAFD.hh"
#include "base/stressbalance/ssa/SSA_diagnostics.hh"
#include "base/stressbalance/ssa/SSAFD_diagnostics.hh"
%}

%include "base/stressbalance/ssa/node_types.hh"
%include "base/stressbalance/ShallowStressBalance.hh"
%include "base/stressbalance/ssa/SSA.hh"
%include "base/stressbalance/ssa/SSAFD.hh"
%include "base/stressbalance/ssa/SSAFEM.hh"

%template(SSADiag) pism::Diag<pism::stressbalance::SSA>;
%include "base/stressbalance/ssa/SSA_diagnostics.hh"
%template(SSAFDDiag) pism::Diag<pism::stressbalance::SSAFD>;
%include "base/stressbalance/ssa/SSAFD_diagnostics.hh"
