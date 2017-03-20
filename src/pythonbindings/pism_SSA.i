%{
#include "base/stressbalance/ssa/SSAFEM.hh"
#include "base/stressbalance/ssa/SSAFD.hh"
#include "base/stressbalance/ssa/SSA_diagnostics.hh"
#include "base/stressbalance/ssa/SSAFD_diagnostics.hh"
%}

%shared_ptr(pism::stressbalance::ShallowStressBalance)
%shared_ptr(pism::stressbalance::ZeroSliding)
%shared_ptr(pism::stressbalance::PrescribedSliding)
%include "base/stressbalance/ShallowStressBalance.hh"

%shared_ptr(pism::stressbalance::SSA)
%include "base/stressbalance/ssa/SSA.hh"
%shared_ptr(pism::stressbalance::SSAFD)
%include "base/stressbalance/ssa/SSAFD.hh"
%shared_ptr(pism::stressbalance::SSAFEM)
%include "base/stressbalance/ssa/SSAFEM.hh"

%template(SSADiag) pism::Diag<pism::stressbalance::SSA>;
%include "base/stressbalance/ssa/SSA_diagnostics.hh"
%template(SSAFDDiag) pism::Diag<pism::stressbalance::SSAFD>;
%include "base/stressbalance/ssa/SSAFD_diagnostics.hh"
