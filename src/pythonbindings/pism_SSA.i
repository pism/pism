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

%shared_ptr(pism::Diag<pism::stressbalance::SSA>)
%template(SSADiag) pism::Diag<pism::stressbalance::SSA>;
%shared_ptr(pism::stressbalance::SSA_taud_mag)
%shared_ptr(pism::stressbalance::SSA_taud)
%shared_ptr(pism::stressbalance::SSA_calving_front_pressure_difference)
%include "base/stressbalance/ssa/SSA_diagnostics.hh"
%shared_ptr(pism::Diag<pism::stressbalance::SSAFD>)
%template(SSAFDDiag) pism::Diag<pism::stressbalance::SSAFD>;
%shared_ptr(pism::stressbalance::SSAFD_nuH)
%include "base/stressbalance/ssa/SSAFD_diagnostics.hh"
