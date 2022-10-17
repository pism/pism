%{
#include "stressbalance/ssa/SSAFEM.hh"
#include "stressbalance/ssa/SSAFD.hh"
#include "stressbalance/ssa/SSA_diagnostics.hh"
#include "stressbalance/ssa/SSAFD_diagnostics.hh"
#include "stressbalance/StressBalance.hh"
%}

%shared_ptr(pism::stressbalance::ShallowStressBalance)
%shared_ptr(pism::stressbalance::SSB_Modifier)

%shared_ptr(pism::stressbalance::StressBalance)

%ignore pism::array::Array2D< pism::stressbalance::DeviatoricStresses >::array;
%ignore pism::array::Array2D< pism::stressbalance::DeviatoricStresses >::add;
%shared_ptr(pism::array::Array2D< pism::stressbalance::DeviatoricStresses >)
%template(ArrayDeviatoricStresses) pism::array::Array2D<pism::stressbalance::DeviatoricStresses>;

%ignore pism::array::Array2D< pism::stressbalance::PrincipalStrainRates >::array;
%ignore pism::array::Array2D< pism::stressbalance::PrincipalStrainRates >::add;
%shared_ptr(pism::array::Array2D< pism::stressbalance::PrincipalStrainRates >)
%template(ArrayPrincipalStrainRates) pism::array::Array2D<pism::stressbalance::PrincipalStrainRates>;
%include "stressbalance/StressBalance.hh"

%shared_ptr(pism::stressbalance::ZeroSliding)
%shared_ptr(pism::stressbalance::PrescribedSliding)
%include "stressbalance/ShallowStressBalance.hh"

%shared_ptr(pism::stressbalance::SSA)
%include "stressbalance/ssa/SSA.hh"
%shared_ptr(pism::stressbalance::SSAFD)
%include "stressbalance/ssa/SSAFD.hh"
%shared_ptr(pism::stressbalance::SSAFEM)
%include "stressbalance/ssa/SSAFEM.hh"

%shared_ptr(pism::Diag<pism::stressbalance::SSA>)
%template(_SSADiag) pism::Diag<pism::stressbalance::SSA>;
%shared_ptr(pism::stressbalance::SSA_taud_mag)
%shared_ptr(pism::stressbalance::SSA_taud)
%shared_ptr(pism::stressbalance::SSA_calving_front_pressure_difference)
%include "stressbalance/ssa/SSA_diagnostics.hh"
%shared_ptr(pism::Diag<pism::stressbalance::SSAFD>)
%template(_SSAFDDiag) pism::Diag<pism::stressbalance::SSAFD>;
%shared_ptr(pism::stressbalance::SSAFD_nuH)
%include "stressbalance/ssa/SSAFD_diagnostics.hh"
