%{
#include "stressbalance/sia/SIAFD.hh"
#include "stressbalance/sia/BedSmoother.hh"
%}

%shared_ptr(pism::stressbalance::ConstantInColumn)
%include "stressbalance/SSB_Modifier.hh"
%shared_ptr(pism::stressbalance::SIAFD)
%include "stressbalance/sia/SIAFD.hh"

%shared_ptr(pism::stressbalance::BedSmoother)
%include "stressbalance/sia/BedSmoother.hh"
