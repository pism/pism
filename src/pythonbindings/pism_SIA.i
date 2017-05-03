%{
#include "base/stressbalance/sia/SIAFD.hh"
#include "base/stressbalance/sia/PISMBedSmoother.hh"
%}

%shared_ptr(pism::stressbalance::SSB_Modifier)
%shared_ptr(pism::stressbalance::ConstantInColumn)
%include "base/stressbalance/SSB_Modifier.hh"
%shared_ptr(pism::stressbalance::SIAFD)
%include "base/stressbalance/sia/SIAFD.hh"

%shared_ptr(pism::stressbalance::BedSmoother)
%include "base/stressbalance/sia/PISMBedSmoother.hh"
