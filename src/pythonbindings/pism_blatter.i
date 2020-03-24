%{
#include "stressbalance/blatter/BlatterStressBalance.hh"
%}

%shared_ptr(pism::stressbalance::BlatterStressBalance)
%include "stressbalance/blatter/BlatterStressBalance.hh"

pism_class(pism::stressbalance::Poisson3, "pism/stressbalance/blatter/Poisson3.hh")
