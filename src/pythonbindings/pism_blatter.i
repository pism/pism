
pism_class(pism::stressbalance::Blatter,
           "pism/stressbalance/blatter/Blatter.hh")

pism_class(pism::stressbalance::BlatterTestXY,
           "pism/stressbalance/blatter/verification/BlatterTestXY.hh")

pism_class(pism::stressbalance::BlatterTestXZ,
           "pism/stressbalance/blatter/verification/BlatterTestXZ.hh")

pism_class(pism::stressbalance::BlatterTestCFBC,
           "pism/stressbalance/blatter/verification/BlatterTestCFBC.hh")

pism_class(pism::stressbalance::BlatterTestHalfar,
           "pism/stressbalance/blatter/verification/BlatterTestHalfar.hh")

pism_class(pism::stressbalance::BlatterTestvanderVeen,
           "pism/stressbalance/blatter/verification/BlatterTestvanderVeen.hh")

pism_class(pism::stressbalance::BlatterISMIPHOM,
           "pism/stressbalance/blatter/ismip-hom/BlatterISMIPHOM.hh")

/* BlatterMod has to be wrapped after Blatter*/
pism_class(pism::stressbalance::BlatterMod,
           "pism/stressbalance/blatter/BlatterMod.hh")

%{
#include "pism/stressbalance/blatter/verification/manufactured_solutions.hh"
%}
%include "pism/stressbalance/blatter/verification/manufactured_solutions.hh"
