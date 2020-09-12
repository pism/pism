
pism_class(pism::stressbalance::Blatter,
           "pism/stressbalance/blatter/Blatter.hh")

pism_class(pism::stressbalance::BlatterTest1,
           "pism/stressbalance/blatter/verification/BlatterTest1.hh")

pism_class(pism::stressbalance::BlatterISMIPHOM,
           "pism/stressbalance/blatter/ismip-hom/BlatterISMIPHOM.hh")

/* BlatterMod has to be wrapped after Blatter*/
pism_class(pism::stressbalance::BlatterMod,
           "pism/stressbalance/blatter/BlatterMod.hh")
