import PISM

com = PISM.PETSc.COMM_WORLD
sys = PISM.UnitSystem("")
log = PISM.StringLogger(com, 1)

pism_config = PISM.config_from_options(com, log, sys)

config = PISM.ConfigJSON(sys)
config.init_from_string("""
{
 "constants" : {"ice" : {}, "fresh_water" : {},  "sea_water" : {}},
 "bootstrapping" : {"defaults" : {}},
 "calving" : {"eigen_calving" : {}, "thickness_calving" : {}, "float_kill" : {}},
 "enthalpy_converter" : {},
 "flow_law" : {
   "gpbld" : {},
   "Hooke" : {},
   "isothermal_Glen": {},
   "Paterson_Budd" : {}},
 "grid" : {},
 "hydrology" : {},
 "fracture_density" : {},
 "stress_balance" : {
   "ssa" : {"strength_extension" : {}, "fd" : {"lateral_drag" : {}}},
   "sia" : {}
 },
 "time" : {},
 "geometry" : {"update" : {}, "part_grid" : {}},
 "climate_forcing" : {},
 "inverse" : {"design" : {}, "ssa" : {}, "tikhonov" : {}},
 "bed_deformation" : {},
 "ocean" : {},
 "time_stepping" : {"skip" : {}},
 "energy" : {"basal_melt" : {}},
 "age" : {},
 "output" : {"runtime" : {}, "sizes" : {}},
 "run_info" : {},
 "surface" : {"force_to_thickness" : {}, "pdd" : {"fausto" : {}}},
 "atmosphere" : {"fausto_air_temp" : {}},
 "regional" : {},
 "basal_yield_stress" : {"constant" : {}, "mohr_coulomb" : {"topg_to_phi" : {}}},
 "basal_resistance" : {"plastic" : {}, "pseudo_plastic" : {}}
}""")
config.import_from(pism_config)

print(config.dump())
