
# Verification tests.
add_library (verif OBJECT
  TemperatureModel_Verification.cc
  iCMthermo.cc
  iceCompModel.cc
  PSVerification.cc
  BTU_Verification.cc
  tests/exactTestH.c
  tests/exactTestK.c
  tests/exactTestO.c
  tests/exactTestL.cc
  tests/exactTestsABCD.c
  tests/exactTestsFG.cc
  tests/exactTestsIJ.c
  tests/exactTestM.c
  tests/exactTestN.c
  tests/exactTestP.cc
)
