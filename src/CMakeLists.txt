add_library (pismnum
  num/cubature.c
  num/extrasGSL.cc)

set_source_files_properties (base/iMutil.cc
  PROPERTIES COMPILE_FLAGS -DPISM_REVISION='\"${Pism_REVISION}\"')

add_library (pismbase
  base/beddefLC.cc
  base/forcing.cc
  base/grid.cc
  base/iMIO.cc
  base/iMadaptive.cc
  base/iMbasal.cc
  base/iMbeddef.cc
  base/iMbootstrap.cc
  base/iMdefaults.cc
  base/iMforcing.cc
  base/iMgrainsize.cc
  base/iMinverse.cc
  base/iMmatlab.cc
  base/iMnames.cc
  base/iMoptions.cc
  base/iMpdd.cc
  base/iMreport.cc
  base/iMsia.cc
  base/iMssa.cc
  base/iMtemp.cc
  base/iMtests.cc
  base/iMutil.cc
  base/iMvelocity.cc
  base/iMviewers.cc
  base/iceModel.cc
  base/iceModelVec.cc
  base/iceModelVec2.cc
  base/iceModelVec3.cc
  base/materials.cc
  base/nc_util.cc
  base/pism_const.cc
  base/pism_signal.c)

add_library (pismverif
  verif/exactTestH.c
  verif/exactTestK.c
  verif/exactTestL.c
  verif/exactTestsABCDE.c
  verif/exactTestsFG.c
  verif/exactTestsIJ.c)

# All Pism executables must link against this
set (Pism_LIBS pismbase pismnum ${Pism_EXTERNAL_LIBS})

add_executable (pismr pismr.cc)
target_link_libraries (pismr ${Pism_LIBS})

add_executable (pismd pismd.cc
  eismint/iceROSSModel.cc)
target_link_libraries (pismd ${Pism_LIBS})

add_executable (pisms pisms.cc 
  eismint/iceEISModel.cc
  eismint/icePSTexModel.cc
  ismip/iceMISMIPModel.cc)
target_link_libraries (pisms ${Pism_LIBS})

add_executable (pismv pismv.cc
  verif/iceCompModel.cc
  verif/iCMthermo.cc
  verif/iceExactSSAModel.cc)
target_link_libraries (pismv pismverif ${Pism_LIBS})

add_executable (pgrn pgrn.cc
  eismint/iceGRNModel.cc)
target_link_libraries (pgrn ${Pism_LIBS})

install (TARGETS
  pismr pismd pisms pismv pgrn # executables
  pismverif pismbase pismnum # libraries
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)



if (Pism_BUILD_EXTRA_EXECS)
  # These don't require Pism at all
  set (EXTRA_EXECS simpleABCD simpleE simpleFG simpleH simpleI simpleJ simpleL gridL)
  foreach (EXEC ${EXTRA_EXECS})
    add_executable (${EXEC} verif/${EXEC}.c)
    target_link_libraries (${EXEC} pismverif ${GSL_LIBRARIES} ${PETSC_LIBRARIES}) # Link against PETSC because it's guaranteed to have BLAS
  endforeach (EXEC)
  # These use Pism in a trivial way and could be made to link without Pism
  add_executable (flowTable flowTable.cc)
  target_link_libraries (flowTable ${Pism_LIBS})
  add_executable (tryLCbd verif/tryLCbd.cc)
  target_link_libraries (tryLCbd ${Pism_LIBS})
  install (TARGETS
    ${EXTRA_EXECS} flowTable tryLCbd
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib)
endif (Pism_BUILD_EXTRA_EXECS)
