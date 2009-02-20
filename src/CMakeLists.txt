add_library (pismnum
  num/cubature.c
  num/extrasGSL.cc)

set_source_files_properties (base/pism_revision.cc
  PROPERTIES COMPILE_FLAGS -DPISM_REVISION='\"${Pism_REVISION}\"')
include_directories (base udunits)

add_library (pismbase
  base/pism_revision.cc
  base/beddefLC.cc
  base/grid.cc
  base/iMIO.cc
  base/iMadaptive.cc
  base/iMbasal.cc
  base/iMbeddef.cc
  base/iMbootstrap.cc
  base/iMdefaults.cc
  base/iMgeometry.cc
  base/iMgrainsize.cc
  base/iMinverse.cc
  base/iMinverseMat.cc
  base/iMmatlab.cc
  base/iMnames.cc
  base/iMoptions.cc
  base/iMreport.cc
  base/iMsia.cc
  base/iMssa.cc
  base/iMssaSNES.cc
  base/iMtemp.cc
  base/iMtests.cc
  base/iMutil.cc
  base/iMvelocity.cc
  base/iMviewers.cc
  base/iceModel.cc
  base/iceModelVec.cc
  base/iceModelVec2.cc
  base/iceModelVec3.cc
  base/LocalInterpCtx.cc
  base/materials.cc
  base/nc_util.cc
  base/pism_const.cc
  base/pism_signal.c
  coupler/forcing.cc
  coupler/pccoupler.cc
  coupler/pPDDcoupler.cc
  )

add_library (pismverif
  verif/exactTestH.c
  verif/exactTestK.c
  verif/exactTestL.c
  verif/exactTestM.c
  verif/exactTestsABCDE.c
  verif/exactTestsFG.c
  verif/exactTestsIJ.c)

add_library (pismudunits
  udunits/udalloc.c
  udunits/utlib.c
  udunits/utparse.c
  udunits/utscan.c)
set_target_properties (pismudunits
  PROPERTIES COMPILE_DEFINITIONS UT_DEFAULT_PATH="${CMAKE_INSTALL_PREFIX}/lib/udunits.dat")

# All Pism executables must link against this
set (Pism_LIBS pismbase pismnum pismudunits ${Pism_EXTERNAL_LIBS})

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
  verif/iceCalvBCModel.cc
  verif/iceExactSSAModel.cc)
target_link_libraries (pismv pismverif ${Pism_LIBS})

add_executable (pgrn eismint/pgrn.cc
  eismint/iceGRNModel.cc)
target_link_libraries (pgrn ${Pism_LIBS})

add_executable (pcctest pcctest.cc)
target_link_libraries (pcctest ${Pism_LIBS})

install (TARGETS
  pismr pisms pismv pismd pgrn ## executables
  pismverif pismbase pismnum pismudunits # libraries
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install (FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/udunits/udunits.dat"
  DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")


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
