include_directories (base coupler eismint ismip num verif udunits)

add_library (pismudunits
  udunits/udalloc.c
  udunits/utlib.c
  udunits/utparse.c
  udunits/utscan.c
)
set_source_files_properties (udunits/utlib.c
  PROPERTIES COMPILE_FLAGS
  "-DUT_INSTALL_PATH='\"${CMAKE_INSTALL_PREFIX}/lib/pism/pismudunits.dat\"' -DUT_SOURCE_PATH='\"${CMAKE_SOURCE_DIR}/src/udunits/pismudunits.dat\"'")

add_library (pismbase
  base/LocalInterpCtx.cc
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
  base/materials.cc
  base/nc_util.cc
  base/pism_const.cc
  base/pism_revision.cc
  base/pism_signal.c
  base/ssaJed/iceShelfExtension.cc
  base/ssaJed/ssa.cc
  base/ssaJed/ssafe.cc
  num/cubature.c
  num/extrasGSL.cc
  coupler/forcing.cc
  coupler/pPDDcoupler.cc
  coupler/pccoupler.cc
  )
target_link_libraries (pismbase pismudunits ${Pism_EXTERNAL_LIBS})

set_source_files_properties (base/pism_revision.cc
  PROPERTIES COMPILE_FLAGS -DPISM_REVISION='\"${Pism_REVISION_TAG}\"')

add_library (pismverif
  verif/exactTestH.c
  verif/exactTestK.c
  verif/exactTestL.c
  verif/exactTestM.c
  verif/exactTestsABCDE.c
  verif/exactTestsFG.c
  verif/exactTestsIJ.c
  verif/iCMthermo.cc
  verif/iceCalvBCModel.cc
  verif/iceCompModel.cc
  verif/iceExactSSAModel.cc)
target_link_libraries (pismverif pismbase)


# Main executables:
add_executable (pismr pismr.cc)

add_executable (pisms pisms.cc 
  eismint/iceEISModel.cc
  eismint/icePSTexModel.cc
  ismip/iceMISMIPModel.cc)

add_executable (pgrn eismint/pgrn.cc eismint/iceGRNModel.cc)

add_executable (pismd pismd.cc eismint/iceROSSModel.cc)

# All of the following are linked against pismbase
foreach (EXEC pismr pisms pgrn pismd)
  target_link_libraries (${EXEC} pismbase)
endforeach (EXEC)

add_executable (pismv pismv.cc
  verif/iceCompModel.cc
  verif/iCMthermo.cc
  verif/iceExactSSAModel.cc
  verif/iceCalvBCModel.cc)
target_link_libraries (pismv pismverif pismbase)

add_executable (pcctest pcctest.cc)
target_link_libraries (pcctest pismbase)

install (TARGETS
  pismr pisms pismv pgrn pismd pcctest ## executables
  pismbase pismverif pismudunits # libraries
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib/pism
  ARCHIVE DESTINATION lib)

install (FILES
  "udunits/pismudunits.dat"
  DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pism")

if (Pism_BUILD_EXTRA_EXECS)
  set (SIMPLE_EXECS simpleABCD simpleE simpleFG simpleH simpleI simpleJ simpleL)
  foreach (EXEC ${SIMPLE_EXECS})
    add_executable (${EXEC} verif/${EXEC}.c)
    target_link_libraries (${EXEC} pismverif)
  endforeach (EXEC)

  add_executable (tryLCbd verif/tryLCbd.cc base/beddefLC.cc base/materials.cc base/pism_const.cc)
  target_link_libraries (tryLCbd pismbase)

  add_executable (flowTable flowTable.cc base/materials.cc)
  target_link_libraries (flowTable pismbase)

  install (TARGETS
    ${SIMPLE_EXECS} tryLCbd flowTable
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib/pism)
endif (Pism_BUILD_EXTRA_EXECS)
