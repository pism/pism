# Find latex and create doc targets:
find_package(LATEX)
if (LATEX_COMPILER)
  set (PISM_DOCS userman cheatsheet installation)

  # install (FILES
  #   ${CMAKE_CURRENT_SOURCE_DIR}/installation.pdf
  #   ${CMAKE_CURRENT_SOURCE_DIR}/cheatsheet.pdf
  #   ${CMAKE_CURRENT_SOURCE_DIR}/userman/manual.pdf
  #   DESTINATION doc/)
endif (LATEX_COMPILER)

find_package(Doxygen)
if (LATEX_COMPILER AND DOXYGEN_FOUND)
  list (APPEND PISM_DOCS browser browser.tgz)

  # install (FILES
  #   ${CMAKE_CURRENT_SOURCE_DIR}/browser/browser.tgz
  #   DESTINATION doc/browser/)
endif ()

foreach (DOC ${PISM_DOCS})
  add_custom_target (${DOC}
    COMMAND ${CMAKE_MAKE_PROGRAM} -C ${CMAKE_CURRENT_SOURCE_DIR} ${DOC} &> /dev/null
    COMMENT "Building ${DOC}...")
endforeach (DOC)

add_custom_target (docs
  DEPENDS ${PISM_DOCS})
