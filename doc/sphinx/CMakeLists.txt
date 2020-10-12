find_program(SPHINX_EXECUTABLE sphinx-build)
mark_as_advanced(SPHINX_EXECUTABLE)

# HTML documentation requires Sphinx
if (SPHINX_EXECUTABLE)
  find_program (NCGEN_PROGRAM "ncgen" REQUIRED)
  mark_as_advanced(NCGEN_PROGRAM)

  add_custom_command (OUTPUT pism_config.nc
    COMMAND ${NCGEN_PROGRAM} -o pism_config.nc ${Pism_SOURCE_DIR}/src/pism_config.cdl
    DEPENDS ${Pism_SOURCE_DIR}/src/pism_config.cdl)

  add_custom_target (manual_html ALL
    COMMAND ${SPHINX_EXECUTABLE} -W -n -b html -D pism_config_file=${CMAKE_CURRENT_BINARY_DIR}/pism_config.nc ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/html
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS pism_config.nc)

  add_custom_target (manual_linkcheck
    COMMAND ${SPHINX_EXECUTABLE} -n -b linkcheck -D pism_config_file=${CMAKE_CURRENT_BINARY_DIR}/pism_config.nc ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS pism_config.nc)

  find_program(SPHINX_AUTOBUILD_EXECUTABLE sphinx-autobuild)
  mark_as_advanced(SPHINX_AUTOBUILD_EXECUTABLE)

  if (SPHINX_AUTOBUILD_EXECUTABLE)
    add_custom_target (manual_autobuild
      COMMAND ${SPHINX_AUTOBUILD_EXECUTABLE} --port 8080 -B ${CMAKE_CURRENT_SOURCE_DIR} -D pism_config_file=${CMAKE_CURRENT_BINARY_DIR}/pism_config.nc ${CMAKE_CURRENT_BINARY_DIR}/manual-tmp
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDS pism_config.nc)
  endif()

  if (Pism_DOC_DIR)
    install (DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html"
      DESTINATION ${Pism_DOC_DIR}
      PATTERN ".doctrees" EXCLUDE
      PATTERN ".buildinfo" EXCLUDE
      PATTERN "objects.inv" EXCLUDE)
  endif ()                      # if Pism_DOC_DIR

endif()

# PDF documentation requires Sphinx and pdflatex (and latexmk)
if (SPHINX_EXECUTABLE AND PDFLATEX_COMPILER)

  find_program (MAKE_EXECUTABLE NAMES make gmake)

  add_custom_command (OUTPUT pism_manual.pdf
    COMMAND ${SPHINX_EXECUTABLE} -E -b latex -t latex -D pism_config_file=${CMAKE_CURRENT_BINARY_DIR}/pism_config.nc ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${MAKE_EXECUTABLE} -C ${CMAKE_CURRENT_BINARY_DIR} all
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_custom_target(manual_pdf ALL
    DEPENDS pism_manual.pdf
    DEPENDS pism_config.nc)

  if (Pism_DOC_DIR)
    install (FILES ${CMAKE_CURRENT_BINARY_DIR}/pism_manual.pdf
      DESTINATION ${Pism_DOC_DIR})
  endif ()                      # if Pism_DOC_DIR

endif()
