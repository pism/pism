find_program(SPHINX_EXECUTABLE sphinx-build)
mark_as_advanced(SPHINX_EXECUTABLE)

# HTML documentation requires Sphinx
if (SPHINX_EXECUTABLE)
  find_program (NCGEN_PROGRAM "ncgen" REQUIRED)
  mark_as_advanced(NCGEN_PROGRAM)

  # copy conf.py, changing some things along the way:
  configure_file(conf.py conf.py)

  add_custom_command (OUTPUT pism_config.nc
    COMMAND ${NCGEN_PROGRAM} -o pism_config.nc ${Pism_SOURCE_DIR}/src/pism_config.cdl
    DEPENDS ${Pism_SOURCE_DIR}/src/pism_config.cdl)

  add_custom_target (manual_html ALL
    COMMAND ${SPHINX_EXECUTABLE} -c ${CMAKE_CURRENT_BINARY_DIR} -W -n -b html ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/manual_html
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS pism_config.nc)

  add_custom_target (manual_linkcheck
    COMMAND ${SPHINX_EXECUTABLE} -c ${CMAKE_CURRENT_BINARY_DIR} -n -b linkcheck ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDS pism_config.nc)

  find_program(SPHINX_AUTOBUILD_EXECUTABLE sphinx-autobuild)
  mark_as_advanced(SPHINX_AUTOBUILD_EXECUTABLE)

  if (SPHINX_AUTOBUILD_EXECUTABLE)
    add_custom_target (manual_autobuild
      COMMAND ${SPHINX_AUTOBUILD_EXECUTABLE} --open-browser --port 8080 -c ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/manual-tmp
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDS pism_config.nc)
  endif()

  install (DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/manual_html"
    DESTINATION ${CMAKE_INSTALL_DOCDIR}
    PATTERN ".doctrees" EXCLUDE
    PATTERN ".buildinfo" EXCLUDE
    PATTERN "objects.inv" EXCLUDE)

endif()

# PDF documentation requires Sphinx and pdflatex (and latexmk)
if (SPHINX_EXECUTABLE AND PDFLATEX_COMPILER)

  find_program (MAKE_EXECUTABLE NAMES make gmake)

  add_custom_command (OUTPUT pism_manual.pdf
    COMMAND ${SPHINX_EXECUTABLE} -T -W -E -b latex -t latex -c ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${MAKE_EXECUTABLE} -C ${CMAKE_CURRENT_BINARY_DIR} all
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  add_custom_target(manual_pdf ALL
    DEPENDS pism_config.nc
    DEPENDS pism_manual.pdf
  )

  install (FILES ${CMAKE_CURRENT_BINARY_DIR}/pism_manual.pdf
    DESTINATION ${CMAKE_INSTALL_DOCDIR})

endif()
