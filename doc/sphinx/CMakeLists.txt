find_package(Sphinx)

# HTML documentation requires Sphinx
if (SPHINX_EXECUTABLE)

  add_custom_target (manual_html ALL
    COMMAND ${SPHINX_EXECUTABLE} -n -b html ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/html
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  add_custom_target (manual_linkcheck
    COMMAND ${SPHINX_EXECUTABLE} -n -b linkcheck ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

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
    COMMAND ${SPHINX_EXECUTABLE} -b latex -t latex ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${MAKE_EXECUTABLE} all
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )

  add_custom_target(manual_pdf ALL
    DEPENDS pism_manual.pdf)

  if (Pism_DOC_DIR)
    install (FILES ${CMAKE_CURRENT_BINARY_DIR}/pism_manual.pdf
      DESTINATION ${Pism_DOC_DIR})
  endif ()                      # if Pism_DOC_DIR

endif()
