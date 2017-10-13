find_package(Sphinx)
find_package(Python)

if (SPHINX_EXECUTABLE)

  find_program (MAKE_EXECUTABLE NAMES make gmake)

  add_custom_target (manual_html
    COMMAND ${SPHINX_EXECUTABLE} -b html ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )

  if (PDFLATEX_COMPILER)
    add_custom_target (manual_pdf
      COMMAND ${SPHINX_EXECUTABLE} -b latex -t latex ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/pdf
      COMMAND ${MAKE_EXECUTABLE} -C ${CMAKE_BINARY_DIR}/pdf all
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
  endif()

endif()
