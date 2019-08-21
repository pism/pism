find_package(Sphinx)

if (SPHINX_EXECUTABLE)

  add_custom_target (pismpython_docs ALL
    COMMAND ${CMAKE_COMMAND} -E env PISM_BUILD_DIR="${CMAKE_BINARY_DIR}" ${SPHINX_EXECUTABLE} -b html ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/html
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  add_custom_target (pismpython_linkcheck
    COMMAND ${CMAKE_COMMAND} -E env PISM_BUILD_DIR="${CMAKE_BINARY_DIR}" ${SPHINX_EXECUTABLE} -b linkcheck ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  if (Pism_DOC_DIR)
    install (DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html"
      DESTINATION ${Pism_DOC_DIR}
      PATTERN ".doctrees" EXCLUDE
      PATTERN ".buildinfo" EXCLUDE
      PATTERN "objects.inv" EXCLUDE)
  endif ()                      # if Pism_DOC_DIR

endif()
