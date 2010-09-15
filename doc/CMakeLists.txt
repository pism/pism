# Find latex and create doc targets:
find_package(LATEX)
if (LATEX_COMPILER)
  set (PISM_DOCS userman cheatsheet installation)
endif (LATEX_COMPILER)

find_package(Doxygen)
if (LATEX_COMPILER AND DOXYGEN_FOUND)

  list (APPEND PISM_DOCS browser)

endif (LATEX_COMPILER AND DOXYGEN_FOUND)

foreach (DOC ${PISM_DOCS})
  add_custom_target (${DOC}
    COMMAND ${CMAKE_MAKE_PROGRAM} -C ${CMAKE_CURRENT_SOURCE_DIR} ${DOC} &> /dev/null
    COMMENT "Building ${DOC}..."
    )
endforeach (DOC)

add_custom_target (docs
  DEPENDS ${PISM_DOCS}
  )
