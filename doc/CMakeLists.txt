# Find latex and create doc targets:
find_package(LATEX)
find_package(Doxygen)

if (PDFLATEX_COMPILER)
  # Only offer to build docs if LaTeX is found.
  option (Pism_BUILD_DOCS "Build PISM PDF documentation" OFF)

  file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/revision.tex "${Pism_REVISION_TAG}")

  # Installation manual and the cheatsheet:
  foreach (NAME installation cheatsheet)
    configure_file (${NAME}.tex ${NAME}.tex COPYONLY)

    add_custom_command (OUTPUT ${NAME}.pdf
      COMMAND ${PDFLATEX_COMPILER} ${NAME}.tex
      COMMAND ${PDFLATEX_COMPILER} ${NAME}.tex
      COMMAND ln -sf ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.pdf ${PROJECT_BINARY_DIR}/
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    # A custom target to be able to build PDFs by "make installation", etc
    add_custom_target(${NAME} ALL
      DEPENDS ${NAME}.pdf)
  endforeach ()

  add_subdirectory(userman)

  if (Pism_BUILD_DOCS)
    install (FILES
      "${CMAKE_CURRENT_BINARY_DIR}/userman/manual.pdf"
      "${CMAKE_CURRENT_BINARY_DIR}/installation.pdf"
      "${CMAKE_CURRENT_BINARY_DIR}/cheatsheet.pdf"
      DESTINATION ${Pism_DOC_DIR})

  endif ()                      # if Pism_BUILD_DOCS
endif ()                        # if PDFLATEX_COMPILER

