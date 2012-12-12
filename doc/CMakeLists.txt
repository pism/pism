# Find latex and create doc targets:
find_package(LATEX)
find_package(Doxygen)

# a macro creating a custom command generating a PDF document using general and command-line option indices.
macro(pism_latex_document name src figs)
  # copy files into the build directory:
  foreach (filename ${src} ${figs})
    if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${filename})
      configure_file(${filename} ${filename} COPYONLY)
    else ()
      set(missing "${missing} ${filename}")
    endif()
  endforeach()

  if (missing)
    message(STATUS "Will try to generate ${missing}")
  endif()

  add_custom_command (OUTPUT ${name}.pdf
    COMMAND ${PDFLATEX_COMPILER} ${name}.tex
    COMMAND ${BIBTEX_COMPILER} ${name}
    COMMAND ${MAKEINDEX_COMPILER} -l ${name}
    COMMAND ${MAKEINDEX_COMPILER} -l -o ${name}.ond ${name}.odx
    COMMAND ${PDFLATEX_COMPILER} ${name}.tex
    COMMAND ${PDFLATEX_COMPILER} ${name}.tex
    COMMAND ln -sf ${CMAKE_CURRENT_BINARY_DIR}/${name}.pdf ${PROJECT_BINARY_DIR}/
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${src} ${figs})

  add_custom_command (OUTPUT ${name}_draft.pdf
    COMMAND ${PDFLATEX_COMPILER} ${name}.tex
    COMMAND mv ${CMAKE_CURRENT_BINARY_DIR}/${name}.pdf ${CMAKE_CURRENT_BINARY_DIR}/${name}_draft.pdf
    COMMAND ln -sf ${CMAKE_CURRENT_BINARY_DIR}/${name}_draft.pdf ${PROJECT_BINARY_DIR}/
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${src} ${figs})

  add_custom_target(${name} ALL
    DEPENDS ${name}.pdf)

  add_custom_target(${name}_draft
    DEPENDS ${name}_draft.pdf)

endmacro()

if (PDFLATEX_COMPILER)
  # Only offer to build docs if LaTeX is found.
  option (Pism_BUILD_PDFS "Build PISM PDF documentation" OFF)

  # copy PISM-related macros, creating the revision tag
  configure_file(pism-macros.tex pism-macros.tex @ONLY)

  # Installation manual:
  foreach (NAME installation)
    configure_file (${NAME}.tex ${NAME}.tex COPYONLY)

    add_custom_command (OUTPUT ${NAME}.pdf
      COMMAND ${PDFLATEX_COMPILER} ${NAME}.tex
      COMMAND ${PDFLATEX_COMPILER} ${NAME}.tex
      COMMAND ln -sf ${CMAKE_CURRENT_BINARY_DIR}/${NAME}.pdf ${PROJECT_BINARY_DIR}/
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      DEPENDS ${NAME}.tex)

    # A custom target to be able to build PDFs by "make installation", etc
    add_custom_target(${NAME} ALL
      DEPENDS ${NAME}.pdf)
  endforeach ()

  add_subdirectory(userman)

  add_subdirectory(forcing)

  if (Pism_BUILD_PDFS)
    install (FILES
      "${CMAKE_CURRENT_BINARY_DIR}/userman/manual.pdf"
      "${CMAKE_CURRENT_BINARY_DIR}/forcing/forcing.pdf"
      "${CMAKE_CURRENT_BINARY_DIR}/installation.pdf"
      DESTINATION ${Pism_DOC_DIR})

  endif ()                      # if Pism_BUILD_PDFS
endif ()                        # if PDFLATEX_COMPILER

add_subdirectory(manpages)