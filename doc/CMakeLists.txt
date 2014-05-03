# Find latex and create doc targets:
project(PISM_Documentation "NONE")
cmake_minimum_required (VERSION 2.8.0)

find_package(LATEX)
find_package(Doxygen)

if (NOT Pism_SOURCE_DIR)
  # This assumes that the current directory ("doc") is right in the
  # PISM source code directory, i.e. pism/doc.
  get_filename_component(Pism_SOURCE_DIR "${PROJECT_SOURCE_DIR}/.." REALPATH)
  message(STATUS "Setting Pism_SOURCE_DIR to ${Pism_SOURCE_DIR}")

  add_subdirectory(${Pism_SOURCE_DIR}/site-packages ${CMAKE_BINARY_DIR}/site-packages)
endif()

list (APPEND CMAKE_MODULE_PATH "${Pism_SOURCE_DIR}/CMake")

# a macro creating a custom command to generate a PDF document using "general" and "command-line option" indices.
macro(pism_latex_document name src figs run_bibtex)
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

  if (${run_bibtex})
    add_custom_command (OUTPUT ${name}.pdf
      COMMAND ${PDFLATEX_COMPILER} ${name}.tex
      COMMAND ${BIBTEX_COMPILER} ${name}
      COMMAND ${MAKEINDEX_COMPILER} -l ${name}
      COMMAND ${MAKEINDEX_COMPILER} -l -o ${name}.ond ${name}.odx
      COMMAND ${PDFLATEX_COMPILER} ${name}.tex
      COMMAND ${PDFLATEX_COMPILER} ${name}.tex
      COMMAND ln -sf ${CMAKE_CURRENT_BINARY_DIR}/${name}.pdf ${CMAKE_BINARY_DIR}/
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      DEPENDS ${src} ${figs})
  else()
    add_custom_command (OUTPUT ${name}.pdf
      COMMAND ${PDFLATEX_COMPILER} ${name}.tex
      COMMAND ${MAKEINDEX_COMPILER} -l ${name}
      COMMAND ${MAKEINDEX_COMPILER} -l -o ${name}.ond ${name}.odx
      COMMAND ${PDFLATEX_COMPILER} ${name}.tex
      COMMAND ${PDFLATEX_COMPILER} ${name}.tex
      COMMAND ln -sf ${CMAKE_CURRENT_BINARY_DIR}/${name}.pdf ${CMAKE_BINARY_DIR}/
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      DEPENDS ${src} ${figs})
  endif()

  add_custom_command (OUTPUT ${name}_draft.pdf
    COMMAND ${PDFLATEX_COMPILER} ${name}.tex
    COMMAND mv ${CMAKE_CURRENT_BINARY_DIR}/${name}.pdf ${CMAKE_CURRENT_BINARY_DIR}/${name}_draft.pdf
    COMMAND ln -sf ${CMAKE_CURRENT_BINARY_DIR}/${name}_draft.pdf ${CMAKE_BINARY_DIR}/
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

  # Installation manual:
  add_subdirectory(installation)

  # User's Manual:
  add_subdirectory(userman)

  # Climate forcing components manual:
  add_subdirectory(forcing)

  # Source code browser:
  add_subdirectory (browser EXCLUDE_FROM_ALL)

  # PISM's Python bindings docs:
  add_subdirectory (pismpython EXCLUDE_FROM_ALL)

  if (Pism_BUILD_PDFS AND Pism_DOC_DIR)
    install (FILES
      "${CMAKE_CURRENT_BINARY_DIR}/userman/manual.pdf"
      "${CMAKE_CURRENT_BINARY_DIR}/forcing/forcing.pdf"
      "${CMAKE_CURRENT_BINARY_DIR}/installation/installation.pdf"
      DESTINATION ${Pism_DOC_DIR})

  endif ()                      # if Pism_BUILD_PDFS
endif ()                        # if PDFLATEX_COMPILER
