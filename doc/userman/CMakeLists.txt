# User's Manual:

# Create the revision stamp:
file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/revision.tex "${Pism_REVISION_TAG}")

# Copy the bibliography file (note the change in the file name)
configure_file (../ice_bib.bib ice-bib.bib COPYONLY)

# Find all the .tex files in userman and copy them to the build directory

set (latex_sources
  eismint-greenland.tex
  eismint-ross.tex
  getting-started.tex
  hard-choices.tex
  manual.tex
  modeling-choices.tex
  pism-pik.tex
  practical-usage.tex
  simplified-geometry.tex
  storglaciaren.tex
  verification.tex
  storglaciaren.tex
  )

foreach (filename ${latex_sources})
  configure_file(${filename} ${filename} COPYONLY)
endforeach()

# Find all the figures and copy them to the build directory
foreach (format jpg png pdf)
  file (GLOB ${format}
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/figs
    figs/*.${format})
endforeach()
foreach (filename ${jpg} ${png} ${pdf})
  list(APPEND figures figs/${filename})
  configure_file(figs/${filename} figs/${filename} COPYONLY)
endforeach()

# Commands to build the manual:
add_custom_command (OUTPUT manual.pdf
  COMMAND ${PDFLATEX_COMPILER} manual.tex
  COMMAND ${BIBTEX_COMPILER} manual
  COMMAND ${MAKEINDEX_COMPILER} -l manual
  COMMAND ${MAKEINDEX_COMPILER} -l -o manual.ond manual.odx
  COMMAND ${PDFLATEX_COMPILER} manual.tex
  COMMAND ${PDFLATEX_COMPILER} manual.tex
  COMMAND ln -sf ${CMAKE_CURRENT_BINARY_DIR}/manual.pdf ${PROJECT_BINARY_DIR}/
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${latex_sources} ${figures})

add_custom_command (OUTPUT manual_draft.pdf
  COMMAND ${PDFLATEX_COMPILER} manual.tex
  COMMAND mv ${CMAKE_CURRENT_BINARY_DIR}/manual.pdf ${CMAKE_CURRENT_BINARY_DIR}/manual_draft.pdf
  COMMAND ln -sf ${CMAKE_CURRENT_BINARY_DIR}/manual.pdf ${PROJECT_BINARY_DIR}/
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${latex_sources} ${figures})

add_custom_target(userman ALL
  DEPENDS manual.pdf)

add_custom_target(manual_draft
  DEPENDS manual_draft.pdf)
