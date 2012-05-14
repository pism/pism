# User's Manual:

# copy PISM-related macros, creating the revision tag
configure_file(../pism-macros.tex pism-macros.tex @ONLY)

# Copy the bibliography file (note the change in the file name)
configure_file (../ice_bib.bib ice-bib.bib COPYONLY)

# The list of all the source files:
set (latex_sources
  ross.tex
  getting-started.tex
  manual.tex
  modeling-choices.tex
  pism-pik.tex
  practical-usage.tex
  simplified-geometry.tex
  storglaciaren.tex
  verification.tex
  storglaciaren.tex
  )

# Find all the figures in figs/
foreach (format jpg png pdf)
  file (GLOB ${format}
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    figs/*.${format})
endforeach()
set (figures ${jpg} ${png} ${pdf})

pism_latex_document(manual "${latex_sources}" "${figures}")
