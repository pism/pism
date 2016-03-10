# Installation Manual:

# copy PISM-related macros, creating the revision tag
configure_file(../pism-macros.tex pism-macros.tex @ONLY)

# The list of all the source files for this manual
set (latex_sources
  pism_installation.tex
  )

pism_latex_document(pism_installation "${latex_sources}"
  ""                            # no figures
  "NO"                          # no bibtex
  )
