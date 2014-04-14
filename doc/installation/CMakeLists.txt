# Installation Manual:

# copy PISM-related macros, creating the revision tag
configure_file(../pism-macros.tex pism-macros.tex @ONLY)

# The list of all the source files for the User's Manual, with section number:
set (latex_sources
  installation.tex
  )

pism_latex_document(installation "${latex_sources}"
  ""                            # no figures
  "NO"                          # no bibtex
  )
