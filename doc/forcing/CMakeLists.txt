# Description of forcing options

# copy PISM-related macros, creating the revision tag
configure_file(../pism-macros.tex pism-macros.tex @ONLY)

# Copy the bibliography file (note the change in the file name)
configure_file (../ice-bib.bib ice-bib.bib COPYONLY)

pism_latex_document(pism_forcing "pism_forcing.tex"
  "pdd-movie.png;dT-movie.png;data-flow.pdf" "YES")
add_dependencies(pism_forcing pism_manual)
