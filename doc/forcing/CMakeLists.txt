# Description of forcing options

# copy PISM-related macros, creating the revision tag
configure_file(../pism-macros.tex pism-macros.tex @ONLY)

# Copy the bibliography file (note the change in the file name)
configure_file (../ice_bib.bib ice-bib.bib COPYONLY)

pism_latex_document(forcing "forcing.tex" "pdd-movie.png;dT-movie.png;data-flow.pdf")
