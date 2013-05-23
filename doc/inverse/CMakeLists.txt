# Description of inverse bolt-on design

# copy PISM-related macros, creating the revision tag
configure_file(../pism-macros.tex pism-macros.tex @ONLY)

# Copy the bibliography file (note the change in the file name)
configure_file (../ice_bib.bib ice-bib.bib COPYONLY)

pism_latex_document(ip_design "ip_design.tex" "" )
add_dependencies(ip_design manual)
