# Description of forcing options

# Create the revision stamp:
file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/revision.tex "${Pism_REVISION_TAG}")

# Copy the bibliography file (note the change in the file name)
configure_file (../ice_bib.bib ice-bib.bib COPYONLY)

pism_latex_document(forcing "forcing.tex" "")
