# User's Manual:

# copy PISM-related macros, creating the revision tag
configure_file(../pism-macros.tex pism-macros.tex @ONLY)

# Copy the bibliography file (note the change in the file name)
configure_file (../ice_bib.bib ice-bib.bib COPYONLY)

# The list of all the source files, with section number:
set (latex_sources
  manual.tex                   # preamble and structure
  getting-started.tex          # 2
  highlevelview.tex            # 3
  init-boot.tex                # 4
  modeling-computational.tex   # 5
  modeling-dynamics.tex        # 6
  modeling-marine.tex          # 7
  practical-usage.tex          # 8
  verification.tex             # 9
  simplified-geometry.tex      # 10
  ross.tex                     # 11
  storglaciaren.tex            # 12
  jako.tex                     # 13
  )

# Find all the figures in figs/
foreach (format jpg png pdf)
  file (GLOB ${format}
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    figs/*.${format})
endforeach()
set (figures ${jpg} ${png} ${pdf})

pism_latex_document(manual "${latex_sources}" "${figures}")
