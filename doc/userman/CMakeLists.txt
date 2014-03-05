# User's Manual:

# copy PISM-related macros, creating the revision tag
configure_file(../pism-macros.tex pism-macros.tex @ONLY)

# Copy the bibliography file (note the change in the file name)
configure_file (../ice-bib.bib ice-bib.bib COPYONLY)

# The list of all the source files for the User's Manual, with section number:
set (latex_sources
  manual.tex                   # preamble and structure
  getting-started.tex          # 2
  highlevelview.tex            # 3
  init-boot.tex                # 4
  modeling-computational.tex   # 5
  modeling-dynamics.tex        # 6
  modeling-subglacier.tex      # 7
  modeling-marine.tex          # 8
  practical-usage.tex          # 9
  verification.tex             # 10
  simplified-geometry.tex      # 11
  validation.tex               # 12
  storglaciaren.tex            # 13
  jako.tex                     # 14
  )

# Find all the figures in figs/
foreach (format jpg png pdf)
  file (GLOB ${format}
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    figs/*.${format})
endforeach()
set (figures ${jpg} ${png} ${pdf})

pism_latex_document(manual "${latex_sources}" "${figures}")
