# Find latex and create doc targets:
find_package(LATEX)
find_package(Doxygen)

# Only offer to build docs if LaTeX is found.
if (PDFLATEX_COMPILER)
  option (Pism_BUILD_DOCS "Build PISM documentation" OFF)
endif()

if (PDFLATEX_COMPILER AND Pism_BUILD_DOCS)
  file (WRITE ${CMAKE_CURRENT_BINARY_DIR}/revision.tex "${Pism_REVISION}")

  # Installation manual and the cheatsheet:
  foreach (NAME installation cheatsheet)

    configure_file (${NAME}.tex ${NAME}.tex COPYONLY)

    add_custom_command (OUTPUT ${NAME}.pdf
      COMMAND ${PDFLATEX_COMPILER} ${NAME}.tex
      COMMAND ${PDFLATEX_COMPILER} ${NAME}.tex
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

    install (FILES "${CMAKE_CURRENT_BINARY_DIR}/${NAME}.pdf"
      DESTINATION share/doc/pism)
  endforeach ()

  # User's Manual:

  # Copy the bibliography file (note the change in the file name)
  configure_file (ice_bib.bib ice-bib.bib COPYONLY)

  # Find all the .tex files in userman and copy them to the build directory
  file (GLOB parts
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/userman
    userman/*.tex)
  foreach (filename ${parts})
    configure_file(userman/${filename} ${filename} COPYONLY)
  endforeach()

  # Find all the figures and copy them to the build directory
  foreach (format jpg png pdf)
    file (GLOB ${format}
      RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/userman/figs
      userman/figs/*.${format})
  endforeach()
  foreach (filename ${jpg} ${png} ${pdf})
    configure_file(userman/figs/${filename} figs/${filename} COPYONLY)
  endforeach()

  # Commands to build the manual:
  add_custom_command (OUTPUT manual.pdf
    COMMAND ${PDFLATEX_COMPILER} manual.tex
    COMMAND ${BIBTEX_COMPILER} manual
    COMMAND ${MAKEINDEX_COMPILER} -l manual
    COMMAND ${MAKEINDEX_COMPILER} -l -o manual.ond manual.odx
    COMMAND ${PDFLATEX_COMPILER} manual.tex
    COMMAND ${PDFLATEX_COMPILER} manual.tex
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  install (FILES "${CMAKE_CURRENT_BINARY_DIR}/manual.pdf"
    DESTINATION share/doc/pism)

  # A global target that gets built with "make all". Ensures that PISM
  # documentation is built for "make install"
  add_custom_target(docs ALL
    DEPENDS installation.pdf cheatsheet.pdf manual.pdf)
  

  # PISM Source Code Browser

  # Only build if both doxygen and dot are found
  if (DOXYGEN_EXECUTABLE AND DOXYGEN_DOT_EXECUTABLE)
    set (browser_helpers
      browser/config_doc.py
      browser/doxybib.bst
      browser/doxybib.py
      ice_bib.bib

      browser/bedrough.txt
      browser/bombproof_enth.txt
      browser/pages.txt
      browser/std_names.txt
      browser/vertchange.txt
      )

    # copy all the files needed to generate browser sources
    foreach (filename ${browser_helpers})
      configure_file(${filename} ${filename} COPYONLY)
    endforeach()

    # copy doxyfile, changing some things along the way:
    configure_file(browser/doxyfile browser/doxyfile)
    
    # Generate the list of references
    add_custom_command (OUTPUT browser/references.txt
      COMMAND ./doxybib.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/browser)

    # Generate the config file
    add_custom_command (OUTPUT browser/pism_config.nc
      COMMAND ncgen -o ${CMAKE_CURRENT_BINARY_DIR}/browser/pism_config.nc ${PROJECT_SOURCE_DIR}/src/pism_config.cdl)

    # Generate the config documentation
    add_custom_command (OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/browser/pism_config.txt
      COMMAND ./config_doc.py > pism_config.txt
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/browser/pism_config.nc
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/browser)

    # Generate the browser
    add_custom_command (OUTPUT browser/html/index.html
      COMMAND ${DOXYGEN_EXECUTABLE} doxyfile
      DEPENDS browser/pism_config.txt browser/references.txt browser/doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/browser)

    add_custom_target (browser ALL
      DEPENDS browser/html/index.html)

    add_custom_target (browser.tgz
      DEPENDS browser/html/index.html
      COMMAND tar -czf browser.tgz html
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/browser)

    install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/browser/html
      DESTINATION share/doc/pism)

  endif()
endif ()

