# Find latex and create doc targets:
find_package(LATEX)
find_package(Doxygen)

# PISM Source Code Browsers
if (LATEX_COMPILER AND DOXYGEN_EXECUTABLE AND DOXYGEN_DOT_EXECUTABLE)

  option(Pism_DOXYGEN_DOCUMENT_ALL "Extract all (even undocumented) entities." ON)

  if (Pism_DOXYGEN_DOCUMENT_ALL)
    set(Pism_DOXYGEN_EXTRACT_ALL "YES")
  else()
    set(Pism_DOXYGEN_EXTRACT_ALL "NO")
  endif()

  set (browser_helpers
    initialization-sequence.dot
    config_doc.py
    doxybib.bst
    doxybib.py
    bedrough.md
    bombproof_enth.md
    mainpage.md
    std_names.md
    vertchange.md
    howto.md
    petsc.md
    pism.sty
    )

  configure_file (../ice-bib.bib ice-bib.bib COPYONLY)

  # copy all the files needed to generate browser sources
  foreach (filename ${browser_helpers})
    configure_file(${filename} ${filename} COPYONLY)
  endforeach()

  # copy doxyfile, changing some things along the way:
  configure_file(doxyfile doxyfile)

  # Generate the list of references
  add_custom_command (OUTPUT references.md
    COMMAND ./doxybib.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # Generate the config file
  add_custom_command (OUTPUT pism_config.nc
    COMMAND ncgen -o ${CMAKE_CURRENT_BINARY_DIR}/pism_config.nc ${PROJECT_SOURCE_DIR}/src/pism_config.cdl)

  # Generate the config documentation
  add_custom_command (OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pism_config.md
    COMMAND ./config_doc.py > pism_config.md
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/pism_config.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # Generate the browser
  add_custom_target (browser
    COMMAND ${DOXYGEN_EXECUTABLE} doxyfile
    DEPENDS pism_config.md references.md doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # Create the compressed browser
  add_custom_target (browser.tgz
    COMMAND tar -czf browser.tgz html
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  add_dependencies(browser.tgz browser)

else ()
  if (NOT LATEX_COMPILER)
    message (STATUS "LaTeX was not found. Please install it if you want to re-build PISM documentation.")
  endif()

  if (NOT DOXYGEN_EXECUTABLE)
    message (STATUS "Doxygen was not found. Please install it if you want to re-build source code browsers.")
  endif()

  if (NOT DOXYGEN_DOT_EXECUTABLE)
    message (STATUS "'dot' was not found. Please install graphviz if you want to re-build source code browsers.")
  endif()
endif ()

