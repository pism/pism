# Find latex and create doc targets:
find_package(LATEX)
find_package(Doxygen)

# PISM Source Code Browsers
if (LATEX_COMPILER AND DOXYGEN_EXECUTABLE AND DOXYGEN_DOT_EXECUTABLE)
  option (Pism_BUILD_BROWSER "Build source code browsers if LaTeX, doxygen and dot are installed" OFF)

  set (browser_helpers
    initialization-sequence.dot
    config_doc.py
    doxybib.bst
    doxybib.py
    bedrough.txt
    bombproof_enth.txt
    pages.txt
    std_names.txt
    vertchange.txt
    howto.txt
    pism.sty
    )

  configure_file (../ice_bib.bib ice_bib.bib COPYONLY)

  # copy all the files needed to generate browser sources
  foreach (filename ${browser_helpers})
    configure_file(${filename} ${filename} COPYONLY)
  endforeach()

  # copy doxyfile, changing some things along the way:
  configure_file(doxyfile doxyfile)
  
  # Generate the list of references
  add_custom_command (OUTPUT references.txt
    COMMAND ./doxybib.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # Generate the config file
  add_custom_command (OUTPUT pism_config.nc
    COMMAND ncgen -o ${CMAKE_CURRENT_BINARY_DIR}/pism_config.nc ${PROJECT_SOURCE_DIR}/src/pism_config.cdl)

  # Generate the config documentation
  add_custom_command (OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pism_config.txt
    COMMAND ./config_doc.py > pism_config.txt
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/pism_config.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # Generate the browser
  add_custom_command (OUTPUT base/html/index.html
    COMMAND ${DOXYGEN_EXECUTABLE} doxyfile
    DEPENDS pism_config.txt references.txt doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # Create a compressed "base" browser
  add_custom_target (browser_base.tgz
    DEPENDS base/html/index.html
    COMMAND tar -czf browser_base.tgz html
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/base)

  ### Util browser
  configure_file(doxyfile_util doxyfile_util)
  configure_file (util_main.txt util_main.txt)

  # Generate the browser
  add_custom_command (OUTPUT util/html/index.html
    COMMAND ${DOXYGEN_EXECUTABLE} doxyfile_util
    DEPENDS doxyfile_util
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # Create a compressed Util browser
  add_custom_target (browser_util.tgz
    DEPENDS util/html/index.html
    COMMAND tar -czf browser_util.tgz html
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/util)

  # Add custom targets so that "make browser_base" works.
  add_custom_target (browser_base ALL
    DEPENDS base/html/index.html)
  add_custom_target (browser_util ALL
    DEPENDS util/html/index.html)

  # Install if Pism_BUILD_BROWSER is set
  if (Pism_BUILD_BROWSER)
    install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/base/html
      DESTINATION ${Pism_DOC_DIR}/base)
    install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/util/html
      DESTINATION ${Pism_DOC_DIR}/util)
  endif ()

endif ()

