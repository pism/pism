# Find latex and create doc targets:
find_package(LATEX)
find_package(Doxygen)

# PISM Source Code Browsers
if (LATEX_COMPILER AND DOXYGEN_EXECUTABLE AND DOXYGEN_DOT_EXECUTABLE)
  set (browser_helpers
    initialization-sequence.dot
    config_doc.py
    doxybib.bst
    doxybib.py
    bedrough.txt
    hydrology.txt
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

  # copy doxyfile_base, changing some things along the way:
  configure_file(doxyfile_base doxyfile_base)

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

  # Generate the "base" browser
  add_custom_target (browser_base
    COMMAND ${DOXYGEN_EXECUTABLE} doxyfile_base
    DEPENDS pism_config.txt references.txt doxyfile_base
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # Create the compressed "base" browser
  add_custom_target (browser_base.tgz
    COMMAND tar -czf browser_base.tgz html
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/base)

  add_dependencies(browser_base.tgz browser_base)

  ### Util browser
  configure_file(doxyfile_util doxyfile_util)
  configure_file (util_main.txt util_main.txt)

  # Generate the "util" browser
  add_custom_target (browser_util
    COMMAND ${DOXYGEN_EXECUTABLE} doxyfile_util
    DEPENDS doxyfile_util
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # Create the compressed "util" browser
  add_custom_target (browser_util.tgz
    COMMAND tar -czf browser_util.tgz html
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/util)

  add_dependencies(browser_util.tgz browser_util)

  add_custom_target (browsers)
  add_dependencies (browsers browser_base browser_util)

endif ()

