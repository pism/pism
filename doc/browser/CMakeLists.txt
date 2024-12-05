# Find doxygen
find_package(Doxygen)

# PISM Source Code Browsers
if (DOXYGEN_EXECUTABLE AND DOXYGEN_DOT_EXECUTABLE)
  find_package (Python3 COMPONENTS Interpreter)

  option(Pism_DOXYGEN_DOCUMENT_ALL "Extract all (even undocumented) entities." ON)

  if (Pism_DOXYGEN_DOCUMENT_ALL)
    set(Pism_DOXYGEN_EXTRACT_ALL "YES")
  else()
    set(Pism_DOXYGEN_EXTRACT_ALL "NO")
  endif()

  # copy doxyfile, changing some things along the way:
  configure_file(doxyfile doxyfile)

  configure_file(list_references.py list_references.py COPYONLY)
  configure_file (../ice-bib.bib ice-bib.bib COPYONLY)
  configure_file(mainpage.md mainpage.md COPYONLY)

  # Generate the list of references
  add_custom_command (OUTPUT references.md
    COMMAND ${Python3_EXECUTABLE} list_references.py ice-bib.bib > references.md
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # Generate the browser
  add_custom_target (browser
    COMMAND ${DOXYGEN_EXECUTABLE} doxyfile
    DEPENDS references.md doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

else ()
  add_custom_target(browser
    COMMAND ${CMAKE_COMMAND} -E echo "Please install Doxygen and Graphviz to build the source code browser.")
endif ()
